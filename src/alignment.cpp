#include "centrolign/alignment.hpp"

#include <unordered_map>
#include <unordered_set>
#include <stdexcept>

namespace centrolign {

using namespace std;

const uint64_t AlignedPair::gap = -1;

AlignedPair::AlignedPair(uint64_t node_id1, uint64_t node_id2) :
    node_id1(node_id1),
    node_id2(node_id2) {
    // nothing else to do
}

bool AlignedPair::operator==(const AlignedPair& other) const {
    return node_id1 == other.node_id1 && node_id2 == other.node_id2;
}
bool AlignedPair::operator<(const AlignedPair& other) const {
    return node_id1 < other.node_id1 || (node_id1 == other.node_id1 && node_id2 < other.node_id2);
}

void translate(Alignment& alignment,
               const vector<uint64_t>& back_translation1,
               const vector<uint64_t>& back_translation2) {
    
    for (size_t i = 0; i < alignment.size(); ++i) {
        auto& aln_pair = alignment[i];
        if (aln_pair.node_id1 != AlignedPair::gap) {
            aln_pair.node_id1 = back_translation1[aln_pair.node_id1];
        }
        if (aln_pair.node_id2 != AlignedPair::gap) {
            aln_pair.node_id2 = back_translation2[aln_pair.node_id2];
        }
    }
}

void swap_graphs(Alignment& alignment) {
    for (auto& aln_pair : alignment) {
        std::swap(aln_pair.node_id1, aln_pair.node_id2);
    }
}

std::string cigar(const Alignment& alignment) {
    // TODO: copypasta from explicit version
    std::stringstream strm;
    
    int curr_len = 0;
    char curr_op = '\0';
    for (const auto& aln_pair : alignment) {
        char op;
        if (aln_pair.node_id1 == AlignedPair::gap) {
            op = 'I';
        }
        else if (aln_pair.node_id2 == AlignedPair::gap) {
            op = 'D';
        }
        else {
            op = 'M';
        }
        
        if (op == curr_op) {
            ++curr_len;
        }
        else {
            if (curr_len != 0) {
                strm << curr_len << curr_op;
            }
            curr_len = 1;
            curr_op = op;
        }
    }
    
    if (curr_len != 0) {
        strm << curr_len << curr_op;
    }
    
    return strm.str();
}

std::string explicit_cigar(const Alignment& alignment,
                           const std::string& seq1, const std::string& seq2) {
    // TODO: copypasta from graph version
    std::stringstream strm;
    
    int curr_len = 0;
    char curr_op = '\0';
    for (const auto& aln_pair : alignment) {
        char op;
        if (aln_pair.node_id1 == AlignedPair::gap) {
            op = 'I';
        }
        else if (aln_pair.node_id2 == AlignedPair::gap) {
            op = 'D';
        }
        else if (seq1.at(aln_pair.node_id1) == seq2.at(aln_pair.node_id2)) {
            op = '=';
        }
        else {
            op = 'X';
        }
        
        if (op == curr_op) {
            ++curr_len;
        }
        else {
            if (curr_len != 0) {
                strm << curr_len << curr_op;
            }
            curr_len = 1;
            curr_op = op;
        }
    }
    
    if (curr_len != 0) {
        strm << curr_len << curr_op;
    }
    
    return strm.str();
}

// the maximum size mismatch we will infer from a double gap
// TODO: magic number
// TODO: could also align the underlying base sequence, but then why weren't they already aligned in the graph?
static const size_t max_mismatch_size = 4;

Alignment induced_pairwise_alignment(const BaseGraph& graph, uint64_t path_id1, uint64_t path_id2) {
    
    unordered_map<uint64_t, uint64_t> index_in_path1;
    
    const auto& path1 = graph.path(path_id1);
    const auto& path2 = graph.path(path_id2);
    
    for (uint64_t i = 0; i < path1.size(); ++i) {
        if (index_in_path1.count(path1[i])) {
            throw std::runtime_error("Cannot induce a colinear pairwise alignment from a sequence (" + graph.path_name(path_id1) + ") that follows cycles in the graph");
        }
        index_in_path1[path1[i]] = i;
    }
    
    Alignment alignment;
    
    // scan path 2
    uint64_t j = 0;
    for (uint64_t i = 0; i < path2.size(); ++i) {
        auto node_id = path2[i];
        auto it = index_in_path1.find(node_id);
        if (it == index_in_path1.end()) {
            // no aligned base, must be a gap
            alignment.emplace_back(AlignedPair::gap, i);
        }
        else {
            // clear out earlier gaps in path1
            while (j < it->second) {
                alignment.emplace_back(j++, AlignedPair::gap);
            }
            // record the aligned pair
            alignment.emplace_back(j++, i);
        }
    }
    
    // finish off path 1, if necessary
    while (j < path1.size()) {
        alignment.emplace_back(j++, AlignedPair::gap);
    }
    
    // consolidate mismatches equal length gap runs runs as a mismatch
    size_t removed = 0;
    for (size_t i = 0; i < alignment.size();) {
        auto& aln_pair = alignment[i];
        if (aln_pair.node_id1 != AlignedPair::gap && aln_pair.node_id2 != AlignedPair::gap) {
            alignment[i - removed] = aln_pair;
            ++i;
        }
        else {
            // find the extent of the gap run
            size_t j = i;
            size_t gaps1 = 0;
            size_t gaps2 = 0;
            while (j < alignment.size() && (alignment[j].node_id1 == AlignedPair::gap ||
                                            alignment[j].node_id2 == AlignedPair::gap)) {
                gaps1 += (alignment[j].node_id1 == AlignedPair::gap);
                gaps2 += (alignment[j].node_id2 == AlignedPair::gap);
                ++j;
            }
            
            
            // find the final index for each sequence in this run
            uint64_t last1 = -1;
            uint64_t last2 = -1;
            if (i != 0) {
                // this node is guaranteed to be a match or it would have
                // ended up in this run
                auto& last_pair = alignment[i - removed - 1];
                last1 = last_pair.node_id1;
                last2 = last_pair.node_id2;
            }
            
            if (gaps1 == gaps2 && gaps1 <= max_mismatch_size) {
                // they have the same length, so we'll call this a mismatch
                for (uint64_t n = 0; n < gaps1; ++n) {
                    alignment[i - removed + n] = AlignedPair(last1 + n + 1, last2 + n + 1);
                }
                removed += gaps1;
            }
            else {
                // since they're not the same length, we interpret this as a separate insertion
                // and deletion
                
                // consolidate into a single gap for each sequence
                for (uint64_t n = 0; n < gaps2; ++n) {
                    alignment[i - removed + n] = AlignedPair(last1 + n + 1, AlignedPair::gap);
                }
                for (uint64_t n = 0; n < gaps1; ++n) {
                    alignment[i - removed + n + gaps2] = AlignedPair(AlignedPair::gap, last2 + n + 1);
                }
            }
            i = j;
        }
    }
    alignment.resize(alignment.size() - removed);
    
    return alignment;
}

void induced_cyclic_pairwise_alignment_internal(const std::vector<uint64_t>& path1,
                                                const std::vector<uint64_t>& path2,
                                                const std::pair<size_t, size_t>& coord_begin,
                                                const std::pair<size_t, size_t>& coord_end,
                                                std::vector<Alignment>& alignments) {
    
    static const bool debug = false;
    if (debug) {
        std::cerr << "entering induced pairwise alignment algorithm between coordinates " << coord_begin.first << ", " << coord_begin.second << " and " << coord_end.first << ", " << coord_end.second << '\n';
    }
        
    SliceView<std::vector<uint64_t>> subpath1(path1, coord_begin.first, coord_end.first);
    SliceView<std::vector<uint64_t>> subpath2(path2, coord_begin.second, coord_end.second);
    
    Alignment aln = long_common_subsequence_nonrepeating(subpath1, subpath2);
    
    if (aln.empty()) {
        // base case, no more match possible
        if (debug) {
            std::cerr << "no common subsequences left, exiting\n";
        }
        return;
    }
    
    // adjust the indexes
    for (auto& aln_pair : aln) {
        if (aln_pair.node_id1 != AlignedPair::gap) {
            aln_pair.node_id1 += coord_begin.first;
        }
        if (aln_pair.node_id2 != AlignedPair::gap) {
            aln_pair.node_id2 += coord_begin.second;
        }
    }
    
    // convert short double gaps into mismatches
    size_t removed = 0;
    for (size_t i = 0; i < aln.size();) {
        if (aln[i].node_id1 != AlignedPair::gap && aln[i].node_id2 != AlignedPair::gap) {
            // an aligned pair of node IDs
            aln[i - removed] = aln[i];
            ++i;
        }
        else {
            // walk out the next run of gaps
            size_t j = i, gap1 = 0, gap2 = 0;
            while (j < aln.size() && (aln[j].node_id1 == AlignedPair::gap || aln[j].node_id2 == AlignedPair::gap)) {
                if (aln[j].node_id1 == AlignedPair::gap) {
                    ++gap1;
                }
                else {
                    ++gap2;
                }
                ++j;
            }
            if (gap1 == gap2 && gap1 <= max_mismatch_size) {
                // we infer a mismatch here
                size_t g1 = i - removed;
                size_t g2 = g1;
                for (size_t k = i; k < j; ++k) {
                    if (aln[k].node_id1 == AlignedPair::gap) {
                        aln[g2++].node_id2 = aln[k].node_id2;
                    }
                    else {
                        aln[g1++].node_id1 = aln[k].node_id1;
                    }
                }
                
                removed += gap1;
            }
            else {
                // these sequences are unaligned
                for (size_t k = i; k < j; ++k) {
                    aln[k - removed] = aln[k];
                }
            }
            i = j;
        }
    }
    
    aln.resize(aln.size() - removed);
    
    if (debug) {
        std::cerr << "got LCS alignment:\n";
        for (auto ap : aln) {
            std::cerr << (int64_t) ap.node_id1 << '\t' << (int64_t) ap.node_id2 << '\n';
        }
    }
    
    std::pair<size_t, size_t> aln_coord_begin(aln.front().node_id1, aln.front().node_id2);
    std::pair<size_t, size_t> aln_coord_end(aln.back().node_id1 + 1, aln.back().node_id2 + 1);
    
    alignments.emplace_back(std::move(aln));
    
    if (aln_coord_begin.first != coord_begin.first && aln_coord_begin.second != coord_begin.second)  {
        // recursive call into the front of the sequence
        induced_cyclic_pairwise_alignment_internal(path1, path2, coord_begin, aln_coord_begin, alignments);
    }
    if (aln_coord_end.first != coord_end.first && aln_coord_end.second != coord_end.second) {
        // recursive call into the back of the sequence
        induced_cyclic_pairwise_alignment_internal(path1, path2, aln_coord_end, coord_end, alignments);
    }
}

std::vector<std::pair<size_t, size_t>> maximum_noncyclic_extension(const std::vector<uint64_t>& path,
                                                                   const std::vector<std::pair<size_t, size_t>>& covered_intervals) {
    
    static const bool debug = false;
    
    if (debug) {
        std::cerr << "extending intervals:\n";
        for (auto interval : covered_intervals) {
            std::cerr << '\t' << interval.first << '\t' << interval.second << '\n';
        }
        std::cerr << "in path:\n";
        for (auto n : path) {
            std::cerr << n << ' ';
        }
        std::cerr << '\n';
    }
    
    std::vector<std::vector<std::pair<size_t, size_t>>> maximal_interval_extensions(covered_intervals.size());
    
    // get the left-to-right order
    auto lex_order = range_vector(covered_intervals.size());
    std::sort(lex_order.begin(), lex_order.end(), [&](size_t i, size_t j) {
        return covered_intervals[i] < covered_intervals[j];
    });
    
    for (size_t i = 0; i < lex_order.size(); ++i) {
        
        auto& extensions = maximal_interval_extensions[lex_order[i]];
        
        const auto& interval = covered_intervals[lex_order[i]];
        if (debug) {
            std::cerr << "finding extensions for interval " << interval.first << ", " << interval.second << '\n';
        }
        
        // the interval we will try to cover
        size_t left_lim = (i == 0 ? 0 : covered_intervals[lex_order[i - 1]].second);
        size_t right_lim = (i + 1 == lex_order.size() ? path.size() : covered_intervals[lex_order[i + 1]].first);
        
        // the nodes already inside the interval
        std::unordered_set<uint64_t> interval_nodes(path.begin() + interval.first,
                                                    path.begin() + interval.second);
        
        // record positions of nodes in the left flank
        std::unordered_map<uint64_t, size_t> left_flank_positions;
        for (size_t j = interval.first; j > left_lim; --j) {
            uint64_t node_id = path[j - 1];
            if (left_flank_positions.count(node_id) || interval_nodes.count(node_id)) {
                // we can't extend any further to the left without including a cycle
                break;
            }
            left_flank_positions[node_id] = j - 1;
        }
        
        if (debug) {
            std::cerr << "left flank positions:\n";
            for (auto r : left_flank_positions) {
                std::cerr << '\t' << r.first << '\t' << r.second << '\n';
            }
        }
        
        std::pair<size_t, size_t> current_interval(interval.first - left_flank_positions.size(), interval.second);
        if (debug) {
            std::cerr << "left boundary starts at " << current_interval.first << '\n';
        }
        for (size_t j = interval.second; j < right_lim; ++j) {
            uint64_t node_id = path[j];
            if (interval_nodes.count(node_id)) {
                // we can't extend any further without including a cycle
                if (debug) {
                    std::cerr << "breaking at requisite loop in position " << j << '\n';
                }
                break;
            }
            auto it = left_flank_positions.find(node_id);
            if (it != left_flank_positions.end() && it->second >= current_interval.first) {
                // we need to pull in the left end of the interval to include the current node
                
                // yield a maximal interval
                extensions.push_back(current_interval);
                if (debug) {
                    std::cerr << "found optional loop in position " << j << ", yielding interval " << current_interval.first << ", " << current_interval.second << '\n';
                }
                
                current_interval.first = it->second + 1;
                if (debug) {
                    std::cerr << "new left boundary is " << current_interval.first << '\n';
                }
            }
            // add the current node to the interval
            ++current_interval.second;
            if (debug) {
                std::cerr << "extend right boundary to " << current_interval.second << '\n';
            }
            
            interval_nodes.insert(node_id);
        }
        // yield the final interval
        extensions.push_back(current_interval);
    }
    
    if (debug) {
        std::cerr << "computed extension sets:\n";
        for (size_t i = 0; i < maximal_interval_extensions.size(); ++i) {
            std::cerr << covered_intervals[i].first << ", " << covered_intervals[i].second << '\n';
            for (auto extension : maximal_interval_extensions[i]) {
                std::cerr << '\t' << extension.first << ", " << extension.second << '\n';
            }
        }
    }
    
    // records of (gaps closed to left, positions covered to left, backpointer)
    std::vector<std::vector<std::tuple<size_t, size_t, size_t>>> dp(covered_intervals.size());
    
    // main dynamic programming iterations
    for (size_t i = 0; i < lex_order.size(); ++i) {
        
        const auto& extensions = maximal_interval_extensions[lex_order[i]];
        const auto& interval = covered_intervals[lex_order[i]];
        auto& dp_col = dp[lex_order[i]];
        
        dp_col.resize(extensions.size(), std::tuple<size_t, size_t, size_t>(0, 0, -1));
        
        if (i == 0) {
            // first interval, consider flank to the start of the path
            for (size_t j = 0; j < dp_col.size(); ++j) {
                auto& dp_entry = dp_col[j];
                const auto& extension = extensions[j];
                std::get<0>(dp_entry) = (interval.first != 0 && extension.first == 0) ? 1 : 0;
                std::get<1>(dp_entry) = interval.first - extension.first;
            }
        }
        else {
            // middle interval, consider flank up to the previous extensions
            
            const auto& prev_interval = covered_intervals[lex_order[i - 1]];
            const auto& prev_extensions = maximal_interval_extensions[lex_order[i - 1]];
            const auto& prev_dp_col = dp[lex_order[i - 1]];
            
            assert(prev_interval.second <= interval.first);
            
            for (size_t j = 0; j < dp_col.size(); ++j) {
                auto& dp_entry = dp_col[j];
                const auto& extension = extensions[j];
                // try all extensions of previous interval as a predecessors
                for (size_t k = 0; k < prev_dp_col.size(); ++k) {
                    const auto& prev_dp_entry = prev_dp_col[k];
                    const auto& prev_extension = prev_extensions[k];
                    
                    size_t gaps_covered = std::get<0>(prev_dp_entry);
                    size_t bases_covered = std::get<1>(prev_dp_entry);
                    if (prev_interval.second != interval.first && prev_extension.second >= extension.first) {
                        // this maximal extension combo covers the gap
                        gaps_covered += 1;
                        bases_covered += (interval.first - prev_interval.second);
                    }
                    else {
                        // it doesn't span the whole gap
                        bases_covered += (prev_extension.second - prev_interval.second +
                                          interval.first - extension.first);
                    }
                    
                    if (gaps_covered > std::get<0>(dp_entry) ||
                        (gaps_covered == std::get<0>(dp_entry) && bases_covered >= std::get<1>(dp_entry))) {
                        // new opt for this entry
                        dp_entry = std::make_tuple(gaps_covered, bases_covered, k);
                    }
                }
            }
        }
    }
    
    if (debug) {
        std::cerr << "final DP state:\n";
        for (size_t i = 0; i < dp.size(); ++i) {
            std::cerr << i << ": " << covered_intervals[i].first << ", " << covered_intervals[i].second << '\n';
            for (size_t j = 0; j < dp[i].size(); ++j) {
                std::cerr << '\t' << j << " (" << maximal_interval_extensions[i][j].first << ", " << maximal_interval_extensions[i][j].second << "): " << std::get<0>(dp[i][j]) << ", "  << std::get<1>(dp[i][j]) << ", "  << std::get<2>(dp[i][j]) << '\n';
            }
        }
    }
    
    // find the traceback opt
    size_t opt_idx = -1;
    size_t opt_gaps_covered = 0;
    size_t opt_bases_covered = 0;
    if (!lex_order.empty()) {
        const auto& final_dp_col = dp[lex_order.back()];
        const auto& final_interval = covered_intervals[lex_order.back()];
        const auto& final_extensions = maximal_interval_extensions[lex_order.back()];
        for (size_t j = 0; j < final_dp_col.size(); ++j) {
            
            // add the final gap to the right of the last interval
            const auto& dp_entry = final_dp_col[j];
            const auto& extension = final_extensions[j];
            size_t gaps_covered = std::get<0>(dp_entry);
            size_t bases_covered = std::get<1>(dp_entry);
            if (final_interval.second != path.size() && extension.second == path.size()) {
                gaps_covered += 1;
            }
            bases_covered += (extension.second - final_interval.second);
            
            if (gaps_covered > opt_gaps_covered ||
                (gaps_covered == opt_gaps_covered && bases_covered >= opt_bases_covered)) {
                opt_idx = j;
                opt_gaps_covered = gaps_covered;
                opt_bases_covered = bases_covered;
            }
        }
        
        if (debug) {
            std::cerr << "choose final column (" << lex_order.back() << ") opt at " << opt_idx << " with " << opt_gaps_covered << " gaps closed and " << opt_bases_covered << " bases covered\n";
        }
    }
    
    // follow the backpointers to traceback best set of extensions
    size_t tb_row = opt_idx;
    std::vector<std::pair<size_t, size_t>> maximum_coverage_extensions(covered_intervals.size());
    for (int64_t i = lex_order.size() - 1; i >= 0; --i) {
        maximum_coverage_extensions[lex_order[i]] = maximal_interval_extensions[lex_order[i]][tb_row];
        tb_row = std::get<2>(dp[lex_order[i]][tb_row]);
        if (i + 1 != lex_order.size()) {
            // trim back this extension so it doesn't overlap the next one
            maximum_coverage_extensions[lex_order[i]].second = std::min(maximum_coverage_extensions[lex_order[i]].second,
                                                                        maximum_coverage_extensions[lex_order[i + 1]].first);
        }
    }
    
    if (debug) {
        std::cerr << "final extensions:\n";
        for (auto interval : maximum_coverage_extensions) {
            std::cerr << '\t' << interval.first << '\t' << interval.second << '\n';
        }
    }
    
    return maximum_coverage_extensions;
}

std::vector<Alignment> induced_cyclic_pairwise_alignment(const BaseGraph& graph, uint64_t path_id1, uint64_t path_id2) {
    
    static const bool debug = false;
    
    if (debug) {
        std::cerr << "start induced cyclic pairwise alignment for paths " << path_id1 << " (" << graph.path_name(path_id1) << ") and " << path_id2 << " (" << graph.path_name(path_id2) << ") of " << graph.path_size() << "\n";
    }
    
    const auto& path1 = graph.path(path_id1);
    const auto& path2 = graph.path(path_id2);
    
    // get the non-overlapping matching intervals
    std::vector<Alignment> alignments;
    induced_cyclic_pairwise_alignment_internal(path1, path2,
                                               std::pair<size_t, size_t>(0, 0),
                                               std::pair<size_t, size_t>(path1.size(), path2.size()),
                                               alignments);
    
    // pull out the covered intervals on each of the two paths
    std::vector<std::pair<size_t, size_t>> covered_intervals1, covered_intervals2;
    for (const auto& aln : alignments) {
        covered_intervals1.emplace_back(aln.front().node_id1, aln.back().node_id1 + 1);
        covered_intervals2.emplace_back(aln.front().node_id2, aln.back().node_id2 + 1);
    }
    
    if (debug) {
        std::cerr << "covered intervals on path 1:\n";
        for (const auto interval : covered_intervals1) {
            std::cerr << '\t' << interval.first << '\t' << interval.second << '\n';
        }
        std::cerr << "covered intervals on path 2:\n";
        for (const auto interval : covered_intervals2) {
            std::cerr << '\t' << interval.first << '\t' << interval.second << '\n';
        }
    }
    
    // try to to merge blocks together if they are adjacent and don't create a cycle
    {
        auto order1 = range_vector(covered_intervals1.size());
        auto order2 = range_vector(covered_intervals2.size());
        std::sort(order1.begin(), order1.end(), [&](size_t i, size_t j) {
            return covered_intervals1[i] < covered_intervals1[j];
        });
        std::sort(order2.begin(), order2.end(), [&](size_t i, size_t j) {
            return covered_intervals2[i] < covered_intervals2[j];
        });
        auto index2 = invert(order2);
        
        std::vector<bool> keep(covered_intervals1.size(), true);
        size_t merge_run = 0;
        std::unordered_set<uint64_t> node_set1, node_set2;
        for (size_t i = 1; i < order1.size(); ++i) {
            bool did_merge = false;
            if (index2[order1[i]] == index2[order1[i - 1]] + 1) {
                // these two blocks are adjacent in both paths
                
                if (node_set1.empty()) {
                    // we can't recycle a node set, make a new one for the predecessor
                    for (size_t j = covered_intervals1[order1[i - 1 - merge_run]].first; j < covered_intervals1[order1[i - 1 - merge_run]].second; ++j) {
                        node_set1.insert(path1[j]);
                    }
                    for (size_t j = covered_intervals2[order1[i - 1 - merge_run]].first; j < covered_intervals2[order1[i - 1 - merge_run]].second; ++j) {
                        node_set2.insert(path2[j]);
                    }
                }
                
                // check that no nodes repeat in between or in the next interval
                bool compatible = true;
                for (size_t j = covered_intervals1[order1[i - 1 - merge_run]].second; j < covered_intervals1[order1[i]].second && compatible; ++j) {
                    compatible = node_set1.insert(path1[j]).second;
                }
                for (size_t j = covered_intervals2[order1[i - 1 - merge_run]].second; j < covered_intervals2[order1[i]].second && compatible; ++j) {
                    compatible = node_set2.insert(path2[j]).second;
                }
                
                if (compatible) {
                    // these blocks are mergeable
                    
                    if (debug) {
                        std::cerr << "intervals " << order1[i - 1 - merge_run] << " and " << order1[i] << " are mergeable\n";
                    }
                    
                    auto& aln = alignments[order1[i - 1 - merge_run]];
                    
                    // add the intervening region as a double deletion
                    for (size_t j = covered_intervals1[order1[i - 1 - merge_run]].second; j < covered_intervals1[order1[i]].first; ++j) {
                        aln.emplace_back(j, AlignedPair::gap);
                    }
                    for (size_t j = covered_intervals2[order1[i - 1 - merge_run]].second; j < covered_intervals2[order1[i]].first; ++j) {
                        aln.emplace_back(AlignedPair::gap, j);
                    }
                    // copy the second block
                    for (const auto& aln_pair : alignments[order1[i]]) {
                        aln.emplace_back(aln_pair);
                    }
                    alignments[order1[i]].clear();
                    keep[order1[i]] = false;
                    
                    covered_intervals1[order1[i - 1 - merge_run]].second = covered_intervals1[order1[i]].second;
                    covered_intervals2[order1[i - 1 - merge_run]].second = covered_intervals2[order1[i]].second;
                    
                    did_merge = true;
                }
            }
            
            if (did_merge) {
                ++merge_run;
            }
            else {
                // we don't want to recycle a merged node set, so dump it
                node_set1.clear();
                node_set2.clear();
                // reset the run of merges
                merge_run = 0;
            }
            
        }
        
        // move items up into the prefix of the vector
        size_t removed = 0;
        for (size_t i = 0; i < alignments.size(); ++i) {
            if (!keep[i]) {
                ++removed;
            }
            else if (removed != 0) {
                covered_intervals1[i - removed] = covered_intervals1[i];
                covered_intervals2[i - removed] = covered_intervals2[i];
                alignments[i - removed] = std::move(alignments[i]);
            }
        }
        if (removed != 0) {
            // get rid of the emptied positions
            covered_intervals1.resize(covered_intervals1.size() - removed);
            covered_intervals2.resize(covered_intervals2.size() - removed);
            alignments.resize(alignments.size() - removed);
        }
    }
    
    if (debug) {
        std::cerr << "intervals after merging on path 1:\n";
        for (const auto interval : covered_intervals1) {
            std::cerr << '\t' << interval.first << '\t' << interval.second << '\n';
        }
        std::cerr << "intervals after merging on path 2:\n";
        for (const auto interval : covered_intervals2) {
            std::cerr << '\t' << interval.first << '\t' << interval.second << '\n';
        }
    }
    
    // extend them as much as possible to eliminate gaps
    auto extended_intervals1 = maximum_noncyclic_extension(path1, covered_intervals1);
    auto extended_intervals2 = maximum_noncyclic_extension(path2, covered_intervals2);
    
    if (debug) {
        std::cerr << "extended intervals on path 1:\n";
        for (const auto interval : extended_intervals1) {
            std::cerr << '\t' << interval.first << '\t' << interval.second << '\n';
        }
        std::cerr << "extended intervals on path 2:\n";
        for (const auto interval : extended_intervals2) {
            std::cerr << '\t' << interval.first << '\t' << interval.second << '\n';
        }
    }
    
    for (size_t i = 0; i < alignments.size(); ++i) {
        
        auto& aln = alignments[i];
        
        // initially add the left inserts on the right
        size_t num_added_left = 0;
        for (size_t j = extended_intervals1[i].first; j < covered_intervals1[i].first; ++j) {
            aln.emplace_back(j, AlignedPair::gap);
            ++num_added_left;
        }
        for (size_t j = extended_intervals2[i].first; j < covered_intervals2[i].first; ++j) {
            aln.emplace_back(AlignedPair::gap, j);
            ++num_added_left;
        }
        // shift them onto the left side
        std::rotate(aln.rbegin(), aln.rbegin() + num_added_left, aln.rend());
        
        // add the right inserts
        for (size_t j = covered_intervals1[i].second; j < extended_intervals1[i].second; ++j) {
            aln.emplace_back(j, AlignedPair::gap);
        }
        for (size_t j = covered_intervals2[i].second; j < extended_intervals2[i].second; ++j) {
            aln.emplace_back(AlignedPair::gap, j);
        }
    }
    
    if (debug) {
        std::cerr << "extended alignments:\n";
        for (size_t i = 0; i < alignments.size(); ++i) {
            std::unordered_map<uint64_t, size_t> node_set1, node_set2;
            std::cerr << "alignment " << i << '\n';
            for (size_t j = 0; j < alignments[i].size(); ++j) {
                auto ap = alignments[i][j];
                if (ap.node_id1 != AlignedPair::gap) {
                    if (node_set1.count(path1[ap.node_id1])) {
                        std::cerr << "error: repeat node in path 1 at positions " << alignments[i][node_set1[path1[ap.node_id1]]].node_id1 << " and " << ap.node_id1 << '\n';
                        exit(1);
                    }
                    node_set1[path1[ap.node_id1]] = j;
                }
                if (ap.node_id2 != AlignedPair::gap) {
                    if (node_set2.count(path2[ap.node_id2])) {
                        std::cerr << "error: repeat node in path 2 at positions " << alignments[i][node_set2[path2[ap.node_id2]]].node_id2 << " and " << ap.node_id2 << '\n';
                        exit(1);
                    }
                    node_set2[path2[ap.node_id2]] = j;
                }
                std::cerr << '\t' << (int64_t) ap.node_id1 << '\t' << (int64_t) ap.node_id2 << '\n';
            }
        }
    }
        
    auto order = range_vector(alignments.size());
    
    // make dangling insertions for path 1
    std::sort(order.begin(), order.end(), [&](size_t i, size_t j) {
        return extended_intervals1[i].first < extended_intervals1[j].first;
    });
    for (size_t i = 0; i <= order.size(); ++i) {
        
        size_t l = (i == 0 ? 0 : extended_intervals1[order[i - 1]].second);
        size_t r = (i == order.size() ? path1.size() : extended_intervals1[order[i]].first);
        
        if (l != r) {
            if (debug) {
                std::cerr << "adding dangling inserts for path 1 interval " << l << ", " << r << '\n';
            }
            std::unordered_set<uint64_t> nodes_seen;
            alignments.emplace_back();
            for (size_t j = l; j < r; ++j) {
                if (nodes_seen.count(path1[j])) {
                    // we have to break this into another block to avoid having a cycle
                    alignments.emplace_back();
                    nodes_seen.clear();
                }
                // extend the alignment by this node
                alignments.back().emplace_back(j, AlignedPair::gap);
                nodes_seen.insert(path1[j]);
            }
        }
    }
    
    // make dangling insertions for path 2
    std::sort(order.begin(), order.end(), [&](size_t i, size_t j) {
        return extended_intervals2[i].first < extended_intervals2[j].first;
    });
    for (size_t i = 0; i <= order.size(); ++i) {
        size_t l = (i == 0 ? 0 : extended_intervals2[order[i - 1]].second);
        size_t r = (i == order.size() ? path2.size() : extended_intervals2[order[i]].first);
        if (l != r) {
            if (debug) {
                std::cerr << "adding dangling inserts for " << i << "-th of " << order.size() << " path 2 interval " << l << ", " << r << '\n';
            }
            std::unordered_set<uint64_t> nodes_seen;
            alignments.emplace_back();
            for (size_t j = l; j < r; ++j) {
                if (nodes_seen.count(path2[j])) {
                    // we have to break this into another block to avoid having a cycle
                    alignments.emplace_back();
                    nodes_seen.clear();
                }
                // extend the alignment by this node
                alignments.back().emplace_back(AlignedPair::gap, j);
                nodes_seen.insert(path2[j]);
            }
        }
    }
    
    if (debug) {
        std::cerr << "completed alignment set has " << alignments.size() << " blocks\n";
    }
    
    return alignments;
}

void output_maf(std::ostream& out, const std::vector<Alignment>& blocks,
                const BaseGraph& graph, uint64_t path_id1, uint64_t path_id2) {
    
    const auto& path1 = graph.path(path_id1);
    const auto& path2 = graph.path(path_id2);
    
    // header
    out << "track name=" << graph.path_name(path_id1) << "_vs_" << graph.path_name(path_id2) << "_induced\n";
    out << "##maf version=1\n";
    for (const auto& block : blocks) {
        
        // measure the sequences
        size_t start1 = -1, start2 = -1, size1 = 0, size2 = 0;
        for (const auto& aln_pair : block) {
            if (aln_pair.node_id1 != AlignedPair::gap) {
                if (start1 == -1) {
                    start1 = aln_pair.node_id1;
                }
                ++size1;
            }
            if (aln_pair.node_id2 != AlignedPair::gap) {
                if (start2 == -1) {
                    start2 = aln_pair.node_id2;
                }
                ++size2;
            }
        }
        out << "\na\n";
        out << "s\t" << graph.path_name(path_id1) << '\t' << (start1 == -1 ? path1.size() : start1) << '\t' << size1 << "\t+\t";
        for (const auto& aln_pair : block) {
            if (aln_pair.node_id1 != AlignedPair::gap) {
                char base = graph.label(path1[aln_pair.node_id1]);
                if (base <= 4) {
                    base = decode_base(base);
                }
                out << base;
            }
            else {
                out << '-';
            }
        }
        out << '\n';
        out << "s\t" << graph.path_name(path_id2) << '\t' << (start2 == -1 ? path2.size() : start2) << '\t' << size2 << "\t+\t";
        for (const auto& aln_pair : block) {
            if (aln_pair.node_id2 != AlignedPair::gap) {
                char base = graph.label(path2[aln_pair.node_id2]);
                if (base <= 4) {
                    base = decode_base(base);
                }
                out << base;
            }
            else {
                out << '-';
            }
        }
        out << '\n';
    }
    
    
}



}
