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
                                                std::vector<Alignment>& alignments,
                                                std::vector<std::pair<std::pair<size_t, size_t>, std::pair<size_t, size_t>>>& intervals) {
    
    static const size_t max_mismatch_size = 4;
    
    
    SliceView<std::vector<uint64_t>> subpath1(path1, coord_begin.first, coord_end.first);
    SliceView<std::vector<uint64_t>> subpath2(path2, coord_begin.second, coord_end.second);
    
    Alignment aln;
    std::pair<size_t, size_t> aln_coord_begin, aln_coord_end;
    std::tie(aln, aln_coord_begin, aln_coord_end) = longest_common_subsequence_nonrepeating(subpath1, subpath2);
    
    if (aln.empty()) {
        // base case, no more match possible
        return;
    }
    
    // convert these to past-the-last
    ++aln_coord_end.first;
    ++aln_coord_end.second;
    
    // adjust the coordinates from relative to absolute
    aln_coord_begin.first += coord_begin.first;
    aln_coord_end.first += coord_begin.first;
    aln_coord_begin.second += coord_begin.second;
    aln_coord_end.second += coord_begin.second;
    
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
    
    alignments.emplace_back(std::move(aln));
    intervals.emplace_back(aln_coord_begin, aln_coord_begin);
    
    if (aln_coord_begin.first != coord_begin.first && aln_coord_begin.second != coord_begin.second)  {
        // recursive call into the front of the sequence
        induced_cyclic_pairwise_alignment_internal(path1, path2, coord_begin, aln_coord_begin,
                                                   alignments, intervals);
    }
    if (aln_coord_end.first != coord_end.first && aln_coord_end.second != coord_end.second) {
        // recursive call into the back of the sequence
        induced_cyclic_pairwise_alignment_internal(path1, path2, aln_coord_end, coord_end,
                                                   alignments, intervals);
    }
}

std::vector<std::pair<size_t, size_t>> maximum_noncyclic_extension(const std::vector<uint64_t>& path,
                                                                   const std::vector<std::pair<size_t, size_t>>& covered_intervals) {
    
    std::vector<std::vector<std::pair<size_t, size_t>>> maximal_interval_extensions(covered_intervals.size());
    
    // get the left-to-right order
    auto lex_order = range_vector(covered_intervals.size());
    std::sort(lex_order.begin(), lex_order.end(), [&](size_t i, size_t j) {
        return covered_intervals[i] < covered_intervals[j];
    });
    
    for (size_t i = 0; i < lex_order.size(); ++i) {
        
        auto& extensions = maximal_interval_extensions[lex_order[i]];
        
        const auto& interval = covered_intervals[lex_order[i]];
        
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
        
        std::pair<size_t, size_t> current_interval(interval.first - left_flank_positions.size(), interval.second);
        for (size_t j = interval.second; j < right_lim; ++j) {
            uint64_t node_id = path[j];
            if (interval_nodes.count(node_id)) {
                // we can't extend any further without including a cycle
                break;
            }
            auto it = left_flank_positions.find(node_id);
            if (it != left_flank_positions.end() && it->second >= current_interval.first) {
                // we need to pull in the left end of the interval to include the current node
                
                // yield a maximal interval
                extensions.push_back(current_interval);
                
                current_interval.first = it->second + 1;
            }
            // add the current node to the interval
            ++current_interval.second;
            
            interval_nodes.insert(node_id);
        }
        // yield the final interval
        extensions.push_back(current_interval);
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
    
    // find the traceback opt
    size_t opt_idx = -1;
    size_t opt_gaps_covered = 0;
    size_t opt_bases_covered = 0;
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
    
    return maximum_coverage_extensions;
}

std::vector<Alignment> induced_cyclic_pairwise_alignment(const BaseGraph& graph, uint64_t path_id1, uint64_t path_id2) {
    
    
    const auto& path1 = graph.path(path_id1);
    const auto& path2 = graph.path(path_id2);
    
    // get the non-overlapping matching intervals
    std::vector<Alignment> alignments;
    std::vector<std::pair<std::pair<size_t, size_t>, std::pair<size_t, size_t>>> intervals;
    induced_cyclic_pairwise_alignment_internal(path1, path2,
                                               std::pair<size_t, size_t>(0, 0),
                                               std::pair<size_t, size_t>(path1.size(), path2.size()),
                                               alignments, intervals);
    
    // pull out the covered intervals on each of the two paths
    std::vector<std::pair<size_t, size_t>> covered_intervals1, covered_intervals2;
    covered_intervals1.reserve(intervals.size());
    covered_intervals2.reserve(intervals.size());
    for (const auto& range : intervals) {
        covered_intervals1.emplace_back(range.first.first, range.second.first);
        covered_intervals2.emplace_back(range.first.second, range.second.second);
    }
    
    // extend them as much as possible to eliminate gaps
    auto extended_intervals1 = maximum_noncyclic_extension(path1, covered_intervals1);
    auto extended_intervals2 = maximum_noncyclic_extension(path2, covered_intervals2);
    
    for (size_t i = 0; i < alignments.size(); ++i) {
        
        auto& aln_intervals = intervals[i];
        auto& aln = alignments[i];
        
        // initially add the left inserts on the right
        size_t num_added_left = 0;
        for (size_t j = extended_intervals1[i].first; j < aln_intervals.first.first; ++j) {
            aln.emplace_back(path1[j], AlignedPair::gap);
            ++num_added_left;
        }
        for (size_t j = extended_intervals2[i].first; j < aln_intervals.first.second; ++j) {
            aln.emplace_back(AlignedPair::gap, path2[j]);
            ++num_added_left;
        }
        // shift them onto the left side
        std::rotate(aln.rbegin(), aln.rbegin() + num_added_left, aln.rend());
        
        // add the right inserts
        for (size_t j = aln_intervals.second.first; j < extended_intervals1[i].second; ++j) {
            aln.emplace_back(path1[j], AlignedPair::gap);
        }
        for (size_t j = aln_intervals.second.second; j < extended_intervals2[i].second; ++j) {
            aln.emplace_back(AlignedPair::gap, path2[j]);
        }
            
        // update the alignment interval
        aln_intervals.first.first = extended_intervals1[i].first;
        aln_intervals.second.first = extended_intervals1[i].second;
        aln_intervals.first.second = extended_intervals2[i].first;
        aln_intervals.second.second = extended_intervals2[i].second;
    }
    
    // TODO: add any remaining gaps as dangling blocks
    
    // TODO: covert to sequence indexes
    
    
    return alignments;
}


}
