#include "centrolign/stitcher.hpp"

#include <chrono>
#include <stdexcept>

namespace centrolign {

using namespace std;

const bool Stitcher::debug = false;
const bool Stitcher::instrument = false;

Stitcher::Stitcher() {
    alignment_params.match = 20;
    alignment_params.mismatch = 80;
    alignment_params.gap_open[0] = 60;
    alignment_params.gap_extend[0] = 30;
    alignment_params.gap_open[1] = 800;
    alignment_params.gap_extend[1] = 5;
    alignment_params.gap_open[2] = 2500;
    alignment_params.gap_extend[2] = 1;
}

void Stitcher::subalign(const SubGraphInfo& extraction1,
                        const SubGraphInfo& extraction2,
                        Alignment& stitched, bool only_deletion_alns) const {
    
    // TODO: it's wasteful doing this every time, but i expect it not to dominate
    
    // find out what gap size makes you need the larger gap parameters
    std::vector<size_t> cutoffs(alignment_params.gap_extend.size() - 1);
    for (size_t i = 1; i < alignment_params.gap_open.size(); ++i) {
        // the algebra is much easier if we get to assume a fixed order
        if (alignment_params.gap_open[i - 1] > alignment_params.gap_open[i] ||
            alignment_params.gap_extend[i - 1] < alignment_params.gap_extend[i]) {
            throw std::runtime_error("Affine gap parameters must be increasing in gap open penalty and decreasing in gap extend penalty");
        }
        
        auto diff_open = alignment_params.gap_open[i] - alignment_params.gap_open[i - 1];
        auto diff_extend = alignment_params.gap_extend[i - 1] - alignment_params.gap_extend[i];
        
        // round up to the nearest integer
        cutoffs[i - 1] = (diff_open + diff_extend - 1) / diff_extend;
    }
    
    Alignment inter_aln;
    size_t c = 0;
    while (c < cutoffs.size() &&
           extraction1.subgraph.node_size() > cutoffs[c] &&
           extraction2.subgraph.node_size() > cutoffs[c]) {
        ++c;
    }
    
    if (c == 0) {
        auto trunc_params = truncate_parameters<3, 1>(alignment_params);
        inter_aln = std::move(do_alignment(extraction1, extraction2, only_deletion_alns, trunc_params));
    }
    else if (c == 1) {
        auto trunc_params = truncate_parameters<3, 2>(alignment_params);
        inter_aln = std::move(do_alignment(extraction1, extraction2, only_deletion_alns, trunc_params));
    }
    else {
        inter_aln = std::move(do_alignment(extraction1, extraction2, only_deletion_alns, alignment_params));
    }
    
    translate(inter_aln, extraction1.back_translation, extraction2.back_translation);
    
    if (debug) {
        std::cerr << "inter-anchor alignment:\n";
        for (auto& ap : inter_aln) {
            std::cerr << '\t' << (int64_t) ap.node_id1 << '\t' << (int64_t) ap.node_id2 << '\n';
        }
    }
    
    for (auto aln_pair : inter_aln) {
        stitched.push_back(aln_pair);
    }
}

void Stitcher::do_instrument(const SubGraphInfo& extraction1,
                             const SubGraphInfo& extraction2,
                             int64_t score1, int64_t score2,
                             double dur1, double dur2) const {
    
    int64_t min1 = -1;
    int64_t max1 = -1;
    int64_t min2 = -1;
    int64_t max2 = -1;
    if (!extraction1.back_translation.empty()) {
        min1 = *std::min_element(extraction1.back_translation.begin(),
                                 extraction1.back_translation.end());
        max1 = *std::max_element(extraction1.back_translation.begin(),
                                 extraction1.back_translation.end());
    }
    if (!extraction2.back_translation.empty()) {
        min2 = *std::min_element(extraction2.back_translation.begin(),
                                 extraction2.back_translation.end());
        max2 = *std::max_element(extraction2.back_translation.begin(),
                                 extraction2.back_translation.end());
    }
    int64_t min_penalty = std::numeric_limits<int64_t>::max();
    int64_t approx_gap = std::abs((int64_t) (extraction1.subgraph.node_size() - extraction2.subgraph.node_size()));
    for (size_t i = 0; i < alignment_params.gap_open.size(); ++i) {
        min_penalty = min<int64_t>(min_penalty,
                                   alignment_params.gap_extend[i] * approx_gap + alignment_params.gap_open[i]);
    }
    int64_t approx_max_match = alignment_params.match * min(extraction1.subgraph.node_size(),
                                                            extraction2.subgraph.node_size());
    int64_t approx_max_score = approx_max_match - min_penalty;
    std::cerr << '#' << '\t' << extraction1.subgraph.node_size() << '\t' << extraction2.subgraph.node_size() << '\t' << ((extraction1.subgraph.node_size() + 1) * (extraction2.subgraph.node_size() + 1)) << '\t' << score1 << '\t' << score2 << '\t' << approx_max_match << '\t' << approx_max_score << '\t' << min1 << '\t' << max1 << '\t' << min2 << '\t' << max2 << '\t' << dur1 << '\t' << dur2 << '\n';
}



std::vector<std::pair<size_t, size_t>> Stitcher::identify_despecification_partition(const std::vector<anchor_t>& anchors) const {
    
    static const bool debug = false;
    
    static const double inf = std::numeric_limits<double>::max();
    static const double mininf = std::numeric_limits<double>::lowest();
    
    // find windows that only include one SV indel
    // TODO: ignores indels that fall on either end of the anchor chain
    std::vector<std::pair<int64_t, int64_t>> search_limit(anchors.size(), std::pair<int64_t, int64_t>(0, 0));
    int64_t prev_indel = -1;
    int64_t before_prev_indel = -1;
    for (size_t i = 0; i < anchors.size(); ++i) {
        if (i != 0 && std::abs(anchors[i].gap_before) >= min_indel_fuzz_length) {
            before_prev_indel = prev_indel;
            prev_indel = i;
        }
        if (before_prev_indel != -1 && prev_indel != -1) {
            search_limit[i].first = before_prev_indel + 1;
            search_limit[i].second = std::min<int64_t>(i, prev_indel + 1);
        }
        else if (prev_indel != -1) {
            search_limit[i].first = std::min<int64_t>(1, i);
            search_limit[i].second = std::min<int64_t>(prev_indel + 1, i);
        }
        if (debug) {
            std::cerr << i << ", prev " << prev_indel << ", before " << before_prev_indel << '\n';
        }
    }
    
    if (debug) {
        std::cerr << "search limit\n";
        for (size_t i = 0; i < search_limit.size(); ++i) {
            std::cerr << '\t' << i << '\t' << search_limit[i].first << '\t' << search_limit[i].second << '\n';
        }
    }
    
    std::vector<double> score_prefix_sum(anchors.size() + 1, 0.0);
    for (size_t i = 0; i < anchors.size(); ++i) {
        score_prefix_sum[i + 1] = score_prefix_sum[i] + anchors[i].score;
    }
    
    if (debug) {
        std::cerr << "score prefix sum\n";
        for (size_t i = 0; i < score_prefix_sum.size(); ++i) {
            std::cerr << '\t' << i << '\t' << score_prefix_sum[i] << '\n';
        }
    }
    
    std::vector<double> score_index_key(anchors.size(), mininf);
    for (size_t i = 1; i < anchors.size(); ++i) {
        score_index_key[i] = score_prefix_sum[i] + indel_fuzz_score_proportion * anchors[i - 1].score;
    }
    
    if (debug) {
        std::cerr << "score index key\n";
        for (size_t i = 0; i < score_index_key.size(); ++i) {
            std::cerr << '\t' << i << '\t' << score_index_key[i] << '\n';
        }
    }
    
    // records of (index, score limit, (number of despecified indels, first index of final despecification))
    // index only included in key to deduplicate key pairs
    std::vector<std::tuple<int64_t, double, std::tuple<int64_t, double, size_t>>> search_tree_data;
    search_tree_data.reserve(anchors.size() + 1);
    for (size_t i = 0; i < anchors.size(); ++i) {
        search_tree_data.emplace_back(i, score_index_key[i], std::tuple<int64_t, double, size_t>(0, 0.0, 0));
    }
    OrthogonalMaxSearchTree<int64_t, double, std::tuple<int64_t, double, size_t>> search_tree(search_tree_data);
    if (debug) {
        std::cerr << "search tree:\n";
        for (auto r : search_tree) {
            std::cerr << '\t' << std::get<0>(r) << '\t' <<  std::get<1>(r) << ":\t" <<  std::get<0>(std::get<2>(r)) << ", " << std::get<1>(std::get<2>(r)) << ", " << std::get<2>(std::get<2>(r)) << '\n';
        }
    }
    
    // records of (num indels despecified, score removed) for excluded, included
    std::vector<std::pair<std::tuple<int64_t, double, size_t>, std::tuple<int64_t, double, size_t>>> dp(anchors.size() + 1,
                                                                                                        std::make_pair(std::tuple<int64_t, double, size_t>(-1, 0.0, 0),
                                                                                                                       std::tuple<int64_t, double, size_t>(-1, 0.0, 0)));
    // back pointers for the score with the current item included
    std::vector<size_t> backpointer(dp.size(), -1);
    
    // boundary conditions
    std::get<0>(dp.front().first) = 0;
    
    size_t opt_idx = 0;
    
    // note: we step early because we require the removed intervals to be brackted by kept anchors
    for (size_t i = 1; i + 1 < dp.size(); ++i) {
        
        dp[i].first = std::max(dp[i - 1].first, dp[i - 1].second);
        
        double score_query_key = score_prefix_sum[i] - indel_fuzz_score_proportion * anchors[i].score;
        
        if (debug) {
            std::cerr << "iter " << i << ", query with index limit " << search_limit[i].first << ", " << search_limit[i].second << ", score key interval " << score_query_key << ", " << inf << '\n';
        }
        
        auto max_it = search_tree.range_max(search_limit[i].first, search_limit[i].second,
                                            score_query_key, inf);
        
        if (max_it != search_tree.end()) {
            if (debug) {
                std::cerr << "got max " << std::get<0>(std::get<2>(*max_it)) << ", " << std::get<1>(std::get<2>(*max_it)) << ", " << std::get<2>(std::get<2>(*max_it)) << " with keys " << std::get<0>(*max_it) << " and " << std::get<1>(*max_it)<< '\n';
            }
            dp[i].second = std::make_tuple(std::get<0>(std::get<2>(*max_it)) + 1,
                                           std::get<0>(std::get<2>(*max_it)) - score_prefix_sum[i] + score_prefix_sum[std::get<0>(*max_it)],
                                           i);
            
            backpointer[i] = std::get<0>(*max_it);
            if (debug) {
                std::cerr << "DP value " << std::get<0>(dp[i].second) << ", " << std::get<1>(dp[i].second) << ", " << std::get<2>(dp[i].second) << " with backpointer to " << backpointer[i] << '\n';
            }
            
            if (dp[i].second > dp[opt_idx].second) {
                opt_idx = i;
            }
        }
        
        auto it = search_tree.find(i, score_index_key[i]);
        search_tree.update(it, std::make_tuple(std::get<0>(dp[i].first),
                                               std::get<1>(dp[i].first),
                                               i));
        
        if (debug) {
            std::cerr << "update key " << std::get<0>(*it) << ", "  << std::get<1>(*it) << " to "<< std::get<0>(std::get<2>(*it)) << ", " << std::get<1>(std::get<2>(*it)) << ", " << std::get<2>(std::get<2>(*it)) << '\n';
            std::cerr << "search tree state:\n";
            for (auto r : search_tree) {
                std::cerr << '\t' << std::get<0>(r) << '\t' <<  std::get<1>(r) << ":\t" <<  std::get<0>(std::get<2>(r)) << ", " << std::get<1>(std::get<2>(r)) << ", " << std::get<2>(std::get<2>(r)) << '\n';
            }
            std::cerr << "DP state:\n";
            for (size_t i = 0; i < dp.size(); ++i) {
                std::cerr << i << "\texcl " << (int64_t) std::get<0>(dp[i].first) << ", " << std::get<1>(dp[i].first) << ", " << std::get<2>(dp[i].first) << "\tincl " << (int64_t) std::get<0>(dp[i].second) << ", " <<  std::get<1>(dp[i].second) << ", " <<  std::get<2>(dp[i].second) << "\tback " << (int64_t) backpointer[i] << '\n';
            }
            std::cerr << "opt occurs at " << opt_idx << '\n';
        }
    }
    
    if (debug) {
        std::cerr << "final DP state:\n";
        for (size_t i = 0; i < dp.size(); ++i) {
            std::cerr << i << "\texcl " << (int64_t) std::get<0>(dp[i].first) << ", " << std::get<1>(dp[i].first) << ", " << std::get<2>(dp[i].first) << "\tincl " << (int64_t) std::get<0>(dp[i].second) << ", " <<  std::get<1>(dp[i].second) << ", " <<  std::get<2>(dp[i].second) << "\tback " << (int64_t) backpointer[i] << '\n';
        }
        std::cerr << "opt occurs at " << opt_idx << '\n';
    }
    
    return traceback(dp, backpointer, opt_idx);
}

void Stitcher::despecify_indel_breakpoints(std::vector<anchor_t>& anchors) const {
    
    auto removal_partition = identify_despecification_partition(anchors);
    
    static const bool instrument = false;
    if (instrument) {
        std::cerr << "removing anchors to despecify SV indel breakpoints:\n";
        for (auto interval : removal_partition) {
            for (size_t i = interval.first - 1; i <= interval.second; ++i) {
                const auto& anchor = anchors[i];
                char rem = (i == interval.first - 1 || i == interval.second) ? 'k' : 'r';
                std::cerr << '\\' << '\t' << i << '\t' << rem << '\t' << anchor.walk1.size() << '\t' << anchor.count1 << '\t' << anchor.count2 << '\t' << (anchor.count1 * anchor.count2) << '\t' << anchor.score << '\t' << anchor.walk1.front() << '\t' << anchor.walk1.back() << '\t' << anchor.walk2.front() << '\t' << anchor.walk2.back() << '\t' << anchor.gap_before << '\t'<< anchor.gap_score_before << '\t' << anchor.gap_after << '\t' << anchor.gap_score_after << '\n';
                
            }
        }
    }
    
    // filter out the anchors identified by the partition
    size_t removed = 0;
    size_t d = 0;
    int64_t gap = 0;
    double gap_score = 0.0;
    for (size_t i = 0; i < anchors.size(); ++i) {
        if (d < removal_partition.size() && i >= removal_partition[d].first && i < removal_partition[d].second) {
            gap += anchors[i].gap_before;
            gap_score += anchors[i].gap_score_before;
            ++removed;
        }
        else if (removed != 0) {
            anchors[i - removed] = std::move(anchors[i]);
        }
        if (d < removal_partition.size() && i == removal_partition[d].second) {
            anchors[i - removed - 1].gap_after = gap;
            anchors[i - removed - 1].gap_score_after = gap_score;
            anchors[i - removed].gap_before = gap;
            anchors[i - removed].gap_score_before = gap_score;
            gap = 0;
            gap_score = 0.0;
            ++d;
        }
    }
    
    if (removed != 0) {
        anchors.resize(anchors.size() - removed);
    }
}

}
