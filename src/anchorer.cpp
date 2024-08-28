#include "centrolign/anchorer.hpp"

#include <algorithm>

#include "centrolign/topological_order.hpp"
#include "centrolign/minmax_distance.hpp"

namespace centrolign {

using namespace std;

const bool Anchorer::debug_anchorer = false;

pair<int64_t, int64_t> Extractor::source_sink_minmax(const SubGraphInfo& extraction) {
    
    auto dists = minmax_distance(extraction.subgraph, &extraction.sources);
    pair<int64_t, int64_t> mm(numeric_limits<int64_t>::max(), -1);
    for (auto node_id : extraction.sinks) {
        mm.first = min(mm.first, dists[node_id].first);
        mm.second = max(mm.second, dists[node_id].second);
    }
    return mm;
}

std::vector<size_t> Extractor::get_logging_indexes(size_t size) {
    std::vector<size_t> logging_indexes;
    for (size_t i = 1; i < 10; ++i) {
        logging_indexes.push_back((size * i) / 10);
    }
    auto end = std::unique(logging_indexes.begin(), logging_indexes.end());
    logging_indexes.resize(end - logging_indexes.begin());
    return logging_indexes;
}

uint64_t Anchorer::AnchorGraph::add_node(size_t set, size_t idx1, size_t idx2, double weight,
                                         double initial_weight, double final_weight) {
    nodes.emplace_back(set, idx1, idx2, weight, initial_weight, final_weight);
    return nodes.size() - 1;
}

void Anchorer::AnchorGraph::add_edge(uint64_t from_id, uint64_t to_id, double weight) {
    nodes[from_id].edges.push_back(to_id);
    nodes[from_id].edge_weights.push_back(weight);
    ++nodes[to_id].in_degree;
}

const vector<size_t>& Anchorer::AnchorGraph::next(uint64_t node_id) const {
    return nodes[node_id].edges;
}

size_t Anchorer::AnchorGraph::next_size(uint64_t node_id) const {
    return nodes[node_id].edges.size();
}

size_t Anchorer::AnchorGraph::previous_size(uint64_t node_id) const {
    return nodes[node_id].in_degree;
}

size_t Anchorer::AnchorGraph::node_size() const {
    return nodes.size();
}

tuple<size_t, size_t, size_t> Anchorer::AnchorGraph::label(uint64_t node_id) const {
    const auto& node = nodes[node_id];
    return make_tuple(node.set, node.idx1, node.idx2);
}

vector<uint64_t> Anchorer::AnchorGraph::heaviest_weight_path(double min_score) const {
    
    // set up DP structures
    std::vector<size_t> backpointer(nodes.size(), -1);
    std::vector<double> dp(nodes.size());
    for (size_t i = 0; i < nodes.size(); ++i) {
        dp[i] = nodes[i].initial_weight;
    }
    
    // do dynamic programming
    uint64_t max_id = -1;
    double max_weight = min_score;
    for (auto node_id : topological_order(*this)) {
        
        // add nodes own weight in
        const auto& node = nodes[node_id];
        auto& dp_here = dp[node_id];
        if (dp_here == numeric_limits<double>::lowest()) {
            continue;
        }
        dp_here += node.weight;
        
        if (debug_anchorer) {
            cerr << "DP at " << node_id << " set to " << dp_here << '\n';
        }
        
        if (node.final_weight != numeric_limits<double>::lowest() &&
            dp_here + node.final_weight > max_weight) {
            max_id = node_id;
            max_weight = dp_here + node.final_weight;
        }
        
        // propagate forward
        for (size_t i = 0; i < node.edges.size(); ++i) {
            auto next_id = node.edges[i];
            auto edge_weight = node.edge_weights[i];
            if (dp_here + edge_weight > dp[next_id]) {
                dp[next_id] = dp_here + edge_weight;
                backpointer[next_id] = node_id;
            }
        }
    }
    
    if (debug_anchorer) {
        cerr << "max DP occurs at " << max_id << " with DP value " << max_weight << '\n';
    }
    
    // traceback
    vector<uint64_t> traceback;
    if (max_id != -1) {
        traceback.push_back(max_id);
        while (backpointer[traceback.back()] != -1) {
            if (debug_anchorer) {
                cerr << "trace back to " << backpointer[traceback.back()] << '\n';
            }
            traceback.push_back(backpointer[traceback.back()]);
        }
        reverse(traceback.begin(), traceback.end());
    }
    
    if (debug_anchorer) {
        cerr << "heaviest weight path consists of " << traceback.size() << " anchors, which have total weight " << max_weight << '\n';
    }
    
    return traceback;
}


std::vector<size_t> Anchorer::assign_reanchor_budget(const std::vector<std::pair<SubGraphInfo, SubGraphInfo>>& fill_in_graphs) const {
    
    // the sum of the matrix sizes
    size_t total = 0;
    for (const auto& stitch_pair : fill_in_graphs) {
        total += (stitch_pair.first.subgraph.node_size() + 1) * (stitch_pair.second.subgraph.node_size() + 1);
    }
    
    std::vector<size_t> budget;
    budget.reserve(fill_in_graphs.size());
    for (const auto& stitch_pair : fill_in_graphs) {
        // proportionate to this matrix's contribution to the sum of matrix sizes
        size_t size = (stitch_pair.first.subgraph.node_size() + 1) * (stitch_pair.second.subgraph.node_size() + 1);
        budget.push_back(ceil(double(max_num_match_pairs) * double(size) / double(total)));
    }
    
    return budget;
}


void Anchorer::merge_fill_in_chains(std::vector<anchor_t>& anchors,
                                    const std::vector<std::vector<anchor_t>> fill_in_chains,
                                    const std::vector<std::pair<SubGraphInfo, SubGraphInfo>>& fill_in_graphs,
                                    const std::vector<std::vector<std::pair<size_t, std::pair<std::vector<size_t>, std::vector<size_t>>>>>& match_origin) const {
    
    std::vector<anchor_t> merged;
    
    assert(anchors.size() + 1 == fill_in_chains.size());
    
    for (size_t i = 0; i < fill_in_chains.size(); ++i) {
        int remember_idx;
        if (i != 0) {
            // copy the original anchor
            auto& anchor = anchors[i - 1];\
            if (!merged.empty()) {
                // update the gap before the original anchor
                anchor.gap_before = merged.back().gap_after;
                anchor.gap_score_before = merged.back().gap_score_after;
            }
            merged.emplace_back(std::move(anchor));
        }
        
        const auto& back_trans1 = fill_in_graphs[i].first.back_translation;
        const auto& back_trans2 = fill_in_graphs[i].second.back_translation;
        
        const auto& to_merge = fill_in_chains[i];
        for (size_t j = 0; j < to_merge.size(); ++j) {
            const auto& anchor = to_merge[j];
            if (j == 0 && !merged.empty()) {
                // update the gap after the original anchor
                merged.back().gap_score_after = anchor.gap_score_before;
                merged.back().gap_after = anchor.gap_before;
            }
            // copy the anchor in this in-between stitching chain
            merged.emplace_back();
            auto& translated = merged.back();
            translated.walk1.reserve(anchor.walk1.size());
            translated.walk2.reserve(anchor.walk2.size());
            // the node IDs must be translated from subgraph to main graph
            for (auto node_id : anchor.walk1) {
                translated.walk1.push_back(back_trans1[node_id]);
            }
            for (auto node_id : anchor.walk2) {
                translated.walk2.push_back(back_trans2[node_id]);
            }
            // we preserve the original queried count from the match index
            translated.count1 = anchor.count1;
            translated.count2 = anchor.count2;
            const auto& origin_set = match_origin[i][anchor.match_set];
            
            // back translate the identity of this anchor
            translated.match_set = origin_set.first;
            translated.idx1 = origin_set.second.first[anchor.idx1];
            translated.idx2 = origin_set.second.second[anchor.idx2];
            translated.gap_before = anchor.gap_before;
            translated.gap_score_before = anchor.gap_score_before;
            translated.gap_after = anchor.gap_after;
            translated.gap_score_after = anchor.gap_score_after;
        }
    }
    
    anchors = std::move(merged);
}


std::vector<std::pair<size_t, size_t>> Anchorer::identify_despecification_partition(const std::vector<anchor_t>& anchors) const {
    
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

void Anchorer::despecify_indel_breakpoints(std::vector<anchor_t>& anchors) const {
    
    auto partition = identify_despecification_partition(anchors);
    
    // filter out the anchors identified by the partition
    size_t removed = 0;
    size_t p = 0;
    for (size_t i = 0; i < anchors.size(); ++i) {
        if (p < partition.size() && i >= partition[p].first && i < partition[p].second) {
            ++removed;
        }
        else if (removed != 0) {
            anchors[i - removed] = std::move(anchors[i]);
        }
        if (p < partition.size() && i == partition[p].second) {
            ++p;
        }
    }
    
    if (removed != 0) {
        anchors.resize(anchors.size() - removed);
    }
}


}
