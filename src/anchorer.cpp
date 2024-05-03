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
        if (i != 0) {
            // copy the original anchor
            merged.emplace_back(std::move(anchors[i - 1]));
        }
        
        const auto& back_trans1 = fill_in_graphs[i].first.back_translation;
        const auto& back_trans2 = fill_in_graphs[i].second.back_translation;
        
        for (const auto& anchor : fill_in_chains[i]) {
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
        }
    }
    
    anchors = std::move(merged);
}

void Anchorer::update_mask(const std::vector<match_set_t>& matches, const std::vector<anchor_t>& chain,
                           std::unordered_set<std::tuple<size_t, size_t, size_t>>& masked_matches) const {
    
    logging::log(logging::Debug, "Updating mask");
    
    // record the node pairings that we're going to mask
    std::unordered_map<uint64_t, uint64_t> paired_node_ids;
    for (const auto& anchor : chain) {
        for (size_t i = 0; i < anchor.walk1.size(); ++i) {
            paired_node_ids[anchor.walk1[i]] = anchor.walk2[i];
        }
    }
    
    for (size_t i = 0; i < matches.size(); ++i) {
        
        const auto& match_set = matches[i];
        
        // for each position in walk, for each node ID, the walk indexes that it matches
        std::vector<std::unordered_map<uint64_t, std::vector<size_t>>> walk2_node(match_set.walks1.front().size());
        for (size_t k = 0; k < match_set.walks2.size(); ++k) {
            const auto& walk2 = match_set.walks2[k];
            for (size_t l = 0; l < walk2.size(); ++l) {
                walk2_node[l][walk2[l]].push_back(k);
            }
        }
        
        // check for walk2s that have a paired node ID in the matching position
        for (size_t j = 0; j < match_set.walks1.size(); ++j) {
            const auto& walk1 = match_set.walks1[j];
            for (size_t l = 0; l < walk1.size(); ++l) {
                auto it = paired_node_ids.find(walk1[l]);
                if (it != paired_node_ids.end()) {
                    auto it2 = walk2_node[l].find(it->second);
                    if (it2 != walk2_node[l].end()) {
                        for (auto k : it2->second) {
                            masked_matches.emplace(i, j, k);
                        }
                    }
                }
            }
        }
    }
}


}
