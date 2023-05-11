#include "centrolign/anchorer.hpp"

#include <algorithm>

#include "centrolign/topological_order.hpp"

namespace centrolign {

using namespace std;

const bool Anchorer::debug_anchorer = false;

std::vector<anchor_t> Anchorer::exhaustive_chain_dp(std::vector<Anchorer::anchor_set_t>& anchor_sets,
                                                    const ChainMerge& chain_merge1,
                                                    const ChainMerge& chain_merge2) const {
    
    // make a graph of match nodes
    AnchorGraph anchor_graph;
    for (size_t i = 0; i < anchor_sets.size(); ++i) {
        
        const auto& anchor_set = anchor_sets[i];
        
        double weight = anchor_weight(anchor_set.walks1.size(), anchor_set.walks2.size(),
                                      anchor_set.walks1.front().size());
        
        for (size_t idx1 = 0; idx1 < anchor_set.walks1.size(); ++idx1) {
            for (size_t idx2 = 0; idx2 < anchor_set.walks2.size(); ++idx2) {
                anchor_graph.add_node(i, idx1, idx2, weight);
            }
        }
    }
    
    // TODO: with no edge costs and positive weights, we can always do DP over the transitive
    // reduction, so we could speed this up by figuring out a better way to reduce transitive edges
    
    // add all possible edges
    for (uint64_t node_id1 = 0; node_id1 < anchor_graph.node_size(); ++node_id1) {
        for (uint64_t node_id2 = 0; node_id2 < anchor_graph.node_size(); ++node_id2) {
            size_t set1, idx11, idx21, set2, idx12, idx22;
            std::tie(set1, idx11, idx21) = anchor_graph.label(node_id1);
            std::tie(set2, idx12, idx22) = anchor_graph.label(node_id2);
            
            auto& anchor_set1 = anchor_sets[set1];
            auto& anchor_set2 = anchor_sets[set2];
            
            if (chain_merge1.reachable(anchor_set1.walks1[idx11].back(),
                                       anchor_set2.walks1[idx12].front()) &&
                chain_merge2.reachable(anchor_set1.walks2[idx21].back(),
                                       anchor_set2.walks2[idx22].front())) {
                
                anchor_graph.add_edge(node_id1, node_id2);
            }
        }
    }
    
    // get heaviest path and convert into a chain
    std::vector<anchor_t> chain;
    for (auto node_id : anchor_graph.heaviest_weight_path()) {
        size_t set, idx1, idx2;
        std::tie(set, idx1, idx2) = anchor_graph.label(node_id);
        chain.emplace_back();
        auto& anchor_set = anchor_sets[set];
        auto& chain_node = chain.back();
        chain_node.walk1 = std::move(anchor_set.walks1[idx1]);
        chain_node.walk2 = std::move(anchor_set.walks2[idx2]);
        chain_node.count1 = anchor_set.walks1.size();
        chain_node.count2 = anchor_set.walks2.size();
    }
    return chain;
}

uint64_t Anchorer::AnchorGraph::add_node(size_t set, size_t idx1, size_t idx2, double weight) {
    nodes.emplace_back(set, idx1, idx2, weight);
    return nodes.size() - 1;
}

void Anchorer::AnchorGraph::add_edge(uint64_t from_id, uint64_t to_id) {
    nodes[from_id].edges.push_back(to_id);
    ++nodes[to_id].in_degree;
}

const vector<size_t> Anchorer::AnchorGraph::next(uint64_t node_id) const {
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

vector<uint64_t> Anchorer::AnchorGraph::heaviest_weight_path() const {
    
    // set up DP structures
    std::vector<size_t> backpointer(nodes.size(), -1);
    std::vector<double> dp(nodes.size(), 0.0);
    
    // do dynamic programming
    uint64_t max_id = -1;
    double max_weight = 0.0;
    for (auto node_id : topological_order(*this)) {
        
        // add nodes own weight in
        const auto& node = nodes[node_id];
        auto& dp_here = dp[node_id];
        dp_here += node.weight;
        
        if (debug_anchorer) {
            cerr << "DP at " << node_id << " set to " << dp_here << '\n';
        }
        
        if (dp_here > max_weight) {
            max_id = node_id;
            max_weight = dp_here;
        }
        
        // propagate forward
        for (auto next_id : node.edges) {
            if (dp_here > dp[next_id]) {
                dp[next_id] = dp_here;
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


}
