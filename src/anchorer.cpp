#include "centrolign/anchorer.hpp"

#include <algorithm>

#include "centrolign/topological_order.hpp"

namespace centrolign {

using namespace std;

const bool Anchorer::debug_anchorer = false;
const bool Anchorer::basic_logging = true;

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
    
    if (debug_anchorer || basic_logging) {
        cerr << "heaviest weight path consists of " << traceback.size() << " anchors, which have total weight " << max_weight << '\n';
    }
    
    return traceback;
}


}
