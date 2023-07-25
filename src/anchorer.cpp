#include "centrolign/anchorer.hpp"

#include <algorithm>

#include "centrolign/topological_order.hpp"

namespace centrolign {

using namespace std;

const bool Anchorer::debug_anchorer = false;


uint64_t Anchorer::AnchorGraph::add_node(size_t set, size_t idx1, size_t idx2, double weight) {
    nodes.emplace_back(set, idx1, idx2, weight);
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

std::vector<anchor_t> Anchorer::traceback_sparse_dp(const std::vector<match_set_t>& match_sets,
                                                    const std::vector<std::vector<std::vector<dp_entry_t>>>& dp,
                                                    bool suppress_verbose_logging) const {
    
    if (debug_anchorer) {
        std::cerr << "finding optimum\n";
    }
    
    // find the optimum dynamic programming values
    dp_entry_t opt(std::numeric_limits<double>::lowest(), -1, -1, -1);
    for (size_t set = 0; set < dp.size(); ++set) {
        const auto& set_dp = dp[set];
        for (size_t i = 0; i < set_dp.size(); ++i) {
            const auto& dp_row = set_dp[i];
            for (size_t j = 0; j < dp_row.size(); ++j) {
                if (std::get<0>(dp_row[j]) > std::get<0>(opt)) {
                    opt = dp_entry_t(std::get<0>(dp_row[j]), set, i, j);
                }
            }
        }
    }
    
    
    if (debug_anchorer) {
        std::cerr << "doing traceback\n";
    }
    
    // traceback into a chain
    std::vector<anchor_t> anchors;
    auto here = opt;
    while (std::get<1>(here) != -1) {
        
        if (debug_anchorer) {
            std::cerr << "following traceback to set " << std::get<1>(here) << ", walk pair " << std::get<2>(here) << " " << std::get<3>(here) << '\n';
        }
        
        // grab the anchors that we used from their set
        auto& match_set = match_sets[std::get<1>(here)];
        anchors.emplace_back();
        auto& anchor = anchors.back();
        anchor.walk1 = match_set.walks1[std::get<2>(here)];
        anchor.count1 = match_set.count1;
        anchor.walk2 = match_set.walks2[std::get<3>(here)];
        anchor.count2 = match_set.count2;
        
        // follow the backpointer from the DP structure
        here = dp[std::get<1>(here)][std::get<2>(here)][std::get<3>(here)];
    }
    
    // take out of reverse order
    std::reverse(anchors.begin(), anchors.end());
    
    if (debug_anchorer) {
        std::cerr << "completed sparse chaining\n";
    }
    
    if (!suppress_verbose_logging) {
        logging::log(logging::Debug, "Optimal chain consists of " + std::to_string(anchors.size()) + " matches with score " + std::to_string(std::get<0>(opt)));
    }
    
    return anchors;
}

}
