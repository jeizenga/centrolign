#ifndef centrolign_minmax_distance_hpp
#define centrolign_minmax_distance_hpp

#include "centrolign/topological_order.hpp"

#include <vector>
#include <utility>
#include <limits>

namespace centrolign {


// return the minimum and maximum distance to each node from either a neighbor
// with in-degree 0 or a set of optional, provided source nodes
template<class Graph>
std::vector<std::pair<int64_t, int64_t>> minmax_distance(const Graph& graph,
                                                         const std::vector<uint64_t>* source_nodes = nullptr) {
    
    
    const bool debug = false;
    
    std::vector<std::pair<int64_t, int64_t>> dp(graph.node_size(),
                                                std::pair<int64_t, int64_t>(std::numeric_limits<int64_t>::max(), -1));
    
    if (source_nodes) {
        for (auto node_id : *source_nodes) {
            dp[node_id] = std::pair<int64_t, int64_t>(0, 0);
        }
    }
    else {
        for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
            if (graph.previous_size(node_id) == 0) {
                dp[node_id] = std::pair<int64_t, int64_t>(0, 0);
            }
        }
    }
    
    if (debug) {
        if (source_nodes) {
            std::cerr << "source nodes:";
            for (auto node_id : *source_nodes) {
                std::cerr << ' ' << node_id;
            }
            std::cerr << '\n';
        }
        std::cerr << "initial conditions:\n";
        for (auto& d : dp) {
            std::cerr << d.first << '\t' << d.second << '\n';
        }
    }
    
    for (uint64_t node_id : topological_order(graph)) {
        auto& dp_here = dp[node_id];
        if (dp_here.first != std::numeric_limits<int64_t>::max()) {
            for (auto next_id : graph.next(node_id)) {
                auto& dp_next = dp[next_id];
                dp_next.first = std::min(dp_next.first, dp_here.first + 1);
                dp_next.second = std::max(dp_next.second, dp_here.second + 1);
            }
        }
    }
    
    
    if (debug) {
        std::cerr << "after DP:\n";
        for (auto& d : dp) {
            std::cerr << d.first << '\t' << d.second << '\n';
        }
    }
    
    return dp;
}

}

#endif /* centrolign_minmax_distance_hpp */
