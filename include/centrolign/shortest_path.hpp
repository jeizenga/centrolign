#ifndef centrolign_shortest_path_hpp
#define centrolign_shortest_path_hpp

#include "centrolign/topological_order.hpp"

#include <vector>
#include <utility>
#include <limits>

namespace centrolign {


// the shortest path between two nodes, including the two nodes themselves,
// linear in the size of the graph
template<class Graph>
std::vector<uint64_t> shortest_path(const Graph& graph,
                                    uint64_t node_id1, uint64_t node_id2);

// the shortest path between any pair of nodes from the two sets
template<class Graph>
std::vector<uint64_t> shortest_path(const Graph& graph,
                                    const std::vector<uint64_t>& node_ids1,
                                    const std::vector<uint64_t>& node_ids2);




/*
 *  Template implementations
 */

template<class Graph>
std::vector<uint64_t> shortest_path(const Graph& graph,
                                    const std::vector<uint64_t>& node_ids1,
                                    const std::vector<uint64_t>& node_ids2)
{
    assert(!node_ids1.empty() && !node_ids2.empty());
    
    const bool debug = false;
    if (debug) {
        std::cerr << "computing shortest path between";
        for (auto n : node_ids1) {
            std::cerr << ' ' << n;
        }
        std::cerr << " and";
        for (auto n : node_ids2) {
            std::cerr << ' ' << n;
        }
        std::cerr << '\n';
    }
    
    std::vector<size_t> dp(graph.node_size(), std::numeric_limits<int64_t>::max());
    
    // dynamic programming
    for (auto node_id1 : node_ids1) {
        dp[node_id1] = 0;
    }
    for (uint64_t node_id : topological_order(graph)) {
        size_t dist_thru = dp[node_id] + graph.label_size(node_id);
        for (auto next_id : graph.next(node_id)) {
            dp[next_id] = std::min(dp[next_id], dist_thru);
        }
    }
    
    if (debug) {
        std::cerr << "dp:\n";
        for (size_t i = 0; i < dp.size(); ++i) {
            std::cerr << i << ": " << dp[i] << '\n';
        }
    }
    
    uint64_t node_id2 = -1;
    for (auto node_id : node_ids2) {
        if (dp[node_id] != std::numeric_limits<int64_t>::max() &&
            (node_id2 == -1 || dp[node_id] < dp[node_id2])) {
            node_id2 = node_id;
        }
    }
    
    // traceback
    std::vector<uint64_t> path;
    if (node_id2 != -1) {
        path.push_back(node_id2);
        while (dp[path.back()] != 0) {
            for (auto prev_id : graph.previous(path.back())) {
                if (dp[prev_id] + graph.label_size(prev_id) == dp[path.back()]) {
                    path.push_back(prev_id);
                    break;
                }
            }
        }
        
        std::reverse(path.begin(), path.end());
    }
    
    return path;
}

template<class Graph>
std::vector<uint64_t> shortest_path(const Graph& graph,
                                    uint64_t node_id1, uint64_t node_id2) {
    return shortest_path(graph, std::vector<uint64_t>(1, node_id1),
                         std::vector<uint64_t>(1, node_id2));
}

}

#endif /* centrolign_shortest_path_hpp */
