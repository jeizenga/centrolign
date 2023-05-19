#ifndef centrolign_count_walks_hpp
#define centrolign_count_walks_hpp

#include <vector>
#include <cstdint>

#include "centrolign/topological_order.hpp"
#include "centrolign/utility.hpp"

namespace centrolign {

// returns the number of source-to-sink walks in the graph, or
// numeric_limit<uint64_t>::max() if the number would overflow a 64 bit int
template<class Graph>
uint64_t count_walks(const Graph& graph) {
    
    uint64_t total = 0;
    
    std::vector<uint64_t> dp(graph.node_size(), 0);
    
    for (auto node_id : topological_order(graph)) {
        
        if (graph.previous_size(node_id) == 0) {
            // base case
            dp[node_id] = 1;
        }
        
        if (graph.next_size(node_id) == 0) {
            total = sat_add(total, dp[node_id]);
        }
        else {
            for (auto next_id : graph.next(node_id)) {
                dp[next_id] = sat_add(dp[next_id], dp[node_id]);
            }
        }
    }
    return total;
}

}

#endif /* centrolign_count_walks_hpp */
