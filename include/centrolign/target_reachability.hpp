#ifndef centrolign_target_reachability_hpp
#define centrolign_target_reachability_hpp

#include "centrolign/topological_order.hpp"
#include "centrolign/utility.hpp"

#include <vector>
#include <utility>
#include <limits>

namespace centrolign {


// return whether each node can reach at least one of the targets
template<class Graph>
std::vector<bool> target_reachability(const Graph& graph, const std::vector<uint64_t>& targets) {
    
    std::vector<bool> reachable(graph.node_size(), false);
    
    for (auto node_id : targets) {
        reachable[node_id] = true;
    }
    
    auto order = topological_order(graph);
    for (uint64_t node_id : ReverseForEachAdapter<std::vector<uint64_t>>(order)) {
        for (auto next_id : graph.next(node_id)) {
            reachable[node_id] = reachable[node_id] || reachable[next_id];
        }
    }
    
    return reachable;
}

}

#endif /* centrolign_target_reachability_hpp */
