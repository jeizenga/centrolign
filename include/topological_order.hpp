#ifndef centrolign_topological_order_hpp
#define centrolign_topological_order_hpp

#include <vector>

namespace centrolign {

// Kahn's algorithm
template<class Graph>
std::vector<uint64_t> topological_order(const Graph& graph) {
    
    std::vector<uint64_t> order;
    order.reserve(graph.node_size());
    
    std::vector<uint64_t> stack;
    std::vector<size_t> in_degree(graph.node_size());
    for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
        in_degree[node_id] = graph.previous_size(i);
        if (in_degree[node_id] == 0) {
            stack.push_back(node_id);
            order.push_back(node_id);
        }
    }
    
    while (!stack.empty()) {
        uint64_t node_id = stack.back();
        stack.pop_back();
        for (uint64_t next_id : graph.next(node_id)) {
            --in_degree[next_id];
            if (in_degree[next_id] == 0) {
                stack.push_back(next_id);
                order.push_back(next_id);
            }
        }
    }
}

}

#endif /* centrolign_topological_order_hpp */
