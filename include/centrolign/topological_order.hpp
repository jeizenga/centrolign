#ifndef centrolign_topological_order_hpp
#define centrolign_topological_order_hpp

#include <vector>
#include <iostream>

namespace centrolign {

// Kahn's algorithm, returns vector of node IDs in topological order
template<class Graph>
std::vector<uint64_t> topological_order(const Graph& graph) {
    
    static bool debug_top_order = false;
    
    std::vector<uint64_t> order;
    order.reserve(graph.node_size());
    
    // find sources
    std::vector<uint64_t> stack;
    std::vector<size_t> in_degree(graph.node_size());
    for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
        in_degree[node_id] = graph.previous_size(node_id);
        if (in_degree[node_id] == 0) {
            stack.push_back(node_id);
            order.push_back(node_id);
        }
    }
    
    if (debug_top_order) {
        std::cerr << "initial in degrees:\n";
        for (size_t i = 0; i < in_degree.size(); ++i) {
            std::cerr << i << ": " << in_degree[i] << '\n';
        }
    }
    
    // (conceptually) remove outward edges and queue up any new sources
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
    
    if (debug_top_order) {
        std::cerr << "final order:\n";
        for (size_t i = 0; i < in_degree.size(); ++i) {
            std::cerr << i << ": " << order[i] << '\n';
        }
    }
    
    // ensure acyclicity
    assert(order.size() == graph.node_size());
    
    return order;
}

}

#endif /* centrolign_topological_order_hpp */
