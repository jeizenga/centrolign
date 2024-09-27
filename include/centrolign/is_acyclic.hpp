#ifndef centrolign_is_acyclic_hpp
#define centrolign_is_acyclic_hpp

#include <vector>

namespace centrolign {

// Kahn's algorithm
// TODO: unify this with topological order?
template<class Graph>
bool is_acyclic(const Graph& graph) {
    
    // find sources
    std::vector<uint64_t> stack;
    std::vector<size_t> in_degree(graph.node_size());
    for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
        in_degree[node_id] = graph.previous_size(node_id);
        if (in_degree[node_id] == 0) {
            stack.push_back(node_id);
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
            }
        }
    }
    
    // did we (conceptually) eliminate all of the edges?
    for (auto deg : in_degree) {
        if (deg) {
            return false;
        }
    }
    
    return true;
}

}

#endif /* centrolign_is_acyclic_hpp */
