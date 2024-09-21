#ifndef centrolign_connected_components_hpp
#define centrolign_connected_components_hpp

#include <vector>
#include <cstdint>

namespace centrolign {

// returns the (weakly) connected components of the graph
template<class Graph>
std::vector<std::vector<uint64_t>> connected_components(const Graph& graph) {
    
    std::vector<std::vector<uint64_t>> components;
    
    std::vector<bool> traversed(graph.node_size(), false);
    
    for (uint64_t start_id = 0; start_id < graph.node_size(); ++start_id) {
        if (traversed[start_id]) {
            continue;
        }
        // we haven't seen this node before, start a new component
        components.emplace_back();
        auto& component = components.back();
        
        // init the stack
        std::vector<uint64_t> stack(1, start_id);
        traversed[start_id] = true;
        
        // do DFS
        while (!stack.empty()) {
            
            auto node_id = stack.back();
            stack.pop_back();
            component.push_back(node_id);
            
            for (bool left : {true, false}) {
                for (auto next_id : (left ? graph.previous(node_id) : graph.next(node_id))) {
                    if (!traversed[next_id]) {
                        traversed[next_id] = true;
                        stack.push_back(next_id);
                    }
                }
            }
        }
    }
    
    return components;
}


}

#endif /* centrolign_connected_components_hpp */
