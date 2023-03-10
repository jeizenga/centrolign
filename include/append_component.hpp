#ifndef centrolign_append_component_hpp
#define centrolign_append_component_hpp

#include <cstdint>

#include "graph.hpp"

namespace centrolign {


// append a graph as a separate connected component onto another
template<class BGraph>
void append_component(BaseGraph& appending, const BGraph& component) {
    
    uint64_t prev_node_size = appending.node_size();
    
    // add the nodes
    for (uint64_t node_id = 0; node_id < component.node_size(); ++node_id) {
        appending.add_node(component.base(node_id));
    }
    
    // add the edges
    for (uint64_t node_id = 0; node_id < component.node_size(); ++node_id) {
        for (uint64_t next_id : component.next(node_id)) {
            // IDs are offset by a constant amount
            appending.add_edge(node_id + prev_node_size, next_id + prev_node_size);
        }
    }
}

}

#endif /* centrolign_append_component_hpp */
