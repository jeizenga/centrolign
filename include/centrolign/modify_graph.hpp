#ifndef centrolign_modify_graph_hpp
#define centrolign_modify_graph_hpp

#include <cstdint>

#include "centrolign/graph.hpp"

namespace centrolign {


// append a graph as a separate connected component onto another
template<class BGraph>
void append_component(BaseGraph& appending, const BGraph& component);

// add sentinel nodes to the end of paths
void add_path_end_sentinels(BaseGraph& graph, uint8_t sentinel_val = 5);


/*
 * Template implementations
 */

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
    
    // add the paths
    for (uint64_t path_id = 0; path_id < component.path_size(); ++path_id) {
        uint64_t new_path_id = component.add_path(component.path_name(path_id));
        for (uint64_t node_id : component.path(path_id)) {
            appending.extend_path(new_path_id, prev_node_size + node_id);
        }
    }
}

}

#endif /* centrolign_modify_graph_hpp */
