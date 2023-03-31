#ifndef centrolign_modify_graph_hpp
#define centrolign_modify_graph_hpp

#include <cstdint>

#include "centrolign/graph.hpp"

namespace centrolign {

// make simple graphs with encoded chars and an embedded path, but no sentinels
SequenceGraph make_sequence_graph(const std::string& name,
                                  const std::string& sequence);
BaseGraph make_base_graph(const std::string& name,
                          const std::string& sequence);

// append a graph as a separate connected component onto another, including
// embedded paths
template<class BGraph>
void append_component(BaseGraph& appending, const BGraph& component);

/*
 * A convenience struct to keep track of sentinel node information
 */
struct SentinelTableau {
    uint64_t src_id = -1;
    uint64_t snk_id = -1;
    char src_sentinel = 0;
    char snk_sentinel = 0;
};

// add sentinel nodes and construct a tableau
SentinelTableau add_sentinels(BaseGraph& graph,
                              char src_sentinel, char snk_sentinel);

// change the sentinel node characters
void reassign_sentinels(BaseGraph& graph, SentinelTableau& tableau,
                        char src_sentinel, char snk_sentinel);

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
        uint64_t new_path_id = appending.add_path(component.path_name(path_id));
        for (uint64_t node_id : component.path(path_id)) {
            appending.extend_path(new_path_id, prev_node_size + node_id);
        }
    }
}

}

#endif /* centrolign_modify_graph_hpp */
