#ifndef centrolign_gesa_hpp
#define centrolign_gesa_hpp

#include <vector>
#include <utility>
#include <algorithm>

#include "graph.hpp"
#include "append_component.hpp"
#include "determinize.hpp"
#include "path_graph.hpp"

namespace centrolign {

class GESA {
public:
    
    // note: these constructors modify the original graphs by adding sentinels
    template<class BGraph>
    GESA(const std::vector<const BGraph*>& graphs);
    template<class BGraph>
    GESA(const BGraph& graph);
    
    
private:
    
    template<class BGraph>
    GESA(const BGraph** const graphs, size_t num_graphs);
    
};

/*
 * Template implementations
 */

template<class BGraph>
GESA::GESA(const std::vector<const BGraph*>& graphs) : GESA(graphs.data(), graph.size()) {
    // nothing to do besides dispatch
}

template<class BGraph>
GESA::GESA(const BGraph& graph) : GESA(&(&graphs), 1) {
    // nothing to do besides dispatch
}

template<class BGraph>
GESA::GESA(const BGraph** const graphs, size_t num_graphs) {
    
    // find the maximum character among all the graphs
    // TODO: in some cases we'll already know this, in which case this feels
    // a bit wasteful
    char max_char = 0;
    for (size_t i = 0; i < num_graphs; ++i) {
        const BGraph& graph = *graphs[i];
        for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
            max_char = std::max(max_char, graph.base(node_id));
        }
    }
    
    // make a single graph with all of the components, with sentinels attached
    std::vector<std::pair<char, char>> component_sentinels(num_graphs);
    BaseGraph joined;
    for (size_t i = 0; i < num_graphs; ++i) {
        // add a new component
        const BGraph& graph = *graphs[i];
        uint64_t first_id = joined.node_size();
        append_component(joined, graph);
        
        // add the sentinels
        uint64_t src_id = joined.add_node(++max_char);
        component_sentinels[i].first = max_char;
        uint64_t snk_id = joined.add_node(++max_char);
        component_sentinels[i].second = max_char;
        
        // add edges to the sentinel nodes
        if (joined.node_size() == first_id + 2) {
            joined.add_edge(src_id, snk_id);
        }
        for (uint64_t node_id = first_id; node_id < src_id; ++node_id) {
            if (joined.next_size(node_id) == 0) {
                joined.add_edge(node_id, snk_id);
            }
            if (joined.previous_size(node_id) == 0) {
                joined.add_edge(src_id, node_id);
            }
        }
    }
    
    // convert the graph into an equivalent reverse deterministic graph
    joined = determinize<BaseGraph>(joined);
    
    // make into a path graph
    PathGraph path_graph(joined);
    // do prefix doubling until prefix sorted
    while (!joined.is_prefix_sorted()) {
        path_graph = PathGraph(path_graph);
    }
    // use the original grpah to add the edges in
    path_graph.construct_edges(joined);
    
    // reassign node IDs to be ordered by prefix rank
    path_graph.order_by_rank();
}

}

#endif /* centrolign_gesa_hpp */
