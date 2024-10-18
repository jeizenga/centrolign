#ifndef centrolign_chain_cycle_graph_hpp
#define centrolign_chain_cycle_graph_hpp

#include <vector>
#include <cstdint>
#include <string>
#include <utility>
#include <unordered_map>

namespace centrolign {

/*
 * Graph that joins the source and sink vertex by an edge, turning the top-
 * level chain into a cycle
 */
template<class Graph>
class ChainCycleGraph {
public:
    
    ChainCycleGraph(const Graph& graph, uint64_t src_id, uint64_t snk_id);
    ChainCycleGraph() = default;
    ~ChainCycleGraph() = default;
    
    size_t node_size() const;
    std::vector<uint64_t> next(uint64_t node_id) const;
    std::vector<uint64_t> previous(uint64_t node_id) const;
    size_t next_size(uint64_t node_id) const;
    size_t previous_size(uint64_t node_id) const;
    
private:
    
    const Graph* graph = nullptr;
    
    uint64_t src_id = -1;
    uint64_t snk_id = -1;
};




/*
 * Template implementations
 */

template<class Graph>
ChainCycleGraph<Graph>::ChainCycleGraph(const Graph& graph, uint64_t src_id, uint64_t snk_id) : graph(&graph), src_id(src_id), snk_id(snk_id) {
    
}

template<class Graph>
size_t ChainCycleGraph<Graph>::node_size() const {
    return graph->node_size();
}

template<class Graph>
std::vector<uint64_t> ChainCycleGraph<Graph>::next(uint64_t node_id) const {
    std::vector<uint64_t> edges;
    if (node_id == snk_id) {
        edges.push_back(src_id);
    }
    else {
        edges = graph->next(node_id);
    }
    return edges;
}

template<class Graph>
std::vector<uint64_t> ChainCycleGraph<Graph>::previous(uint64_t node_id) const {
    std::vector<uint64_t> edges;
    if (node_id == src_id) {
        edges.push_back(snk_id);
    }
    else {
        edges = graph->previous(node_id);
    }
    return edges;
}

template<class Graph>
size_t ChainCycleGraph<Graph>::next_size(uint64_t node_id) const {
    return node_id == snk_id ? 1 : graph->next_size(node_id);
}

template<class Graph>
size_t ChainCycleGraph<Graph>::previous_size(uint64_t node_id) const {
    return node_id == src_id ? 1 : graph->previous_size(node_id);
}


}

#endif /* centrolign_chain_cycle_graph_hpp */
