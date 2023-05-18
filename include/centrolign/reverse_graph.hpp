#ifndef centrolign_reverse_graph_hpp
#define centrolign_reverse_graph_hpp

#include <vector>
#include <cstdint>

namespace centrolign {

/*
 * A topology-only overlay with the reverse topology
 */
template<class Graph>
class ReverseGraph {
public:
    
    ReverseGraph(const Graph& graph) noexcept;
    ReverseGraph() noexcept = default;
    ~ReverseGraph() noexcept = default;
    
    size_t node_size() const;
    const std::vector<uint64_t>& next(uint64_t node_id) const;
    const std::vector<uint64_t>& previous(uint64_t node_id) const;
    size_t next_size(uint64_t node_id) const;
    size_t previous_size(uint64_t node_id) const;
    
private:
    
    const Graph* graph = nullptr;
};

/*
 * Template implementations
 */

template<class Graph>
ReverseGraph<Graph>::ReverseGraph(const Graph& graph) noexcept : graph(&graph)  {
    // only need to initialize
}

template<class Graph>
size_t ReverseGraph<Graph>::node_size() const {
    return graph->node_size();
}

template<class Graph>
const std::vector<uint64_t>& ReverseGraph<Graph>::next(uint64_t node_id) const {
    return graph->previous(node_id);
}

template<class Graph>
const std::vector<uint64_t>& ReverseGraph<Graph>::previous(uint64_t node_id) const {
    return graph->next(node_id);
}

template<class Graph>
size_t ReverseGraph<Graph>::next_size(uint64_t node_id) const {
    return graph->previous_size(node_id);
}

template<class Graph>
size_t ReverseGraph<Graph>::previous_size(uint64_t node_id) const {
    return graph->next_size(node_id);
}

}

#endif /* centrolign_reverse_graph_hpp */
