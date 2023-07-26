#ifndef centrolign_source_sink_graph_hpp
#define centrolign_source_sink_graph_hpp

#include <vector>
#include <cstdint>
#include <unordered_set>

namespace centrolign {

/*
 * A topology-only overlay that adds a single source and single sink
 */
template<class Graph>
class SourceSinkGraph {
public:
    
    SourceSinkGraph(const Graph& graph) noexcept;
    SourceSinkGraph() noexcept = default;
    ~SourceSinkGraph() noexcept = default;
    
    size_t node_size() const;
    const std::vector<uint64_t>& next(uint64_t node_id) const;
    const std::vector<uint64_t>& previous(uint64_t node_id) const;
    size_t next_size(uint64_t node_id) const;
    size_t previous_size(uint64_t node_id) const;
    
    uint64_t source_id() const;
    uint64_t sink_id() const;
    
private:
    
    std::vector<uint64_t> source_nodes;
    std::vector<uint64_t> sink_nodes;
    std::unordered_set<uint64_t> source_node_set;
    std::unordered_set<uint64_t> sink_node_set;
    
    std::vector<uint64_t> source_edge;
    std::vector<uint64_t> sink_edge;
    const std::vector<uint64_t> null_edges;
    
    const Graph* graph = nullptr;
};

/*
 * Template implementations
 */

template<class Graph>
SourceSinkGraph<Graph>::SourceSinkGraph(const Graph& graph) noexcept :
    graph(&graph), source_edge(1, graph.node_size()), sink_edge(1, graph.node_size() + 1)
{
    for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
        if (graph.previous_size(node_id) == 0) {
            source_nodes.push_back(node_id);
            source_node_set.insert(node_id);
        }
        if (graph.next_size(node_id) == 0) {
            sink_nodes.push_back(node_id);
            sink_node_set.insert(node_id);
        }
    }
}

template<class Graph>
size_t SourceSinkGraph<Graph>::node_size() const {
    return graph->node_size() + 2;
}

template<class Graph>
const std::vector<uint64_t>& SourceSinkGraph<Graph>::next(uint64_t node_id) const {
    if (node_id == source_id()) {
        return source_nodes;
    }
    else if (node_id == sink_id()) {
        return null_edges;
    }
    else if (sink_node_set.count(node_id)) {
        return sink_edge;
    }
    else  {
        return graph->next(node_id);
    }
}

template<class Graph>
const std::vector<uint64_t>& SourceSinkGraph<Graph>::previous(uint64_t node_id) const {
    if (node_id == source_id()) {
        return null_edges;
    }
    else if (node_id == sink_id()) {
        return sink_nodes;
    }
    else if (source_node_set.count(node_id)) {
        return source_edge;
    }
    else  {
        return graph->previous(node_id);
    }
}

template<class Graph>
size_t SourceSinkGraph<Graph>::next_size(uint64_t node_id) const {
    if (node_id == source_id()) {
        return source_nodes.size();
    }
    else if (node_id == sink_id()) {
        return 0;
    }
    else if (sink_node_set.count(node_id)) {
        return 1;
    }
    else  {
        return graph->next_size(node_id);
    }
}

template<class Graph>
size_t SourceSinkGraph<Graph>::previous_size(uint64_t node_id) const {
    if (node_id == source_id()) {
        return 0;
    }
    else if (node_id == sink_id()) {
        return sink_nodes.size();
    }
    else if (source_node_set.count(node_id)) {
        return 1;
    }
    else  {
        return graph->previous_size(node_id);
    }
}

template<class Graph>
uint64_t SourceSinkGraph<Graph>::source_id() const {
    return graph->node_size();
}

template<class Graph>
uint64_t SourceSinkGraph<Graph>::sink_id() const {
    return graph->node_size() + 1;
}

}

#endif /* centrolign_source_sink_graph_hpp */
