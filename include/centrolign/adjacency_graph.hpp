#ifndef centrolign_adjacency_graph_hpp
#define centrolign_adjacency_graph_hpp

#include <vector>
#include <cstdint>
#include <string>
#include <utility>
#include <unordered_map>

namespace centrolign {

/*
 * An directed graph with nodes corresponding to adjacency components and
 * edges corresponding to nodes in the underlying graph
 */
class AdjacencyGraph {
public:
    
    template<class Graph>
    AdjacencyGraph(const Graph& graph);
    AdjacencyGraph() = default;
    ~AdjacencyGraph() = default;
    
    struct Edge {
        Edge(uint64_t target, uint64_t label) : target(target), label(label) { }
        Edge() = default;
        ~Edge() = default;
        
        // other node ID pointed to by this edge
        uint64_t target = -1;
        // the node ID of the node in the underlying graph that corresponds to this edge
        uint64_t label = -1;
    };
    
    inline size_t node_size() const;
    inline const std::vector<Edge>& next_edges(uint64_t node_id) const;
    inline size_t next_size(uint64_t node_id) const;
    inline const std::vector<Edge>& previous_edges(uint64_t node_id) const;
    inline size_t previous_size(uint64_t node_id) const;
    
    // to conform the with the standard graph interface in this code base
    inline std::vector<uint64_t> next(uint64_t node_id) const;
    inline std::vector<uint64_t> previous(uint64_t node_id) const;
private:
    
    std::vector<std::pair<std::vector<Edge>, std::vector<Edge>>> adj_list;
};




/*
 * Template and inline implementations
 */

template<class Graph>
AdjacencyGraph::AdjacencyGraph(const Graph& graph) {
    
    // keep track of which adjacency components each of the two sides of each node is in
    std::vector<uint64_t> adjacency_component(2 * graph.node_size(), -1);
    
    // identify the adjacency components
    size_t next_adj_comp = 0;
    for (uint64_t i = 0; i < adjacency_component.size(); ++i) {
        
        if (adjacency_component[i] != -1) {
            // we've already found this node side's adjacency component
            continue;
        }
        
        // bounce back and forth between forward and reverse edges from the seed
        adjacency_component[i] = next_adj_comp;
        std::vector<std::pair<uint64_t, bool>> stack(1, std::pair<uint64_t, bool>(i / 2, i % 2));
        while (!stack.empty()) {
            
            uint64_t node_id;
            bool left;
            std::tie(node_id, left) = stack.back();
            stack.pop_back();
            
            for (auto next_id : (left ? graph.previous(node_id) : graph.next(node_id))) {
                
                size_t j = 2 * next_id + !left;
                if (adjacency_component[j] == -1) {
                    // we haven't seen this node side yet
                    adjacency_component[j] = next_adj_comp;
                    stack.emplace_back(next_id, !left);
                }
            }
        }
        
        ++next_adj_comp;
    }
    
    adj_list.resize(next_adj_comp);
    
    // create edges to the underlying graph's nodes
    for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
        
        uint64_t right_adj_comp = adjacency_component[2 * node_id];
        uint64_t left_adj_comp = adjacency_component[2 * node_id + 1];
        
        adj_list[left_adj_comp].second.emplace_back(right_adj_comp, node_id);
        adj_list[right_adj_comp].first.emplace_back(left_adj_comp, node_id);
        
    }
}

inline size_t AdjacencyGraph::node_size() const {
    return adj_list.size();
}
inline const std::vector<AdjacencyGraph::Edge>& AdjacencyGraph::next_edges(uint64_t node_id) const {
    return adj_list[node_id].second;
}
inline size_t AdjacencyGraph::next_size(uint64_t node_id) const {
    return adj_list[node_id].second.size();
}
inline const std::vector<AdjacencyGraph::Edge>& AdjacencyGraph::previous_edges(uint64_t node_id) const {
    return adj_list[node_id].first;
}
inline size_t AdjacencyGraph::previous_size(uint64_t node_id) const {
    return adj_list[node_id].first.size();
}
inline std::vector<uint64_t> AdjacencyGraph::next(uint64_t node_id) const {
    const auto& edges = next_edges(node_id);
    std::vector<uint64_t> targets(edges.size());
    for (size_t i = 0; i < edges.size(); ++i) {
        targets[i] = edges[i].target;
    }
    return targets;
}
inline std::vector<uint64_t> AdjacencyGraph::previous(uint64_t node_id) const {
    const auto& edges = previous_edges(node_id);
    std::vector<uint64_t> targets(edges.size());
    for (size_t i = 0; i < edges.size(); ++i) {
        targets[i] = edges[i].target;
    }
    return targets;
}
}

#endif /* centrolign_adjacency_graph_hpp */
