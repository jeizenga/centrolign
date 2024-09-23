#ifndef centrolign_cactus_hpp
#define centrolign_cactus_hpp

#include <vector>
#include <cstdint>
#include <string>
#include <utility>
#include <unordered_map>

#include "centrolign/compacted_graph.hpp"
#include "centrolign/adjacency_graph.hpp"
#include "centrolign/utility.hpp"

namespace centrolign {

/**
 * The cactus graph for a sequence graph
 */
template<class BGraph>
class CactusGraph {
public:
    
    CactusGraph(const BGraph& graph);
    CactusGraph() = default;
    ~CactusGraph() = default;
    
    // get the sequence of nodes in the original underlying graph
    inline std::vector<uint64_t> next_edge_label(uint64_t node_id, size_t edge_idx) const;
    inline std::vector<uint64_t> previous_edge_label(uint64_t node_id, size_t edge_idx) const;
    // get the length of nodes in the original underlying graph
    inline size_t next_edge_label_size(uint64_t node_id, size_t edge_idx) const;
    inline size_t previous_edge_label_size(uint64_t node_id, size_t edge_idx) const;
    
    inline size_t node_size() const;
    inline const std::vector<uint64_t>& next(uint64_t node_id) const;
    inline const std::vector<uint64_t>& previous(uint64_t node_id) const;
    inline size_t next_size(uint64_t node_id) const;
    inline size_t previous_size(uint64_t node_id) const;
    
private:
    
    struct Node {
        Node() = default;
        ~Node() = default;
        std::vector<uint64_t> previous;
        std::vector<std::pair<uint64_t, size_t>> previous_origin;
        std::vector<uint64_t> next;
        std::vector<std::pair<uint64_t, size_t>> next_origin;
     };
    
    inline uint64_t edge_to_compacted_id(uint64_t node_id, bool next, size_t edge_idx) const;
    
    inline std::vector<uint64_t> walk_compacted_node(uint64_t compacted_id) const;
    
    CompactedGraph compacted_graph;
    
    AdjacencyGraph adjacency_graph;
    
    const BGraph* underlying_graph = nullptr;
    
    std::vector<Node> nodes;
};



/**
 * Template and inline implementations
 */

template<class BGraph>
CactusGraph<BGraph>::CactusGraph(const BGraph& graph) : compacted_graph(graph), adjacency_graph(compacted_graph), underlying_graph(&graph) {
    // note: graph construction in the initializers
    
    static const bool debug = false;
    
    auto three_edge_comps = three_edge_connected_components(adjacency_graph);
    
    std::vector<size_t> node_to_comp(adjacency_graph.node_size());
    for (size_t i = 0; i < three_edge_comps.size(); ++i) {
        for (auto node_id : three_edge_comps[i]) {
            node_to_comp[node_id] = i;
        }
    }
    
    nodes.resize(three_edge_comps.size());
    for (uint64_t node_id = 0; node_id < adjacency_graph.node_size(); ++node_id) {
        uint64_t comp1 = node_to_comp[node_id];
        const auto& next_edges = adjacency_graph.next_edges(node_id);
        for (size_t i = 0; i < next_edges.size(); ++i) {
            uint64_t comp2 = node_to_comp[next_edges[i].target];
            nodes[comp1].next.emplace_back(comp2);
            nodes[comp1].next_origin.emplace_back(node_id, i);
            nodes[comp2].previous.emplace_back(comp1);
            nodes[comp2].previous_origin.emplace_back(node_id, i);
        }
    }
    if (debug) {
        std::cerr << "compacted graph\n";
        print_topology(compacted_graph, std::cerr);
        std::cerr << "adjacency graph\n";
        print_topology(adjacency_graph, std::cerr);
        std::cerr << "cactus graph\n";
        print_topology(*this, std::cerr);
    }
}

template<class BGraph>
inline uint64_t CactusGraph<BGraph>::edge_to_compacted_id(uint64_t node_id, bool next, size_t edge_idx) const {
    const auto& origin = next ? nodes[node_id].next_origin : nodes[node_id].previous_origin;
    const auto& adj_edge = origin[edge_idx];
    return adjacency_graph.next_edges(adj_edge.first)[adj_edge.second].label;
}

template<class BGraph>
inline std::vector<uint64_t> CactusGraph<BGraph>::previous_edge_label(uint64_t node_id, size_t edge_idx) const {
    return walk_compacted_node(edge_to_compacted_id(node_id, false, edge_idx));
}

template<class BGraph>
inline size_t CactusGraph<BGraph>::previous_edge_label_size(uint64_t node_id, size_t edge_idx) const {
    return compacted_graph.label_size(edge_to_compacted_id(node_id, false, edge_idx));
}

template<class BGraph>
inline std::vector<uint64_t> CactusGraph<BGraph>::next_edge_label(uint64_t node_id, size_t edge_idx) const {
    return walk_compacted_node(edge_to_compacted_id(node_id, true, edge_idx));
}

template<class BGraph>
inline size_t CactusGraph<BGraph>::next_edge_label_size(uint64_t node_id, size_t edge_idx) const {
    return compacted_graph.label_size(edge_to_compacted_id(node_id, true, edge_idx));
}

template<class BGraph>
inline std::vector<uint64_t> CactusGraph<BGraph>::walk_compacted_node(uint64_t compacted_id) const {
    std::vector<uint64_t> path;
    path.reserve(compacted_graph.label_size(compacted_id));
    path.emplace_back(compacted_graph.front(compacted_id));
    while (path.back() != compacted_graph.back(compacted_id)) {
        path.push_back(underlying_graph->next(path.back()).front());
    }
    return path;
}

template<class BGraph>
inline size_t CactusGraph<BGraph>::node_size() const {
    return nodes.size();
}

template<class BGraph>
inline const std::vector<uint64_t>& CactusGraph<BGraph>::next(uint64_t node_id) const {
    return nodes[node_id].next;
}

template<class BGraph>
inline const std::vector<uint64_t>& CactusGraph<BGraph>::previous(uint64_t node_id) const {
    return nodes[node_id].previous;
}

template<class BGraph>
inline size_t CactusGraph<BGraph>::next_size(uint64_t node_id) const {
    return nodes[node_id].next.size();
}

template<class BGraph>
inline size_t CactusGraph<BGraph>::previous_size(uint64_t node_id) const {
    return nodes[node_id].previous.size();
}



}

#endif /* centrolign_cactus_hpp */
