#ifndef centrolign_compacted_graph_hpp
#define centrolign_compacted_graph_hpp

#include <vector>
#include <cstdint>
#include <string>
#include <utility>
#include <unordered_map>
#include <iostream>

namespace centrolign {

/*
 * Expresses the topology and node lengths of a fully compacted
 * graph
 */
class CompactedGraph {
public:
    
    template<class Graph>
    CompactedGraph(const Graph& graph);
    CompactedGraph() = default;
    CompactedGraph(const CompactedGraph& other) noexcept = default;
    CompactedGraph(CompactedGraph&& other) noexcept = default;
    ~CompactedGraph() = default;
    
    CompactedGraph& operator=(const CompactedGraph& other) noexcept = default;
    CompactedGraph& operator=(CompactedGraph&& other) noexcept = default;
    
    size_t node_size() const;
    size_t label_size(uint64_t node_id) const;
    uint64_t front(uint64_t node_id) const;
    uint64_t back(uint64_t node_id) const;
    const std::vector<uint64_t>& next(uint64_t node_id) const;
    const std::vector<uint64_t>& previous(uint64_t node_id) const;
    size_t next_size(uint64_t node_id) const;
    size_t previous_size(uint64_t node_id) const;
    
private:
    
    struct Node {
        
        Node(uint64_t front, uint64_t back, size_t size) : front(front), back(back), size(size) { }
        Node() = default;
        ~Node() = default;
        
        uint64_t front;
        uint64_t back;
        size_t size;
        
        std::vector<uint64_t> next;
        std::vector<uint64_t> prev;
    };
    
    std::vector<Node> nodes;
};




/*
 * Template implementations
 */

template<class Graph>
CompactedGraph::CompactedGraph(const Graph& graph) {
    
    std::unordered_map<uint64_t, uint64_t> front_forward_trans;
    
    // add unipath nodes
    for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
        // TODO: should i handle sentinel nodes differently?
        if (graph.previous_size(node_id) != 1 ||
            (graph.previous_size(node_id) == 1 && graph.next_size(graph.previous(node_id).front()) != 1)) {
            size_t size = graph.label_size(node_id);
            uint64_t back = node_id;
            while (graph.next_size(back) == 1 && graph.previous_size(graph.next(back).front()) == 1) {
                back = graph.next(back).front();
                size += graph.label_size(back);
            }
            front_forward_trans[node_id] = nodes.size();
            nodes.emplace_back(node_id, back, size);
        }
    }
    
    // add edges
    for (uint64_t node_id = 0; node_id < node_size(); ++node_id) {
        auto& node = nodes[node_id];
        for (auto trans_next_id : graph.next(node.back)) {
            auto next_id = front_forward_trans[trans_next_id];
            node.next.push_back(next_id);
            nodes[next_id].prev.push_back(node_id);
        }
    }
}

}

#endif /* centrolign_compacted_graph_hpp */
