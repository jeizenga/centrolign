#ifndef centrolign_path_graph_hpp
#define centrolign_path_graph_hpp

#include <vector>
#include <cstdint>

namespace centrolign {

class PathGraph {
public:
    
    // copy from an initial base graph
    template<class BGraph>
    PathGraph(const BGraph& graph);
    // construct by a prefix-doubling step from another PathGraph
    PathGraph(const PathGraph& graph);
    
    PathGraph() = default;
    ~PathGraph() = default;
    
    size_t node_size() const;
    uint64_t from(uint64_t node_id) const;
    uint64_t to(uint64_t node_id) const;
    size_t rank(uint64_t node_id) const;
    
    bool is_prefix_sorted() const;
    
    // construct edges using the original graph you copied from
    template<class BGraph>
    void construct_edges(const BGraph& graph);
    
    // these functions are only valid after calling construct_edges
    std::vector<uint64_t> next(uint64_t node_id) const;
    std::vector<uint64_t> previous(uint64_t node_id) const;
    size_t next_size(uint64_t node_id) const;
    size_t previous_size(uint64_t node_id) const;
    
    
private:
    
    struct PathGraphNode;
    struct PathGraphEdges;
    
    // used to represent the "next node" after a sink
    static uint64_t null_id = -1;
    
    std::vector<PathGraphNode> nodes;
    std::vector<PathGraphEdges> edges;
    
    // LCP array over the shared prefixes of the unique ranks
    std::vector<size_t> lcp_array;
    
    size_t doubling_step = 0;
    
    struct PathGraphNode {
        
        PathGraphNode(uint64_t from, uint64_t to) : from(from), to(to) { }
        PathGraphNode(uint64_t from, uint64_t to,
                      size_t rank1, size_t rank2) : from(from), to(to)
                                                    rank(rank1), join_rank(rank2) { }
        PathGraphNode() = default;
        ~PathGraphNode() = default;
        
        uint64_t from = null_id;
        uint64_t to = null_id;
        size_t rank = 0;
        size_t join_rank = 0;
    };
    
    // these are held separately since the nodes frequently exist without edges
    struct PathGraphEdges {
        PathGraphEdges() = default;
        ~PathGraphEdges() = default;
        
        std::vector<uint64_t> next;
        std::vector<uint64_t> prev;
    };
};

/*
 * Template implementations
 */

template<class BGraph>
PathGraph::PathGraph(const BGraph& graph) {
    
    // TODO: probably over-engineering the initial sort
    std::vector<uint64_t> cumul;
    
    nodes.reserve(graph.node_size());
    for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
        if (graph.next_size(node_id) == 0) {
            nodes.emplace_back(node_id, null_id);
        }
        else {
            for (uint64_t next_id : graph.next(node_id)) {
                nodes.emplace_back(node_id, next_id);
            }
        }
        
        // record this base as 'seen'
        char base = graph.base(node_id);
        while (cumul.size() <= base) {
            cumul.emplace_back(0);
        }
        cumul[base] = 1;
    }
    
    // convert the 'seen' vector to a cumulative sum vector
    for (size_t i = 1; i < cumul.size(); ++i) {
        cumul[i] += cumul[i - 1];
    }
    // remove offset to convert to 0-based ranks
    for (size_t i = 0; i < cumul.size(); ++i) {
        // note: this goes to -1 for 0s, but we never will read them anyway
        --cumul[i];
    }
    
    // give rank to node based on its base
    for (PathGraphNode& node : nodes) {
        node.rank = cumul[graph.base(node.from)];
    }
}

template<class BGraph>
void PathGraph::construct_edges(const BGraph& graph) {
    // FIXME: implement
}

}

#endif /* centrolign_path_graph_hpp */
