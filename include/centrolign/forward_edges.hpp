#ifndef centrolign_forward_edges_hpp
#define centrolign_forward_edges_hpp


namespace centrolign {

/*
 * The forward edges for sparse dynamic programming as defined by Makinen, et al (2019)
 */
template<typename UIntNode, typename UIntPath>
class ForwardEdges {
public:
    
    template<class XMerge>
    ForwardEdges(const XMerge& xmerge, const std::vector<bool>* to_nodes = nullptr,
                 const std::vector<bool>* from_nodes = nullptr);
    ForwardEdges() = default;
    ~ForwardEdges() = default;
    
    inline const std::vector<std::pair<UIntNode, UIntPath>>& edges(uint64_t node_id) const;
    
    size_t memory_size() const;
    
private:
    
    std::vector<std::vector<std::pair<UIntNode, UIntPath>>> fwd_edges;
};


/*
 * Template implementations
 */

template<typename UIntNode, typename UIntPath>
template<class XMerge>
ForwardEdges<UIntNode, UIntPath>::ForwardEdges(const XMerge& xmerge, const std::vector<bool>* to_nodes, const std::vector<bool>* from_nodes) : fwd_edges(xmerge.node_size()) {
    
    for (uint64_t node_id = 0; node_id < fwd_edges.size(); ++node_id) {
        if (to_nodes && !(*to_nodes)[node_id]) {
            continue;
        }
        for (uint64_t p = 0; p < xmerge.chain_size(); ++p) {
            auto idx = xmerge.predecessor_index(node_id, p);
            if (idx != std::numeric_limits<decltype(idx)>::max()) {
                auto from_id = xmerge.node_at(p, idx);
                if (!from_nodes || (*from_nodes)[from_id]) {
                    fwd_edges[from_id].emplace_back(node_id, p);
                }
            }
        }
    }
    
    for (auto& edge_list : fwd_edges) {
        edge_list.shrink_to_fit();
    }
}


template<typename UIntNode, typename UIntPath>
inline const std::vector<std::pair<UIntNode, UIntPath>>& ForwardEdges<UIntNode, UIntPath>::edges(uint64_t node_id) const {
    return fwd_edges[node_id];
}

template<typename UIntNode, typename UIntPath>
size_t ForwardEdges<UIntNode, UIntPath>::memory_size() const {
    size_t mem = sizeof(fwd_edges) + fwd_edges.capacity() * sizeof(typename decltype(fwd_edges)::value_type);
    for (const auto& edge_list : fwd_edges) {
        mem += edge_list.capacity() * sizeof(typename decltype(fwd_edges)::value_type::value_type);
    }
    return mem;
}

}

#endif /* centrolign_forward_edges_hpp */
