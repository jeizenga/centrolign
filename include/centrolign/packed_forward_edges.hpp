#ifndef centrolign_packed_forward_edges_hpp
#define centrolign_packed_forward_edges_hpp


namespace centrolign {

/*
 * The forward edges for sparse dynamic programming as defined by Makinen, et al (2019)
 * backed by bit-packed data structures
 */
class PackedForwardEdges {
public:
    
    template<class XMerge>
    PackedForwardEdges(const XMerge& xmerge, const std::vector<bool>* to_nodes = nullptr,
                 const std::vector<bool>* from_nodes = nullptr);
    PackedForwardEdges() = default;
    ~PackedForwardEdges() = default;
    
    inline std::vector<std::pair<uint64_t, uint64_t>> edges(uint64_t node_id) const;
    
    inline size_t memory_size() const;
    
private:
    
    PagedVector<128, 127> node_interval;
    PackedVector chain_ids;
    PackedVector node_ids;
};


/*
 * Template implementations
 */

template<class XMerge>
PackedForwardEdges::PackedForwardEdges(const XMerge& xmerge, const std::vector<bool>* to_nodes, const std::vector<bool>* from_nodes) : node_interval(xmerge.node_size() + 1) {
    
    static const bool debug = false;
    
    // figure out the inter
    PackedVector num_edges(xmerge.node_size());
    for (uint64_t node_id = 0; node_id < xmerge.node_size(); ++node_id) {
        if (to_nodes && !(*to_nodes)[node_id]) {
            continue;
        }
        for (uint64_t p = 0; p < xmerge.chain_size(); ++p) {
            auto idx = xmerge.predecessor_index(node_id, p);
            if (idx != std::numeric_limits<decltype(idx)>::max()) {
                auto from_id = xmerge.node_at(p, idx);
                if (!from_nodes || (*from_nodes)[from_id]) {
                    ++num_edges[from_id];
                }
            }
        }
    }
    
    if (debug) {
        std::cerr << "num edges:\n";
        for (size_t i = 0; i < num_edges.size(); ++i) {
            std::cerr << i << '\t' << num_edges.at(i) << '\n';
        }
    }
    
    size_t offset = 0;
    for (size_t i = 0; i < xmerge.node_size(); ++i) {
        offset += size_t(num_edges[i]);
        node_interval[i + 1] = offset;
    }
    {
        // clear the data
        auto dummy = std::move(num_edges);
    }
    
    if (debug) {
        std::cerr << "node interval:\n";
        for (size_t i = 0; i < node_interval.size(); ++i) {
            std::cerr << i << '\t' << node_interval.at(i) << '\n';
        }
    }
    
    chain_ids = std::move(decltype(chain_ids)(offset));
    node_ids = std::move(decltype(node_ids)(offset));
    
    decltype(node_interval) next_idx(node_interval);
    
    for (uint64_t node_id = 0; node_id < xmerge.node_size(); ++node_id) {
        if (to_nodes && !(*to_nodes)[node_id]) {
            continue;
        }
        for (uint64_t p = 0; p < xmerge.chain_size(); ++p) {
            auto idx = xmerge.predecessor_index(node_id, p);
            if (idx != std::numeric_limits<decltype(idx)>::max()) {
                auto from_id = xmerge.node_at(p, idx);
                size_t i = next_idx[from_id];
                next_idx[from_id] = i + 1;
                if (!from_nodes || (*from_nodes)[from_id]) {
                    chain_ids[i] = p;
                    node_ids[i] = node_id;
                }
            }
        }
    }
    
    if (debug) {
        std::cerr << "edges:\n";
        for (size_t i = 0; i < chain_ids.size(); ++i) {
            std::cerr << i << '\t' << chain_ids.at(i) << '\t' << node_ids.at(i) << '\n';
        }
    }
}


inline std::vector<std::pair<uint64_t, uint64_t>> PackedForwardEdges::edges(uint64_t node_id) const {
    std::vector<std::pair<uint64_t, uint64_t>> fwd_edges;
    for (size_t i = node_interval[node_id], n = node_interval[node_id + 1]; i < n; ++i) {
        fwd_edges.emplace_back(node_ids[i], chain_ids[i]);
    }
    return fwd_edges;
}

inline size_t PackedForwardEdges::memory_size() const {
    return node_ids.memory_size() + chain_ids.memory_size() + node_interval.memory_size();
}


}

#endif /* centrolign_packed_forward_edges_hpp */
