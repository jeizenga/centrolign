#ifndef centrolign_path_merge_hpp
#define centrolign_path_merge_hpp

#include <vector>
#include <cstdint>
#include <iostream>
#include <tuple>

#include "centrolign/topological_order.hpp"
#include "centrolign/modify_graph.hpp"
#include "centrolign/graph.hpp"

namespace centrolign {

/*
 * Data structure for reachability using a path cover
 */
class PathMerge {
public:
    // build using only embedded paths
    PathMerge(const BaseGraph& graph);
    // build using embedded paths and sentinel nodes
    PathMerge(const BaseGraph& graph, const SentinelTableau& tableau);
    
    PathMerge() = default;
    ~PathMerge() = default;
    
    // number of chains
    inline size_t chain_size() const;
    
    // an arbitrary chain and index that the node belongs to
    inline std::pair<uint64_t, size_t> chain(uint64_t node_id) const;
    
    // the indexes (within their chain) of the nearest predecessors of this
    // node in each chain
    inline const std::vector<size_t>& predecessor_indexes(uint64_t node_id) const;
    
    // return true if both nodes are on a path and one can reach the other
    inline bool reachable(uint64_t from_id, uint64_t to_id) const;
    
    // the index of the node on the chain, or -1 if it is not on the path
    inline size_t index_on(uint64_t node_id, uint64_t path_id) const;
    
    // the chains that contain this node
    inline std::vector<uint64_t> chains_on(uint64_t node_id) const;
    
    // the node at this index of the chain
    inline uint64_t node_at(uint64_t chain_id, size_t index) const;
    
    // generate edges to nodes that are nearest predecessors within one
    // of the chains
    std::vector<std::vector<std::pair<uint64_t, uint64_t>>> chain_forward_edges() const;
    
    size_t memory_size() const;
    
private:
    PathMerge(const BaseGraph& graph, const SentinelTableau* tableau);
    
    static const bool debug = false;
    
    // an arbitrary choice of a first path
    std::vector<uint64_t> path_head;
    // records of (index on path, next path)
    std::vector<std::vector<std::pair<size_t, uint64_t>>> index_on_path;
    // the last to reach table
    std::vector<std::vector<size_t>> table;
    // the original graph (to avoid re-instantiating the paths)
    const BaseGraph* g = nullptr;
    // sentinels (to avoid reidentifying it for forward edges) TODO: ugly
    uint64_t src_id = -1;
    uint64_t snk_id = -1;
};

/*
 * Template implementations
 */


inline size_t PathMerge::chain_size() const {
    return table.empty() ? 0 : table.front().size();
}

inline std::pair<uint64_t, size_t> PathMerge::chain(uint64_t node_id) const {
    auto path_id = path_head[node_id];
    if (path_id != -1) {
        return std::make_pair(path_id, index_on_path[node_id][path_id].first);
    }
    return std::pair<uint64_t, size_t>(-1, -1);
}

inline const std::vector<size_t>& PathMerge::predecessor_indexes(uint64_t node_id) const {
    return table[node_id];
}

inline bool PathMerge::reachable(uint64_t from_id, uint64_t to_id) const {
    uint64_t chain_from;
    size_t idx_from;
    std::tie(chain_from, idx_from) = chain(from_id);
    if (chain_from == -1) {
        // not on a path
        // TODO: do i need to maintain this anymore?
        return false;
    }
    // is the last index that can reach the target after the source?
    size_t last_idx_to_reach = table[to_id][chain_from];
    return (last_idx_to_reach != -1 && idx_from <= last_idx_to_reach);
}

inline size_t PathMerge::index_on(uint64_t node_id, uint64_t path_id) const {
    return index_on_path[node_id][path_id].first;
}

inline std::vector<uint64_t> PathMerge::chains_on(uint64_t node_id) const {
    std::vector<uint64_t> paths;
    auto p = path_head[node_id];
    while (p != -1) {
        paths.push_back(p);
        p = index_on_path[node_id][p].second;
    }
    return paths;
}

inline uint64_t PathMerge::node_at(uint64_t chain_id, size_t index) const {
    if (chain_id == g->path_size()) {
        return index ? snk_id : src_id;
    }
    else {
        return g->path(chain_id)[index];
    }
}

}

#endif /* centrolign_path_merge_hpp */
