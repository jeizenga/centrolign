#ifndef centrolign_path_merge_hpp
#define centrolign_path_merge_hpp

#include <vector>
#include <cstdint>
#include <iostream>
#include <tuple>

#include "centrolign/topological_order.hpp"
#include "centrolign/modify_graph.hpp"

namespace centrolign {


/*
 * Data structure for reachability using a path cover
 */
class PathMerge {
public:
    // build using only embedded paths
    template<class PGraph>
    PathMerge(const PGraph& graph);
    // build using embedded paths and sentinel nodes
    template<class PGraph>
    PathMerge(const PGraph& graph, const SentinelTableau& tableau);
    
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
    
    // generate edges to nodes that are nearest predecessors within one
    // of the chains
    std::vector<std::vector<uint64_t>> chain_forward_edges() const;
    
private:
    template<class PGraph>
    PathMerge(const PGraph& graph, const SentinelTableau* tableau);
    
    static const bool debug = false;
    
    // an arbitrary choice of a first path
    std::vector<uint64_t> path_head;
    // records of (index on path, next path)
    std::vector<std::vector<std::pair<size_t, uint64_t>>> index_on_path;
    // the last to reach table
    std::vector<std::vector<size_t>> table;
};

/*
 * Template implementations
 */

template<class PGraph>
PathMerge::PathMerge(const PGraph& graph) : PathMerge(graph, nullptr) {
    
}

template<class PGraph>
PathMerge::PathMerge(const PGraph& graph, const SentinelTableau& tableau) : PathMerge(graph, &tableau) {
    
}

template<class PGraph>
PathMerge::PathMerge(const PGraph& graph, const SentinelTableau* tableau) :
    path_head(graph.node_size(), -1),
    index_on_path(graph.node_size(),
                  std::vector<std::pair<size_t, uint64_t>>(graph.path_size() + int(tableau != nullptr),
                                                           std::pair<size_t, uint64_t>(-1, -1))),
    table(graph.node_size(), std::vector<size_t>(graph.path_size() + int(tableau != nullptr), -1))
{
    
    static const bool debug = false;
    
    // identify steps of paths and seed the DP
    for (uint64_t path_id = 0; path_id < graph.path_size(); ++path_id) {
        size_t index = 0;
        for (auto node_id : graph.path(path_id)) {
            for (auto next_id : graph.next(node_id)) {
                // note: we don't need the max here because we're going in increasing order
                table[next_id][path_id] = index;
            }
            index_on_path[node_id][path_id] = std::make_pair(index, path_head[node_id]);
            path_head[node_id] = path_id;
            ++index;
        }
    }
    
    // use DP to fill it in between cross-path edges
    for (auto node_id : topological_order(graph)) {
        auto& row = table[node_id];
        for (auto prev_id : graph.previous(node_id)) {
            const auto& prev_row = table[prev_id];
            for (uint64_t path_id = 0; path_id < graph.path_size(); ++path_id) {
                // cast to int64_t so that -1's get overwritten
                row[path_id] = std::max<int64_t>(prev_row[path_id], row[path_id]);
            }
        }
    }
    
    // add a "pseudo-path" chain for the sentinels
    if (tableau) {
        index_on_path[tableau->src_id][graph.path_size()].first = 0;
        index_on_path[tableau->snk_id][graph.path_size()].first = 1;
        path_head[tableau->src_id] = graph.path_size();
        path_head[tableau->snk_id] = graph.path_size();
        for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
            if (node_id != tableau->src_id) {
                table[node_id][graph.path_size()] = 0;
            }
        }
    }
    
    if (debug) {
        std::cerr << "path lists:\n";
        for (uint64_t n = 0; n < graph.node_size(); ++n) {
            std::cerr << n << ":\t" << path_head[n] << '\n';
        }
        
        std::cerr << "index on paths:\n";
        for (uint64_t n = 0; n < graph.node_size(); ++n) {
            std::cerr << n << ":";
            for (uint64_t p = 0; p < index_on_path[n].size(); ++p) {
                std::cerr << '\t' << (int64_t) index_on_path[n][p].first << ',' << (int64_t) index_on_path[n][p].second;
            }
            std::cerr << '\n';
        }
        
        std::cerr << "table:\n";
        for (uint64_t n = 0; n < graph.node_size(); ++n) {
            std::cerr << n << ":";
            for (uint64_t p = 0; p < table[n].size(); ++p) {
                std::cerr << '\t' << (int64_t) table[n][p];
            }
            std::cerr << '\n';
        }
    }
}

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

}

#endif /* centrolign_path_merge_hpp */
