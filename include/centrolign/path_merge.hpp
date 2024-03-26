#ifndef centrolign_path_merge_hpp
#define centrolign_path_merge_hpp

#include <vector>
#include <cstdint>
#include <iostream>
#include <tuple>
#include <limits>

#include "centrolign/topological_order.hpp"
#include "centrolign/modify_graph.hpp"
#include "centrolign/graph.hpp"

namespace centrolign {

/*
 * Data structure for reachability using a path cover
 */
template<typename UIntSize = size_t, typename UIntChain = uint64_t>
class PathMerge {
public:
    // build using only embedded paths
    PathMerge(const BaseGraph& graph);
    // build using embedded paths and sentinel nodes
    PathMerge(const BaseGraph& graph, const SentinelTableau& tableau);
    
    PathMerge() = default;
    PathMerge(const PathMerge& other) = default;
    PathMerge(PathMerge&& other) = default;
    ~PathMerge() = default;
    PathMerge& operator=(const PathMerge& other) = default;
    PathMerge& operator=(PathMerge&& other) = default;
    
    // number of chains
    inline size_t chain_size() const;
    
    // an arbitrary chain and index that the node belongs to
    inline std::pair<uint64_t, size_t> chain(uint64_t node_id) const;
    
    // the indexes (within their chain) of the nearest predecessors of this
    // node in each chain
    inline const std::vector<UIntSize>& predecessor_indexes(uint64_t node_id) const;
    
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
        
    // an arbitrary choice of a first path for each node
    std::vector<UIntChain> path_head;
    // records of (index on path, next path containing the node)
    std::vector<std::vector<std::pair<UIntSize, UIntChain>>> index_on_path;
    // the last to reach table
    std::vector<std::vector<UIntSize>> table;
    // the original graph (to avoid re-instantiating the paths)
    const BaseGraph* g = nullptr;
    // sentinels (to avoid reidentifying it for forward edges) TODO: ugly
    uint64_t src_id = -1;
    uint64_t snk_id = -1;
};

/*
 * Template implementations
 */


template<typename UIntSize, typename UIntChain>
PathMerge<UIntSize, UIntChain>::PathMerge(const BaseGraph& graph) : PathMerge(graph, nullptr) {
    
}

template<typename UIntSize, typename UIntChain>
PathMerge<UIntSize, UIntChain>::PathMerge(const BaseGraph& graph, const SentinelTableau& tableau) : PathMerge(graph, &tableau) {
    src_id = tableau.src_id;
    snk_id = tableau.snk_id;
}

template<typename UIntSize, typename UIntChain>
PathMerge<UIntSize, UIntChain>::PathMerge(const BaseGraph& graph, const SentinelTableau* tableau) :
    path_head(graph.node_size(), std::numeric_limits<UIntChain>::max()),
    index_on_path(graph.path_size() + int(tableau != nullptr),
                  std::vector<std::pair<UIntSize, UIntChain>>(graph.node_size(),
                                                              std::make_pair(std::numeric_limits<UIntSize>::max(),
                                                                             std::numeric_limits<UIntChain>::max()))),
    table(graph.node_size(), std::vector<UIntSize>(graph.path_size() + int(tableau != nullptr),
                                                   std::numeric_limits<UIntSize>::max())),
    g(&graph)
{
    
    static const bool debug = false;
    
    if (debug) {
        std::cerr << "recording indexes on paths\n";
    }
        
    // identify steps of paths and seed the DP
    for (UIntChain path_id = 0; path_id < graph.path_size(); ++path_id) {
        UIntSize index = 0;
        for (auto node_id : graph.path(path_id)) {
            for (auto next_id : graph.next(node_id)) {
                // note: we don't need the max here because we're going in increasing order,
                // so the max is the last to write
                table[next_id][path_id] = index;
            }
            index_on_path[path_id][node_id] = std::make_pair(index, path_head[node_id]);
            path_head[node_id] = path_id;
            ++index;
        }
    }
    
    if (debug) {
        std::cerr << "filling the last to reach table\n";
    }
    
    // use DP to fill it in between cross-path edges
    for (auto node_id : topological_order(graph)) {
        auto& row = table[node_id];
        for (auto prev_id : graph.previous(node_id)) {
            const auto& prev_row = table[prev_id];
            for (UIntChain path_id = 0; path_id < graph.path_size(); ++path_id) {
                if (row[path_id] == std::numeric_limits<UIntSize>::max()) {
                    row[path_id] = prev_row[path_id];
                }
                else if (prev_row[path_id] != std::numeric_limits<UIntSize>::max()) {
                    row[path_id] = std::max(prev_row[path_id], row[path_id]);
                }
            }
        }
    }
    
    // add a "pseudo-path" chain for the sentinels
    if (tableau) {
        if (debug) {
            std::cerr << "creating the tableau pseudo-path\n";
        }
        index_on_path[graph.path_size()][tableau->src_id].first = 0;
        index_on_path[graph.path_size()][tableau->snk_id].first = 1;
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
            std::cerr << n << ":\t" << (int64_t) path_head[n] << '\n';
        }
        
        std::cerr << "index on paths:\n";
        if (!index_on_path.empty()) {
            for (uint64_t n = 0; n < index_on_path[0].size(); ++n) {
                std::cerr << n << ":";
                for (uint64_t p = 0; p < index_on_path.size(); ++p) {
                    std::cerr << '\t' << (int64_t) index_on_path[p][n].first << ',' << (int64_t) index_on_path[p][n].second;
                }
                std::cerr << '\n';
            }
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

template<typename UIntSize, typename UIntChain>
size_t PathMerge<UIntSize, UIntChain>::memory_size() const {
    size_t index_size = 0;
    size_t table_size = 0;
    for (size_t i = 0; i < index_on_path.size(); ++i) {
        index_size += index_on_path[i].capacity();
    }
    for (size_t i = 0; i < table.size(); ++i) {
        table_size += table[i].capacity();
    }
    
    return (index_size * sizeof(typename decltype(index_on_path)::value_type::value_type)
            + table_size * sizeof(typename decltype(table)::value_type::value_type)
            + path_head.capacity() * sizeof(typename decltype(path_head)::value_type)
            + index_on_path.capacity() * sizeof(typename decltype(index_on_path)::value_type)
            + table.capacity() * sizeof(typename decltype(table)::value_type)
            + sizeof(index_on_path) + sizeof(path_head) + sizeof(table));
}

template<typename UIntSize, typename UIntChain>
std::vector<std::vector<std::pair<uint64_t, uint64_t>>> PathMerge<UIntSize, UIntChain>::chain_forward_edges() const {
    
    // identify the forward links
    std::vector<std::vector<std::pair<uint64_t, uint64_t>>> forward(table.size());
    for (uint64_t node_id = 0; node_id < table.size(); ++node_id) {
        const auto& row = table[node_id];
        for (uint64_t p = 0; p < g->path_size(); ++p) {
            if (row[p] != std::numeric_limits<UIntSize>::max()) {
                forward[g->path(p)[row[p]]].emplace_back(node_id, p);
            }
        }
        // TODO: do i really even need forward links on the pseudo-path?
        // i could get rid of the saved src_id if i didn't make these
        if (src_id != -1 && node_id != src_id) {
            forward[src_id].emplace_back(node_id, g->path_size());
        }
    }
    
    return forward;
}

template<typename UIntSize, typename UIntChain>
inline size_t PathMerge<UIntSize, UIntChain>::chain_size() const {
    return table.empty() ? 0 : table.front().size();
}

template<typename UIntSize, typename UIntChain>
inline std::pair<uint64_t, size_t> PathMerge<UIntSize, UIntChain>::chain(uint64_t node_id) const {
    auto path_id = path_head[node_id];
    if (path_id != std::numeric_limits<UIntChain>::max()) {
        return std::pair<uint64_t, size_t>(path_id, index_on_path[path_id][node_id].first);
    }
    return std::pair<uint64_t, size_t>(-1, -1);
}

template<typename UIntSize, typename UIntChain>
inline const std::vector<UIntSize>& PathMerge<UIntSize, UIntChain>::predecessor_indexes(uint64_t node_id) const {
    return table[node_id];
}

template<typename UIntSize, typename UIntChain>
inline bool PathMerge<UIntSize, UIntChain>::reachable(uint64_t from_id, uint64_t to_id) const {
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
    return (last_idx_to_reach != std::numeric_limits<UIntSize>::max() && idx_from <= last_idx_to_reach);
}

template<typename UIntSize, typename UIntChain>
inline size_t PathMerge<UIntSize, UIntChain>::index_on(uint64_t node_id, uint64_t path_id) const {
    return index_on_path[path_id][node_id].first;
}

template<typename UIntSize, typename UIntChain>
inline std::vector<uint64_t> PathMerge<UIntSize, UIntChain>::chains_on(uint64_t node_id) const {
    std::vector<uint64_t> paths;
    auto p = path_head[node_id];
    while (p != std::numeric_limits<UIntChain>::max()) {
        paths.push_back(p);
        p = index_on_path[p][node_id].second;
    }
    return paths;
}

template<typename UIntSize, typename UIntChain>
inline uint64_t PathMerge<UIntSize, UIntChain>::node_at(uint64_t chain_id, size_t index) const {
    if (chain_id == g->path_size()) {
        return index ? snk_id : src_id;
    }
    else {
        return g->path(chain_id)[index];
    }
}

}

#endif /* centrolign_path_merge_hpp */
