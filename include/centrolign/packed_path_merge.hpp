#ifndef centrolign_packed_path_merge_hpp
#define centrolign_packed_path_merge_hpp

#include <vector>
#include <cstdint>
#include <iostream>
#include <tuple>
#include <limits>

#include "centrolign/topological_order.hpp"
#include "centrolign/modify_graph.hpp"
#include "centrolign/graph.hpp"
#include "centrolign/packed_vector.hpp"
#include "centrolign/paged_vector.hpp"

namespace centrolign {

/*
 * Data structure for reachability using a path cover backed by
 * bit-packed vectors
 */
template<typename UIntSize = size_t, typename UIntChain = uint64_t, size_t PageSize = 2048, size_t TiltBias = 128>
class PackedPathMerge {
public:
    // build using only embedded paths
    PackedPathMerge(const BaseGraph& graph);
    // build using embedded paths and sentinel nodes
    PackedPathMerge(const BaseGraph& graph, const SentinelTableau& tableau);
    
    PackedPathMerge() = default;
    PackedPathMerge(const PackedPathMerge& other) = default;
    PackedPathMerge(PackedPathMerge&& other) = default;
    ~PackedPathMerge() = default;
    PackedPathMerge& operator=(const PackedPathMerge& other) = default;
    PackedPathMerge& operator=(PackedPathMerge&& other) = default;
    
    using node_id_t = UIntSize;
    using chain_id_t = UIntChain;
    
    // number of chains
    inline size_t chain_size() const;
    
    // number of nodes
    inline size_t node_size() const;
    
    // an arbitrary chain and index that the node belongs to
    inline std::pair<uint64_t, size_t> chain(uint64_t node_id) const;
    
    // the indexes (within their chain) of the nearest predecessors of this
    // node in each chain
    inline size_t predecessor_index(uint64_t node_id, uint64_t path_id) const;
    
    // return true if both nodes are on a path and one can reach the other
    inline bool reachable(uint64_t from_id, uint64_t to_id) const;
    
    // the index of the node on the chain, or -1 if it is not on the path
    inline size_t index_on(uint64_t node_id, uint64_t path_id) const;
    
    // the chains that contain this node
    inline std::vector<uint64_t> chains_on(uint64_t node_id) const;
    
    // the node at this index of the chain
    inline uint64_t node_at(uint64_t chain_id, size_t index) const;
    
    size_t memory_size() const;
    
private:
    PackedPathMerge(const BaseGraph& graph, const SentinelTableau* tableau);
        
    // an arbitrary choice of a first path for each node
    PackedVector path_head;
    // for each path, the index of the node on that path + 1, or 0 if none
    std::vector<PagedVector<PageSize, TiltBias>> index_on_path;
    // for each path, the next ID + 1 of another path containing that node, or 0 if none
    std::vector<PackedVector> next_path_containing;
    // the last to reach table
    std::vector<PagedVector<PageSize, TiltBias>> table;
    // the original graph (to avoid re-instantiating the paths)
    const BaseGraph* g = nullptr;
    // sentinels (to avoid reidentifying it for forward edges) TODO: ugly
    uint64_t src_id = -1;
    uint64_t snk_id = -1;
};

/*
 * Template implementations
 */

template<typename UIntSize, typename UIntChain, size_t PageSize, size_t TiltBias>
PackedPathMerge<UIntSize, UIntChain, PageSize, TiltBias>::PackedPathMerge(const BaseGraph& graph) : PackedPathMerge(graph, nullptr) {
    
}

template<typename UIntSize, typename UIntChain, size_t PageSize, size_t TiltBias>
PackedPathMerge<UIntSize, UIntChain, PageSize, TiltBias>::PackedPathMerge(const BaseGraph& graph, const SentinelTableau& tableau) : PackedPathMerge(graph, &tableau) {
    src_id = tableau.src_id;
    snk_id = tableau.snk_id;
}

template<typename UIntSize, typename UIntChain, size_t PageSize, size_t TiltBias>
PackedPathMerge<UIntSize, UIntChain, PageSize, TiltBias>::PackedPathMerge(const BaseGraph& graph, const SentinelTableau* tableau) : path_head(graph.node_size()), g(&graph)
{
    static const bool debug = false;
    
    if (debug) {
        std::cerr << "constructing with template params " << typeid(UIntSize).name() << ", " << typeid(UIntChain).name() << ", " << PageSize << ", " << TiltBias << '\n';
        std::cerr << "recording indexes on paths\n";
    }
    
    // identify steps of paths and seed the DP
    for (uint64_t path_id = 0; path_id < graph.path_size(); ++path_id) {
        
        table.emplace_back(graph.node_size());
        index_on_path.emplace_back(graph.node_size());
        next_path_containing.emplace_back(graph.node_size());
        
        // note: start at 1 to reserve 0 as the sentinel and
        size_t index = 1;
        for (auto node_id : graph.path(path_id)) {
            for (auto next_id : graph.next(node_id)) {
                // note: we don't need the max here because we're going in increasing order,
                // so the max is the last to write
                table[path_id].set(next_id, index);
            }
            index_on_path[path_id].set(node_id, index);
            next_path_containing[path_id].set(node_id, path_head.at(node_id));
            path_head.set(node_id, path_id + 1); // start at 1 for 0 to be sentinel
            ++index;
        }
    }
    
    if (debug) {
        std::cerr << "filling the last to reach table\n";
    }
    
    // use DP to fill it in between cross-path edges
    for (auto node_id : topological_order(graph)) {
        for (auto prev_id : graph.previous(node_id)) {
            for (uint64_t path_id = 0; path_id < graph.path_size(); ++path_id) {
                auto& row = table[path_id];
                if (row.at(node_id) == 0) {
                    row.set(node_id, row.at(prev_id));
                }
                else if (row.at(prev_id) != 0) {
                    row.set(node_id, std::max(row.at(prev_id), row.at(node_id)));
                }
            }
        }
    }
    
    // add a "pseudo-path" chain for the sentinels
    if (tableau) {
        if (debug) {
            std::cerr << "creating the tableau pseudo-path\n";
        }
        
        table.emplace_back(graph.node_size());
        index_on_path.emplace_back(graph.node_size());
        next_path_containing.emplace_back(graph.node_size());
        
        index_on_path[graph.path_size()].set(tableau->src_id, 1);
        index_on_path[graph.path_size()].set(tableau->snk_id, 2);
        path_head.set(tableau->src_id, graph.path_size() + 1);
        path_head.set(tableau->snk_id, graph.path_size() + 1);
        for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
            if (node_id != tableau->src_id) {
                table[graph.path_size()].set(node_id, 1);
            }
        }
    }
    
    if (debug) {
        std::cerr << "path lists:\n";
        for (uint64_t n = 0; n < graph.node_size(); ++n) {
            std::cerr << n << ":\t" << (int64_t) path_head.at(n) << '\n';
        }
        
        std::cerr << "index on paths:\n";
        if (!index_on_path.empty()) {
            for (uint64_t n = 0; n < graph.node_size(); ++n) {
                std::cerr << n << ":";
                for (uint64_t p = 0; p < graph.path_size(); ++p) {
                    std::cerr << '\t' << (int64_t) index_on_path[p].at(n) << ',' << (int64_t) next_path_containing[p].at(n);
                }
                std::cerr << '\n';
            }
        }
        
        std::cerr << "table:\n";
        for (uint64_t n = 0; n < graph.node_size(); ++n) {
            std::cerr << n << ":";
            for (uint64_t p = 0; p < graph.path_size(); ++p) {
                std::cerr << '\t' << (int64_t) table[p].at(n);
            }
            std::cerr << '\n';
        }
    }
}

template<typename UIntSize, typename UIntChain, size_t PageSize, size_t TiltBias>
size_t PackedPathMerge<UIntSize, UIntChain, PageSize, TiltBias>::memory_size() const {
    size_t mem = path_head.memory_size() + sizeof(src_id) + sizeof(snk_id) + sizeof(g);
    for (uint64_t path_id = 0; path_id < index_on_path.size(); ++path_id) {
        mem += index_on_path[path_id].memory_size() + next_path_containing[path_id].memory_size() + table[path_id].memory_size();
    }
    return mem;
}

template<typename UIntSize, typename UIntChain, size_t PageSize, size_t TiltBias>
inline size_t PackedPathMerge<UIntSize, UIntChain, PageSize, TiltBias>::chain_size() const {
    return table.size();
}

template<typename UIntSize, typename UIntChain, size_t PageSize, size_t TiltBias>
inline size_t PackedPathMerge<UIntSize, UIntChain, PageSize, TiltBias>::node_size() const {
    return table.empty() ? 0 : table.front().size();
}

template<typename UIntSize, typename UIntChain, size_t PageSize, size_t TiltBias>
inline std::pair<uint64_t, size_t> PackedPathMerge<UIntSize, UIntChain, PageSize, TiltBias>::chain(uint64_t node_id) const {
    uint64_t path_id = path_head.at(node_id);
    if (path_id != 0) {
        return std::pair<uint64_t, size_t>(path_id - 1, index_on_path[path_id - 1].at(node_id) - 1);
    }
    return std::pair<uint64_t, size_t>(-1, -1);
}

template<typename UIntSize, typename UIntChain, size_t PageSize, size_t TiltBias>
inline size_t PackedPathMerge<UIntSize, UIntChain, PageSize, TiltBias>::predecessor_index(uint64_t node_id, uint64_t path_id) const {
    return table[path_id].at(node_id) - 1;
}

template<typename UIntSize, typename UIntChain, size_t PageSize, size_t TiltBias>
inline bool PackedPathMerge<UIntSize, UIntChain, PageSize, TiltBias>::reachable(uint64_t from_id, uint64_t to_id) const {
    uint64_t chain_from;
    size_t idx_from;
    std::tie(chain_from, idx_from) = chain(from_id);
    if (chain_from == -1) {
        // not on a path
        // TODO: do i need to maintain this anymore?
        return false;
    }
    // is the last index that can reach the target after the source?
    size_t last_idx_to_reach = table[chain_from].at(to_id) - 1;
    return (last_idx_to_reach != -1 && idx_from <= last_idx_to_reach);
}

template<typename UIntSize, typename UIntChain, size_t PageSize, size_t TiltBias>
inline size_t PackedPathMerge<UIntSize, UIntChain, PageSize, TiltBias>::index_on(uint64_t node_id, uint64_t path_id) const {
    return index_on_path[path_id].at(node_id) - 1;
}

template<typename UIntSize, typename UIntChain, size_t PageSize, size_t TiltBias>
inline std::vector<uint64_t> PackedPathMerge<UIntSize, UIntChain, PageSize, TiltBias>::chains_on(uint64_t node_id) const {
    std::vector<uint64_t> paths;
    uint64_t p = path_head.at(node_id) - 1;
    while (p != -1) {
        paths.push_back(p);
        p =  next_path_containing[p].at(node_id) - 1;
    }
    return paths;
}

template<typename UIntSize, typename UIntChain, size_t PageSize, size_t TiltBias>
inline uint64_t PackedPathMerge<UIntSize, UIntChain, PageSize, TiltBias>::node_at(uint64_t chain_id, size_t index) const {
    if (chain_id == g->path_size()) {
        return index ? snk_id : src_id;
    }
    else {
        return g->path(chain_id)[index];
    }
}

}

#endif /* centrolign_packed_path_merge_hpp */
