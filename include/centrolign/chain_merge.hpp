#ifndef centrolign_chain_merge_hpp
#define centrolign_chain_merge_hpp

#include <vector>
#include <cstdint>
#include <iostream>
#include <tuple>

#include "centrolign/topological_order.hpp"

namespace centrolign {

class Anchorer; // forward declaration

// TODO: is it possible to merge two of these when combining graphs?

/*
 * Data structure for reachability using a path cover
 */
class ChainMerge {
public:
    template<class PGraph>
    ChainMerge(const PGraph& graph);
    ChainMerge() = default;
    ~ChainMerge() = default;
    
    // number of chains
    inline size_t chain_size() const;
    
    // the chain that the node belongs to
    inline uint64_t chain(uint64_t node_id) const;
    
    // the index of the node within its chain
    inline size_t chain_index(uint64_t node_id) const;
    
    // the indexes (within their chain) of the nearest predecessors of this
    // node in each chain
    inline const std::vector<size_t>& predecessor_indexes(uint64_t node_id) const;
    
    // return true if both nodes are on a path and one can reach the other
    inline bool reachable(uint64_t from_id, uint64_t to_id) const;
    
    // generate edges to nodes that are nearest predecessors within one
    // of the chains
    std::vector<std::vector<uint64_t>> chain_forward_edges() const;
    
private:
    
    static const bool debug = false;
    
    // the DP table: for each node, for each chain, the latest index in the chain
    // that can reach that node
    std::vector<std::vector<size_t>> table;
    // the assignment of nodes to chains
    std::vector<std::pair<uint64_t, size_t>> node_to_chain;

};

/*
 * Template implementations
 */
template<class PGraph>
ChainMerge::ChainMerge(const PGraph& graph) :
    node_to_chain(graph.node_size(), std::pair<uint64_t, size_t>(-1, -1)),
    table(graph.node_size(), std::vector<size_t>(graph.path_size(), -1))
{
    
    // assign nodes to a chain
    for (uint64_t path_id = 0; path_id < graph.path_size(); ++path_id) {
        size_t index = 0;
        for (auto node_id : graph.path(path_id)) {
            if (node_to_chain[node_id].first == -1) {
                node_to_chain[node_id] = std::make_pair(path_id, index++);
            }
        }
    }
    
    // find the last index in each chain that can reach each node
    for (uint64_t node_id : topological_order(graph)) {
        auto& row = table[node_id];
        uint64_t chain_id;
        size_t idx;
        std::tie(chain_id, idx) = node_to_chain[node_id];
        if (chain_id == -1) {
            // not covered by path, probably a sentinel node
            continue;
        }
        
        // propagate to later nodes with DP
        for (uint64_t next_id : graph.next(node_id)) {
            auto& next_row = table[next_id];
            uint64_t next_chain_id;
            int64_t next_idx;
            std::tie(next_chain_id, next_idx) = node_to_chain[next_id];
            if (next_chain_id == -1) {
                // not covered by path, probably a sentinel node
                continue;
            }
            for (uint64_t c = 0; c < row.size(); ++c) {
                // use signed comparison so that -1s get overwritten
                if (c == chain_id) {
                    next_row[c] = std::max<int64_t>(next_row[c], idx);
                }
                else {
                    next_row[c] = std::max<int64_t>(next_row[c], row[c]);
                }
            }
        }
    }
    
    if (debug) {
        std::vector<std::vector<uint64_t>> chains(graph.path_size());
        for (uint64_t i = 0; i < node_to_chain.size(); ++i) {
            uint64_t c;
            int64_t idx;
            std::tie(c, idx) = node_to_chain[i];
            if (c != -1) {
                while (chains[c].size() <= idx) {
                    chains[c].emplace_back(-1);
                }
                chains[c][idx] = i;
            }
        }
        std::cerr << "chains:\n";
        for (size_t i = 0; i < chains.size(); ++i) {
            std::cerr << i << ':';
            for (size_t j = 0; j < chains[i].size(); ++j) {
                std::cerr << ' ' << chains[i][j] << '(' << j << ')';
            }
            std::cerr << '\n';
        }
        std::cerr << "table:\n";
        for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
            std::cerr << node_id << ':';
            for (auto i : table[node_id]) {
                std::cerr << '\t' << (int64_t) i;
            }
            std::cerr << '\n';
        }
        std::cerr << '\n';
    }
}

inline size_t ChainMerge::chain_size() const {
    return table.empty() ? 0 : table.front().size();
}

inline uint64_t ChainMerge::chain(uint64_t node_id) const {
    return node_to_chain[node_id].first;
}

inline size_t ChainMerge::chain_index(uint64_t node_id) const {
    return node_to_chain[node_id].second;
}

inline const std::vector<size_t>& ChainMerge::predecessor_indexes(uint64_t node_id) const {
    return table[node_id];
}

inline bool ChainMerge::reachable(uint64_t from_id, uint64_t to_id) const {
    uint64_t chain_from;
    size_t idx_from;
    std::tie(chain_from, idx_from) = node_to_chain[from_id];
    if (chain_from == -1 && node_to_chain[to_id].first == -1) {
        // one or the other is not on a path
        return false;
    }
    // is the last index that can reach the target after the source?
    size_t last_idx_to_reach = table[to_id][chain_from];
    return (last_idx_to_reach != -1 && idx_from <= last_idx_to_reach);
}

}

#endif /* centrolign_chain_merge_hpp */
