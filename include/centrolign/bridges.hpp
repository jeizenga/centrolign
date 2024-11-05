#ifndef centrolign_bridges_hpp
#define centrolign_bridges_hpp

#include <vector>
#include <cstdint>
#include <unordered_set>

#include "centrolign/labeled_graph.hpp"


namespace centrolign {

// identify edges that, if removed, increase the number of connected components
// in the graph with Schmidt's (2013) algorithm. edges are returned oriented forward
template<class Graph>
std::vector<std::pair<uint64_t, uint64_t>> bridges(const Graph& graph);

// get graphs that are isomorphic to each 2-connected component, where nodes are
// labeled with the corresponding node ID from the original graph
template<class Graph>
std::vector<LabeledGraph<uint64_t>> bridge_components(const Graph& graph);








/**
 * Template implementations
 */

template<class Graph>
std::vector<std::pair<uint64_t, uint64_t>> bridges(const Graph& graph) {
    
    std::vector<std::pair<uint64_t, uint64_t>> bridges_found;
    
    if (graph.node_size() != 0) {
        
        std::vector<uint64_t> dfs_order(graph.node_size(), -1);
        // records of (parent, traversed upward, DFS index, edges downward) with edges as (node id, to left)
        std::vector<std::tuple<uint64_t, bool, size_t, std::vector<std::pair<uint64_t, bool>>>>
        dfs_tree(graph.node_size(), std::tuple<uint64_t, bool, size_t, std::vector<std::pair<uint64_t, bool>>>(-1, false, -1, std::vector<std::pair<uint64_t, bool>>()));
        
        size_t dfs_idx = 0;
        for (uint64_t root_id = 0; root_id < graph.node_size(); ++root_id) {
            
            if (std::get<2>(dfs_tree[root_id]) != -1) {
                // we found this node doing DFS from a different root
                continue;
            }
            
            // records of (node, next edge idx, edges), with edges as (node_id, to left)
            std::vector<std::tuple<uint64_t, size_t, std::vector<std::pair<uint64_t, bool>>>> stack;
            
            // add to stack and DFS tree
            auto enqueue = [&](uint64_t node_id, uint64_t parent) {
                
                stack.emplace_back(node_id, 0, std::vector<std::pair<uint64_t, bool>>());
                // note: this will duplicate self-looping edges, but they don't affect this algorithm
                // so it's okay (self-loops are irrelevant to bridges and to DFS)
                for (auto prev_id : graph.previous(node_id)) {
                    std::get<2>(stack.back()).emplace_back(prev_id, true);
                }
                for (auto next_id : graph.next(node_id)) {
                    std::get<2>(stack.back()).emplace_back(next_id, false);
                }
                std::get<0>(dfs_tree[node_id]) = parent;
                std::get<2>(dfs_tree[node_id]) = dfs_idx;
                dfs_order[dfs_idx++] = node_id;
            };
            
            // initialize the stack at the root
            enqueue(root_id, -1);
            // do DFS
            while (!stack.empty()) {
                
                auto& top = stack.back();
                
                if (std::get<1>(top) == std::get<2>(top).size()) {
                    // exhausted edges from this node
                    stack.pop_back();
                    continue;
                }
                
                auto next_edge = std::get<2>(top)[std::get<1>(top)++];
                
                if (std::get<2>(dfs_tree[next_edge.first]) == -1) {
                    // the next node hasn't been visited yet
                    std::get<3>(dfs_tree[std::get<0>(top)]).emplace_back(next_edge);
                    enqueue(next_edge.first, std::get<0>(top));
                }
            }
        }
        
        // walk the chain loops
        for (size_t i = 0; i < dfs_order.size(); ++i) {
            
            uint64_t node_id = dfs_order[i];
            const auto& tree_record = dfs_tree[node_id];
            std::vector<uint64_t> chain_heads;
            
            size_t j = 0;
            // note: this iteration order is the same as the nodes are added to the stack
            // frame, which is essential to this procedure functioning correctly
            for (bool left : {true, false}) {
                for (uint64_t next_id : (left ? graph.previous(node_id) : graph.next(node_id))) {
                    if (j < std::get<3>(tree_record).size() && std::get<3>(tree_record)[j] == std::make_pair(next_id, left)) {
                        // this edge was used in the DFS tree
                        ++j;
                    }
                    else if (std::get<2>(dfs_tree[next_id]) > i) {
                        // this edge points down the tree rather than up it
                        chain_heads.emplace_back(next_id);
                    }
                }
            }
            
            // mark nodes in the untraversed portion of this chain's cycle as traversed
            for (uint64_t chain_cursor : chain_heads) {
                while (chain_cursor != node_id && !std::get<1>(dfs_tree[chain_cursor])) {
                    std::get<1>(dfs_tree[chain_cursor]) = true;
                    chain_cursor = std::get<0>(dfs_tree[chain_cursor]);
                }
            }
        }
        
        // find the edges that were never part of a cycle
        for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
            
            const auto& tree_record = dfs_tree[node_id];
            for (const auto& edge : std::get<3>(tree_record)) {
                if (!std::get<1>(dfs_tree[edge.first])) {
                    // this edge is a bridge
                    if (edge.second) {
                        bridges_found.emplace_back(edge.first, node_id);
                    }
                    else {
                        bridges_found.emplace_back(node_id, edge.first);
                    }
                }
            }
        }
    }
    
    return bridges_found;
}


template<class Graph>
std::vector<LabeledGraph<uint64_t>> bridge_components(const Graph& graph) {
    
    std::vector<std::unordered_set<uint64_t>> blocked_edges(graph.node_size());
    for (const auto& edge : bridges(graph)) {
        blocked_edges[edge.first].insert(edge.second);
        blocked_edges[edge.second].insert(edge.first);
    }
    
    std::vector<LabeledGraph<uint64_t>> components;
    
    std::vector<bool> traversed(graph.node_size(), false);
    
    std::vector<uint64_t> forward_translation(graph.node_size(), -1);
    
    for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
        
        if (traversed[node_id]) {
            continue;
        }
        
        components.emplace_back();
        auto& component = components.back();
        
        std::vector<uint64_t> stack(1, node_id);
        traversed[node_id] = true;
        uint64_t seed_id = component.add_node(node_id);
        forward_translation[node_id] = seed_id;
        while (!stack.empty()) {
            
            auto here_id = stack.back();
            stack.pop_back();
            
            for (bool left : {true, false}) {
                for (auto adj_id : (left ? graph.previous(here_id) : graph.next(here_id))) {
                    if (blocked_edges[here_id].count(adj_id)) {
                        // this edge crosses a bridge
                        continue;
                    }
                    
                    if (!traversed[adj_id]) {
                        uint64_t new_id = component.add_node(adj_id);
                        forward_translation[adj_id] = new_id;
                        traversed[adj_id] = true;
                    }
                    
                    if (adj_id < here_id || (adj_id == here_id && left)) {
                        if (left) {
                            component.add_edge(forward_translation[adj_id], forward_translation[here_id]);
                        }
                        else {
                            component.add_edge(forward_translation[here_id], forward_translation[adj_id]);
                        }
                    }
                }
            }
        }
    }
    
    return components;
}

}

#endif /* centrolign_bridges_hpp */
