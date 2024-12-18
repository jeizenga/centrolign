#ifndef centrolign_superbubble_distance_oracle_hpp
#define centrolign_superbubble_distance_oracle_hpp

#include <vector>
#include <utility>
#include <tuple>
#include <limits>
#include <unordered_map>

#include "centrolign/superbubbles.hpp"
#include "centrolign/topological_order.hpp"
#include "centrolign/source_sink_graph.hpp"

namespace centrolign {

/*
 * Query structure for efficient (near O(1) in practice) minimum distance queries based on
 * a superbubble and chain decomposition
 */
class SuperbubbleDistanceOracle {
public:
    template<class Graph>
    SuperbubbleDistanceOracle(const Graph& graph);
    SuperbubbleDistanceOracle() noexcept = default;
    SuperbubbleDistanceOracle(const SuperbubbleDistanceOracle& other) noexcept = default;
    SuperbubbleDistanceOracle(SuperbubbleDistanceOracle&& other) noexcept = default;
    ~SuperbubbleDistanceOracle() = default;
    SuperbubbleDistanceOracle& operator=(const SuperbubbleDistanceOracle& other) noexcept = default;
    SuperbubbleDistanceOracle& operator=(SuperbubbleDistanceOracle&& other) noexcept = default;
    
    // return the minimum distance between from node_id1 to node_id2, or -1 if
    // node_id2 can't be reached
    size_t min_distance(uint64_t node_id1, uint64_t node_id2) const;
    
private:
    
    
    // assumes all child chains have been entered
    template<class Graph>
    void enter_net_graph_min_distances(const Graph& graph, uint64_t bub_id, const NetGraph& net_graph);
    
    // pairs of (feature id, is chain)
    std::vector<std::pair<uint64_t, bool>> path_to_root(uint64_t node_id) const;
    
    SuperbubbleTree superbubble_tree;
    
    // the narrowest superbubble that a node belongs to (given to the later superbubble in the case
    // of a shared boundary)
    std::vector<uint64_t> node_to_superbubble;
    // the index of a superbubble within its chain
    std::vector<size_t> superbubble_link_index;
    
    // maps from (feature id, is chain) pairs to distance
    std::vector<std::unordered_map<std::pair<std::pair<uint64_t, bool>, std::pair<uint64_t, bool>>, size_t>> net_graph_tables;
    
    std::vector<std::vector<size_t>> chain_prefix_sums;
};

template<class Graph>
SuperbubbleDistanceOracle::SuperbubbleDistanceOracle(const Graph& graph) :
    node_to_superbubble(graph.node_size())
{
    {
        // somewhat hack-y way to communicate the single source single sink
        // info without modifying the graph
        SourceSinkGraph<Graph> overlay(graph);
        SentinelTableau tableau;
        tableau.src_id = overlay.source_id();
        tableau.snk_id = overlay.sink_id();
        
        // find the superbubbles
        superbubble_tree = std::move(SuperbubbleTree(overlay, tableau));
    }
        
    static const bool debug = false;
    
    chain_prefix_sums.resize(superbubble_tree.chain_size());
    // 1 extra for the outer "superbubble"
    net_graph_tables.resize(superbubble_tree.structure_size() + 1);
    superbubble_link_index.resize(superbubble_tree.structure_size());
    
    for (const auto& feature : superbubble_tree.postorder()) {
        if (feature.second) {
            // it is a chain
            auto chain_id = feature.first;
            
            // fill out prefix sum
            const auto& chain = superbubble_tree.structures_inside(chain_id);
            auto& prefix_sum = chain_prefix_sums[chain_id];
            prefix_sum.resize(chain.size() + 1);
            for (size_t i = 0; i < chain.size(); ++i) {
                auto bub_id = chain[i];
                superbubble_link_index[bub_id] = i;
                uint64_t bub_src_id, bub_snk_id;
                std::tie(bub_src_id, bub_snk_id) = superbubble_tree.structure_boundaries(bub_id);
                prefix_sum[i + 1] = prefix_sum[i] + net_graph_tables[bub_id][std::make_pair(std::make_pair(bub_src_id, false),
                                                                                            std::make_pair(bub_snk_id, false))];
            }
            
            if (debug) {
                std::cerr << "prefix sum for chain " << chain_id << " from nodes " << superbubble_tree.structure_boundaries(superbubble_tree.structures_inside(chain_id).front()).first << " to " << superbubble_tree.structure_boundaries(superbubble_tree.structures_inside(chain_id).back()).second << ":\n";
                for (auto p : prefix_sum) {
                    std::cerr << ' ' << p;
                }
                std::cerr << '\n';
            }
        }
        else {
            // it is a superbubble
            auto bub_id = feature.first;
            
            // fill out distance table
            NetGraph net_graph(graph, superbubble_tree, bub_id);
            enter_net_graph_min_distances(graph, bub_id, net_graph);
        }
    }
    
    // enter the outside-bubble graph material
    NetGraph net_graph(graph, superbubble_tree);
    enter_net_graph_min_distances(graph, superbubble_tree.structure_size(), net_graph);
}


template<class Graph>
void SuperbubbleDistanceOracle::enter_net_graph_min_distances(const Graph& graph, uint64_t bub_id,
                                                              const NetGraph& net_graph) {
    
    static const bool debug = false;
    
    auto& table = net_graph_tables[bub_id];
    
    auto order = topological_order(net_graph);
    for (size_t i = 0; i < order.size(); ++i) {
        // use dynamic programming to find min distance from this node
        auto source_label = net_graph.label(order[i]);
        
        if (!source_label.second) {
            // this node is a graph node
            auto node_id = source_label.first;
            if (bub_id == superbubble_tree.structure_size() ||                    // in the outer "bubble"
                node_id == superbubble_tree.structure_boundaries(bub_id).first || // the start of this bubble
                superbubble_tree.structure_beginning_at(node_id) == -1) {         // not the start of any bubble
                // record the membership of this node to this bubble
                node_to_superbubble[node_id] = bub_id;
            }
        }
        
        // distances are from beginning of node to beginning of node
        std::vector<size_t> dp(net_graph.node_size(), -1);
        dp[order[i]] = 0;
        for (size_t j = i; j < order.size(); ++j) {
            auto net_id = order[j];
            if (dp[net_id] == -1) {
                continue;
            }
            uint64_t feature_id;
            bool is_chain;
            std::tie(feature_id, is_chain) = net_graph.label(net_id);
            size_t len;
            if (is_chain) {
                auto final_node = superbubble_tree.structure_boundaries(superbubble_tree.structures_inside(feature_id).back()).second;
                len = chain_prefix_sums[feature_id].back() + graph.label_size(final_node);
            }
            else {
                len = graph.label_size(feature_id);
            }
            size_t dist_thru = dp[net_id] + len;
            for (auto next_id : net_graph.next(net_id)) {
                dp[next_id] = std::min(dp[next_id], dist_thru);
            }
        }
        
        // enter the finite distances in the table
        for (size_t j = i; j < order.size(); ++j) {
            auto net_id = order[j];
            if (dp[net_id] != -1) {
                table[std::make_pair(source_label, net_graph.label(net_id))] = dp[net_id];
            }
        }
    }
    
    if (debug) {
        std::vector<std::pair<std::pair<uint64_t, bool>, std::pair<uint64_t, bool>>> keys;
        for (const auto& rec : table) {
            keys.push_back(rec.first);
        }
        std::sort(keys.begin(), keys.end());
        std::cerr << "table for superbubble " << bub_id;
        if (bub_id < superbubble_tree.structure_size()) {
            std::cerr << " between nodes " << superbubble_tree.structure_boundaries(bub_id).first << " and " << superbubble_tree.structure_boundaries(bub_id).second << '\n';
        }
        else {
            std::cerr << ", the exterior bubble\n";
        }
        for (auto key : keys) {
            std::cerr << key.first.first << "," << key.first.second << " -> " << key.second.first << "," << key.second.second << ": " << table.at(key) << '\n';
        }
    }
}


}

#endif /* centrolign_superbubble_distance_oracle_hpp */
