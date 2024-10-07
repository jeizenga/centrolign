#ifndef centrolign_structure_distances_hpp
#define centrolign_structure_distances_hpp

#include <vector>
#include <utility>
#include <tuple>
#include <limits>
#include <queue>

#include "centrolign/structure_tree.hpp"
#include "centrolign/topological_order.hpp"
#include "centrolign/is_acyclic.hpp"
#include "centrolign/snarls.hpp"

namespace centrolign {

/*
 * Index of the maximum and minimum distance across superbubble features
 */
template<class StructureType, bool AssumeAcyclic>
class StructureDistances;

using SuperbubbleDistances = StructureDistances<SuperbubbleTree, true>;
using SnarlDistances = StructureDistances<SnarlTree, false>;

template<class StructureType, bool AssumeAcyclic>
class StructureDistances {
public:
    template<class Graph>
    StructureDistances(const StructureType& structures,
                       const Graph& graph);
    
    StructureDistances() = default;
    ~StructureDistances() = default;
    
    // return the minimum and maximum length walk through a feature, including both
    // the start and end nodes
    inline std::pair<size_t, size_t> structure_min_max_dist(uint64_t struct_id) const;
    inline std::pair<size_t, size_t> chain_min_max_dist(uint64_t chain_id) const;
    
private:
    
    // the min and max walk distance walks for superbubbles
    std::vector<std::pair<size_t, size_t>> structure_dists;
    
    // the min and max walk distance walks for chains
    std::vector<std::pair<size_t, size_t>> chain_dists;
    
};

template<class StructureType, bool AssumeAcyclic>
template<class Graph>
StructureDistances<StructureType, AssumeAcyclic>::StructureDistances(const StructureType& structures,
                                                                     const Graph& graph) :
    structure_dists(structures.structure_size(), std::pair<size_t, size_t>(0, 0)),
    chain_dists(structures.chain_size(), std::pair<size_t, size_t>(0, 0))
{
    static const bool debug = false;
    
    for (const auto& feature : structures.postorder()) {
        if (feature.second) {
            // chain
            if (debug) {
                std::cerr << "reached postorder of chain " << feature.first << " containing snarls:";
                for (auto s : structures.structures_inside(feature.first)) {
                    std::cerr << ' ' << s;
                }
                std::cerr << '\n';
            }
            
            // is a chain, sum over superbubbles
            auto& dists_here = chain_dists[feature.first];
            const auto& links = structures.structures_inside(feature.first);
            for (size_t i = 0; i < links.size(); ++i) {
                auto& struct_dists = structure_dists[links[i]];
                dists_here.first += struct_dists.first;
                if (!AssumeAcyclic && dists_here.second == -1 || struct_dists.second == -1) {
                    dists_here.second = -1;
                }
                else {
                    dists_here.second += struct_dists.second;
                }
                if (i != 0) {
                    // correct for the overlap in the two bubbles
                    size_t overlap = graph.label_size(structures.structure_boundaries(links[i]).first);
                    dists_here.first -= overlap;
                    if (AssumeAcyclic || dists_here.second != -1) {
                        dists_here.second -= overlap;
                    }
                }
            }
            if (debug) {
                std::cerr << "got min/max distance " << dists_here.first << ", " << dists_here.second << '\n';
            }
        }
        else {
            // snarl/superbubble
            if (debug) {
                uint64_t s, e;
                std::tie(s, e) = structures.structure_boundaries(feature.first);
                std::cerr << "reached postorder of snarl " << feature.first << " with boundaries " << s << " and " << e << '\n';
            }
            
            // is a superbubble, do DP over the net graph
            
            NetGraph netgraph(graph, structures, feature.first);
            
            // figure out if we're in a cyclic cnarl
            bool acyclic = true;
            if (!AssumeAcyclic) {
                for (auto chain_id : structures.chains_inside(feature.first)) {
                    if (chain_dists[chain_id].second == -1) {
                        acyclic = false;
                        break;
                    }
                }
                if (acyclic) {
                    acyclic = is_acyclic(netgraph);
                }
            }
            
            auto& struct_dists = structure_dists[feature.first];
            
            if (acyclic) {
                // acyclic graph, do linear DP algorithm
                
                std::vector<std::pair<int64_t, int64_t>> dp(netgraph.node_size(),
                                                            std::pair<int64_t, int64_t>(std::numeric_limits<int64_t>::max(), -1));
                
                auto order = topological_order(netgraph);
                
                // initial conditions
                dp[order.front()].first = graph.label_size(order.front());
                dp[order.front()].second = graph.label_size(order.front());
                
                // dynamic programming
                for (size_t i = 0; i < order.size(); ++i) {
                    auto& dists_here = dp[order[i]];
                    for (auto next_id : netgraph.next(order[i])) {
                        uint64_t feature_id;
                        bool is_chain;
                        std::tie(feature_id, is_chain) = netgraph.label(next_id);
                        int64_t min_dist_thru, max_dist_thru;
                        if (is_chain) {
                            auto chain_mm_dist = chain_min_max_dist(feature_id);
                            min_dist_thru = dists_here.first + chain_mm_dist.first;
                            max_dist_thru = dists_here.second + chain_mm_dist.second;
                        }
                        else {
                            min_dist_thru = dists_here.first + graph.label_size(feature_id);
                            max_dist_thru = dists_here.second + graph.label_size(feature_id);
                        }
                        
                        auto& dists_next = dp[next_id];
                        if (min_dist_thru < dists_next.first) {
                            dists_next.first = min_dist_thru;
                        }
                        if (max_dist_thru > dists_next.second) {
                            dists_next.second = max_dist_thru;
                        }
                    }
                }
                
                // record the final conditions in the index
                struct_dists.first = dp[order.back()].first;
                struct_dists.second = dp[order.back()].second;
                
                if (debug) {
                    std::cerr << "got min/max distance " << struct_dists.first << ", " << struct_dists.second << '\n';
                }
            }
            else {
                // cyclic graph, do dijkstra
                
                std::vector<bool> popped(netgraph.node_size(), false);
                
                std::priority_queue<std::pair<size_t, uint64_t>, std::vector<std::pair<size_t, uint64_t>>, std::greater<std::pair<size_t, uint64_t>>> queue;
                for (uint64_t node_id = 0; node_id < netgraph.node_size(); ++node_id) {
                    if (netgraph.previous_size(node_id) == 0) {
                        queue.emplace(graph.label_size(netgraph.label(node_id).first), node_id);
                        break;
                    }
                }
                
                std::vector<size_t> distance(netgraph.node_size(), 0);
                
                while (!queue.empty()) {
                    auto top = queue.top();
                    queue.pop();
                    if (popped[top.second]) {
                        continue;
                    }
                    popped[top.second] = true;
                    
                    distance[top.second] = top.first;
                    
                    for (auto next_id : netgraph.next(top.second)) {
                        size_t size;
                        uint64_t feature_id;
                        bool is_chain;
                        std::tie(feature_id, is_chain) = netgraph.label(next_id);
                        if (is_chain) {
                            size = chain_dists[feature_id].first;
                        }
                        else {
                            size = graph.label_size(next_id);
                        }
                        queue.emplace(top.first + size, next_id);
                    }
                }
                
                // find the sink and record its distances
                for (uint64_t node_id = 0; node_id < netgraph.node_size(); ++node_id) {
                    if (netgraph.next_size(node_id) == 0) {
                        struct_dists.first = distance[node_id];
                        break;
                    }
                }
                // the maximum distance is undefined
                struct_dists.second = -1;
            }
        }
    }
}

template<class StructureType, bool AssumeAcyclic>
inline std::pair<size_t, size_t> StructureDistances<StructureType, AssumeAcyclic>::structure_min_max_dist(uint64_t struct_id) const {
    return structure_dists[struct_id];
}

template<class StructureType, bool AssumeAcyclic>
inline std::pair<size_t, size_t> StructureDistances<StructureType, AssumeAcyclic>::chain_min_max_dist(uint64_t chain_id) const {
    return chain_dists[chain_id];
}

}

#endif /* centrolign_structure_distances_hpp */
