#ifndef centrolign_superbubble_distances_hpp
#define centrolign_superbubble_distances_hpp

#include <vector>
#include <utility>
#include <tuple>
#include <limits>

#include "centrolign/superbubbles.hpp"
#include "centrolign/topological_order.hpp"

namespace centrolign {

class SuperbubbleDistances {
public:
    template<class Graph>
    SuperbubbleDistances(const SuperbubbleTree& superbubbles,
                         const Graph& graph);
    
    SuperbubbleDistances() = default;
    ~SuperbubbleDistances() = default;
    
    // return the minimum and maximum length walk through a feature, including both
    // the start and end nodes
    inline std::pair<size_t, size_t> superbubble_min_max_dist(uint64_t superbubble_id) const;
    inline std::pair<size_t, size_t> chain_min_max_dist(uint64_t chain_id) const;
    
private:
    
    // the min and max walk distance walks for superbubbles
    std::vector<std::pair<size_t, size_t>> superbubble_dists;
    
    // the min and max walk distance walks for chains
    std::vector<std::pair<size_t, size_t>> chain_dists;
    
};

template<class Graph>
SuperbubbleDistances::SuperbubbleDistances(const SuperbubbleTree& superbubbles,
                                           const Graph& graph) :
    superbubble_dists(superbubbles.superbubble_size(), std::pair<size_t, size_t>(0, 0)),
    chain_dists(superbubbles.chain_size(), std::pair<size_t, size_t>(0, 0))
{
    static const bool debug = false;
    
    for (const auto& feature : superbubbles.postorder()) {
        if (feature.second) {
            if (debug) {
                std::cerr << "reached postorder of chain " << feature.first << " containing snarls:";
                for (auto s : superbubbles.superbubbles_inside(feature.first)) {
                    std::cerr << ' ' << s;
                }
                std::cerr << '\n';
            }
            
            // is a chain, sum over superbubbles
            auto& dists_here = chain_dists[feature.first];
            const auto& links = superbubbles.superbubbles_inside(feature.first);
            for (size_t i = 0; i < links.size(); ++i) {
                auto& bub_dists = superbubble_dists[links[i]];
                dists_here.first += bub_dists.first;
                dists_here.second += bub_dists.second;
                if (i != 0) {
                    // correct for the overlap in the two bubbles
                    size_t overlap = graph.label_size(superbubbles.superbubble_boundaries(links[i]).first);
                    dists_here.first -= overlap;
                    dists_here.second -= overlap;
                }
            }
            if (debug) {
                std::cerr << "got min/max distance " << dists_here.first << ", " << dists_here.second << '\n';
            }
        }
        else {
            if (debug) {
                uint64_t s, e;
                std::tie(s, e) = superbubbles.superbubble_boundaries(feature.first);
                std::cerr << "reached postorder of snarl " << feature.first << " with boundaries " << s << " and " << e << '\n';
            }
            
            // is a superbubble, do DP over the net graph
            
            NetGraph netgraph(graph, superbubbles, feature.first);
            
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
            auto& bub_dists = superbubble_dists[feature.first];
            bub_dists.first = dp[order.back()].first;
            bub_dists.second = dp[order.back()].second;
            
            if (debug) {
                std::cerr << "got min/max distance " << bub_dists.first << ", " << bub_dists.second << '\n';
            }
        }
    }
    
    for (uint64_t chain_id = 0; chain_id < superbubbles.chain_size(); ++chain_id) {
        if (superbubbles.superbubble_containing(chain_id) != -1) {
            // this is not a top-level chain
            continue;
        }
        
        if (debug) {
            std::cerr << "indexing tree with top level chain " << chain_id << '\n';
        }
        
        // stack of (id, is chain, children have been added)
        std::vector<std::tuple<uint64_t, bool, bool>> stack;
        stack.emplace_back(chain_id, true, false);
        
        while (!stack.empty()) {
            if (std::get<2>(stack.back())) {
                // the children have already been added, we can assume postorder is complete
                
                // we've completed measuring distances at this feature
                stack.pop_back();
            }
            else {
                // add the children to the stack
                std::get<2>(stack.back()) = true;
                if (std::get<1>(stack.back())) {
                    // is a chain, add superbubble children
                    for (auto child_id : superbubbles.superbubbles_inside(std::get<0>(stack.back()))) {
                        stack.emplace_back(child_id, false, false);
                    }
                }
                else {
                    // is a superbubble, add chain children
                    for (auto child_id : superbubbles.chains_inside(std::get<0>(stack.back()))) {
                        stack.emplace_back(child_id, true, false);
                    }
                }
            }
        }
    }
}

inline std::pair<size_t, size_t> SuperbubbleDistances::superbubble_min_max_dist(uint64_t superbubble_id) const {
    return superbubble_dists[superbubble_id];
}

inline std::pair<size_t, size_t> SuperbubbleDistances::chain_min_max_dist(uint64_t chain_id) const {
    return chain_dists[chain_id];
}

}

#endif /* centrolign_superbubble_distances_hpp */
