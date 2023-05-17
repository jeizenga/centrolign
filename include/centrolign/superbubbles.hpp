#ifndef centrolign_superbubbles_hpp
#define centrolign_superbubbles_hpp

#include <vector>
#include <utility>
#include <cstdint>
#include <limits>
#include <iostream>
#include <stdexcept>

#include "centrolign/topological_order.hpp"
#include "centrolign/modify_graph.hpp"

namespace centrolign {

/*
 * Data structure that finds superbubbles and provides navigation
 * methods for the superbubble/chain tree
 */
class SuperbubbleTree {
public:
    
    // construct from a single-source, single-sink graph
    template<class Graph>
    SuperbubbleTree(const Graph& graph);
    
    // construct from a graph with sentinel sink and source nodes, which
    // will be removed from the superbubble decomposition
    template<class Graph>
    SuperbubbleTree(const Graph& graph, const SentinelTableau& tableau);
    
    SuperbubbleTree() = default;
    ~SuperbubbleTree() = default;
    
    inline size_t chain_size() const;
    inline size_t superbubble_size() const;
    
    // access superbubbles by graph node ID (returns -1 is there is none)
    inline uint64_t superbubble_beginning_at(uint64_t node_id) const;
    inline uint64_t superbubble_ending_at(uint64_t node_id) const;
    
    // returns pair of graph node IDs
    inline const std::pair<uint64_t, uint64_t>& superbubble_boundaries(uint64_t superbubble_id) const;
    
    // returns vector of chain IDs, in no particular order
    inline const std::vector<uint64_t> chains_inside(uint64_t superbubble_id) const;
    
    // returns chain ID
    inline uint64_t chain_containing(uint64_t superbubble_id) const;
    
    // returns vector of superbubble IDs, in order from start to end
    inline const std::vector<uint64_t>& superbubbles_inside(uint64_t chain_id) const;
    
    // returns superbubble ID, or -1 is there is none
    inline uint64_t superbubble_containing(uint64_t chain_id) const;
    
protected:
    
    template<class Graph>
    SuperbubbleTree(const Graph& graph, const SentinelTableau* tableau);
    
    // find superbubbles with algorithm of Gartner, et al. (2018) for DAGs
    template<class Graph>
    static std::vector<std::pair<uint64_t, uint64_t>> find_superbubbles(const Graph& graph);
    
    struct Chain {
        Chain() = default;
        ~Chain() = default;
        
        std::vector<uint64_t> superbubble_ids;
        uint64_t parent = -1;
    };
    
    struct Superbubble {
        Superbubble(const std::pair<uint64_t, uint64_t>& superbubble) : boundaries(superbubble) {}
        Superbubble() = default;
        ~Superbubble() = default;
        
        std::pair<uint64_t, uint64_t> boundaries;
        uint64_t parent = -1;
        std::vector<uint64_t> chain_ids;
    };
    
    std::vector<Superbubble> superbubbles;
    std::vector<Chain> chains;
    std::vector<uint64_t> superbubble_endings;
    std::vector<uint64_t> superbubble_beginnings;
};





/*
 * Template implementations
 */


template<class Graph>
SuperbubbleTree::SuperbubbleTree(const Graph& graph) : SuperbubbleTree(graph, nullptr) {
    // nothing else
}

template<class Graph>
SuperbubbleTree::SuperbubbleTree(const Graph& graph,
                                 const SentinelTableau& tableau) : SuperbubbleTree(graph, &tableau) {
    // nothing else
}


template<class Graph>
SuperbubbleTree::SuperbubbleTree(const Graph& graph, const SentinelTableau* tableau) {
    
    static const bool debug = false;
    
    superbubble_beginnings.resize(graph.node_size(), -1);
    superbubble_endings.resize(graph.node_size(), -1);
    
    // identify and collect the superbubbles
    for (const auto& sup_bub : find_superbubbles(graph)) {
        if (tableau && (tableau->src_id == sup_bub.first || tableau->snk_id == sup_bub.second)) {
            // superbubbles including the sentinels are not interesting to us
            continue;
        }
        
        if (debug) {
            std::cerr << "recording bubble " << sup_bub.first << " " << sup_bub.second << '\n';
        }
        
        superbubble_beginnings[sup_bub.first] = superbubbles.size();
        superbubble_endings[sup_bub.second] = superbubbles.size();
        
        superbubbles.emplace_back(sup_bub);
    }
    
    // form them into chains
    for (uint64_t bub_id = 0; bub_id < superbubbles.size(); ++bub_id) {
        if (superbubbles[bub_id].parent != -1) {
            // we already chained this one
            continue;
        }
        
        uint64_t chain_id = chains.size();
        chains.emplace_back();
        auto& chain = chains.back();
        
        // extend to the left
        chain.superbubble_ids.push_back(bub_id);
        superbubbles[bub_id].parent = chain_id;
        uint64_t bub_id_here = superbubble_ending_at(superbubble_boundaries(bub_id).first);
        while (bub_id_here != -1) {
            chain.superbubble_ids.push_back(bub_id_here);
            superbubbles[bub_id_here].parent = chain_id;
            bub_id_here = superbubble_ending_at(superbubble_boundaries(bub_id_here).first);
        }
        // we build this backwards, now reverse it
        std::reverse(chain.superbubble_ids.begin(), chain.superbubble_ids.end());
        // extend to the right
        bub_id_here = superbubble_beginning_at(superbubble_boundaries(bub_id).second);
        while (bub_id_here != -1) {
            chain.superbubble_ids.push_back(bub_id_here);
            superbubbles[bub_id_here].parent = chain_id;
            bub_id_here = superbubble_beginning_at(superbubble_boundaries(bub_id_here).second);
        }
    }
    
    // identify the chains within superbubbles
    
    // we'll use a global bank for the traversed status so that we don't need to use
    // hash sets in the inner loop
    std::vector<bool> traversed(graph.node_size(), false);
    for (uint64_t bub_id = 0; bub_id < superbubbles.size(); ++bub_id) {
        
        auto& superbubble = superbubbles[bub_id];
                
        // DFS to find contained chains
        std::vector<uint64_t> stack(1, superbubble.boundaries.first);
        while (!stack.empty()) {
            auto node_id = stack.back();
            stack.pop_back();
            for (auto next_id : graph.next(node_id)) {
                if (next_id == superbubble.boundaries.second || traversed[next_id]) {
                    // don't go into the end of the superbubble
                    continue;
                }
                traversed[next_id] = true;
                auto next_bub_id = superbubble_beginning_at(next_id);
                if (next_bub_id != -1) {
                    // record this as a contained chain
                    auto chain_id = chain_containing(next_bub_id);
                    chains[chain_id].parent = bub_id;
                    superbubble.chain_ids.push_back(chain_id);

                    // jump to the end of the chain
                    auto final_bub_id = superbubbles_inside(chain_id).back();
                    auto final_node_id = superbubble_boundaries(final_bub_id).second;
                    traversed[final_node_id] = true;
                    stack.push_back(final_node_id);
                }
                else {
                    // just a regular node
                    stack.push_back(next_id);
                }
            }
        }
        
        
    }
}

template<class Graph>
std::vector<std::pair<uint64_t, uint64_t>> SuperbubbleTree::find_superbubbles(const Graph& graph) {
    
    static const bool debug = false;
    
    std::vector<std::pair<uint64_t, uint64_t>> superbubbles;
    
    // TODO: i'm pretty sure the depth-first Kahn's algorithm also works here,
    // so i'm not using the DFS postorder like in the paper
    auto order = topological_order(graph);
    std::vector<size_t> index(order.size());
    size_t num_sources = 0;
    size_t num_sinks = 0;
    for (size_t i = 0; i < order.size(); ++i) {
        index[order[i]] = i;
        if (graph.previous_size(i) == 0) {
            ++num_sources;
        }
        if (graph.next_size(i) == 0) {
            ++num_sinks;
        }
    }
    
    if (num_sources != 1 || num_sinks != 1) {
        throw std::runtime_error("Can only find superbubbles in single-source, single-sink graphs");
    }
    
    if (debug) {
        std::cerr << "got postorder:\n";
        for (auto n : order) {
            std::cerr << ' ' << n;
        }
        std::cerr << '\n';
    }
    
    // indexes of the candidate ends
    std::vector<int64_t> candidate_stack;
    // the furthest backward reach (filled in as we go at candidate ends)
    std::vector<int64_t> backward_reach(order.size(),
                                        std::numeric_limits<int64_t>::max());
    for (int64_t i = order.size() - 1; i >= 0; --i) {
        
        if (debug) {
            std::cerr << "at index " << i << ", node " << order[i] << ", stack:\n";
            for (auto t : candidate_stack) {
                std::cerr << ' ' << t;
            }
            std::cerr << '\n';
        }
        
        int64_t forward_reach = -1;
        for (auto node_id : graph.next(order[i])) {
            forward_reach = std::max<int64_t>(forward_reach, index[node_id]);
        }
        if (forward_reach == i + 1) {
            candidate_stack.push_back(i + 1);
            if (debug) {
                std::cerr << "identify " << (i + 1) << " (node " << order[i + 1] <<  ") as a candidate end\n";
            }
        }
        
        while (!candidate_stack.empty() && forward_reach > candidate_stack.back()) {
            // this node reaches outside the candidate interval, so the nearest candidate
            // is not a superbubble
            auto invalid_candidate = candidate_stack.back();
            candidate_stack.pop_back();
            if (!candidate_stack.empty()) {
                // give backward reach info to next nearest candidate
                backward_reach[candidate_stack.back()] = std::min(backward_reach[candidate_stack.back()],
                                                                  backward_reach[invalid_candidate]);
            }
            if (debug) {
                std::cerr << "candidate " << invalid_candidate << " (node " << order[invalid_candidate] <<  ") fails forward reach test against " << forward_reach << '\n';
            }
        }
        if (debug){
            if (!candidate_stack.empty()) {
                std::cerr << "comparing position " << i << " to backward reach of " << candidate_stack.back() << ": " << backward_reach[candidate_stack.back()] << '\n';
            }
        }
        
        if (!candidate_stack.empty() && backward_reach[candidate_stack.back()] == i) {
            // we found a superbubble interval
            auto confirmed_candidate = candidate_stack.back();
            superbubbles.emplace_back(order[i], order[confirmed_candidate]);
            candidate_stack.pop_back();
            if (!candidate_stack.empty()) {
                // give backward reach info to candidate parent
                backward_reach[candidate_stack.back()] = std::min(backward_reach[candidate_stack.back()],
                                                                  backward_reach[confirmed_candidate]);
            }
            if (debug) {
                std::cerr << "candidate " << confirmed_candidate << " (node " << order[confirmed_candidate] <<  ") is confirmed as an exit\n";
            }
        }
        for (auto node_id : graph.previous(order[i])) {
            backward_reach[i] = std::min<int64_t>(backward_reach[i], index[node_id]);
        }
        if (!candidate_stack.empty()) {
            backward_reach[candidate_stack.back()] = std::min(backward_reach[candidate_stack.back()],
                                                              backward_reach[i]);
        }
    }
    
    return superbubbles;
}

inline size_t SuperbubbleTree::chain_size() const {
    return chains.size();
}

inline size_t SuperbubbleTree::superbubble_size() const {
    return superbubbles.size();
}

inline uint64_t SuperbubbleTree::superbubble_ending_at(uint64_t node_id) const {
    return superbubble_endings[node_id];
}

inline uint64_t SuperbubbleTree::superbubble_beginning_at(uint64_t node_id) const {
    return superbubble_beginnings[node_id];
}

inline const std::pair<uint64_t, uint64_t>& SuperbubbleTree::superbubble_boundaries(uint64_t superbubble_id) const {
    return superbubbles[superbubble_id].boundaries;
}

inline const std::vector<uint64_t> SuperbubbleTree::chains_inside(uint64_t superbubble_id) const {
    return superbubbles[superbubble_id].chain_ids;
}

inline uint64_t SuperbubbleTree::chain_containing(uint64_t superbubble_id) const {
    return superbubbles[superbubble_id].parent;
}

inline const std::vector<uint64_t>& SuperbubbleTree::superbubbles_inside(uint64_t chain_id) const {
    return chains[chain_id].superbubble_ids;
}

inline uint64_t SuperbubbleTree::superbubble_containing(uint64_t chain_id) const {
    return chains[chain_id].parent;
}

}

#endif /* centrolign_superbubbles_hpp */
