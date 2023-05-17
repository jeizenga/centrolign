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

// find superbubbles with algorithm of Gartner, et al. (2018) for DAGs
template<class Graph>
std::vector<std::pair<uint64_t, uint64_t>> find_superbubbles(const Graph& graph);




/*
 * Template instantiations
 */

template<class Graph>
std::vector<std::pair<uint64_t, uint64_t>> find_superbubbles(const Graph& graph) {
    
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

}

#endif /* centrolign_superbubbles_hpp */
