#ifndef centrolign_determinize_hpp
#define centrolign_determinize_hpp

#include <vector>
#include <map>
#include <algorithm>
#include <utility>
#include <iostream>

#include "centrolign/graph.hpp"
#include "centrolign/utility.hpp"
#include "centrolign/topological_order.hpp"
#include "centrolign/modify_graph.hpp"

namespace centrolign {


// create and return an equivalent base graph that is reverse deterministic
// TODO: in theory this can have exponential complexity, although it doesn't seem likely
// that it will crop up in practice
template<class BGraph>
BaseGraph determinize(const BGraph& graph);

SentinelTableau translate_tableau(const BaseGraph& determinized,
                                  const SentinelTableau& original_tableau);

// add paths from the original graph to the determinized one
// note: tableau should be the translated one
template<class BGraph>
void rewalk_paths(BaseGraph& determinized,
                  const SentinelTableau& translated_tableau,
                  const BGraph& graph);

/*
 * Template implementations
 */

template<class BGraph>
BaseGraph determinize(const BGraph& graph) {
    
    static bool debug_determinize = false;
    
    // make a map of node ID to topological index
    std::vector<uint64_t> top_index = invert(topological_order(graph));
    
    // queue consists of unique sets of node IDs, bucketed by the highest topological
    // index in the set. each set is mapped to the IDs (in the new graph) that it must have an
    // edge to
    // note: the innermost vectors (the keys) are sorted and deduplicated before being inserted
    std::vector<std::map<std::vector<uint64_t>, std::vector<uint64_t>>> queue(top_index.size());
    
    // add sinks into the queue
    for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
        if (graph.next_size(node_id) == 0) {
            size_t index = top_index[node_id];
            queue[index].emplace(std::vector<uint64_t>(1, node_id), std::vector<uint64_t>());
            if (debug_determinize) {
                std::cerr << "initializing queue with sink node " << node_id << " at index " << index << '\n';
            }
        }
    }
    
    
    BaseGraph determinized;
    // work back to front along topological order
    for (int64_t i = queue.size() - 1; i >= 0; --i) {
        if (debug_determinize) {
            std::cerr << "dequeueing sets with final topological index " << i << '\n';
        }
        for (const auto& record : queue[i]) {
            const std::vector<uint64_t>& node_set = record.first;
            const std::vector<uint64_t>& successors = record.second;
            if (debug_determinize) {
                std::cerr << "node set:\n";
                for (auto nid : node_set) {
                    std::cerr << ' ' << nid;
                }
                std::cerr << "\nsuccessors:\n";
                for (auto nid : successors) {
                    std::cerr << ' ' << nid;
                }
                std::cerr << '\n';
            }
            
            // make the new node for this set and add its successors
            uint64_t new_node = determinized.add_node(graph.label(node_set.front()));
            for (uint64_t succ : successors) {
                determinized.add_edge(new_node, succ);
            }
            
            // find all predecessors with equivalent characters
            std::map<char, std::vector<uint64_t>> predecessors;
            for (uint64_t node_id : node_set) {
                for (uint64_t prev_id : graph.previous(node_id)) {
                    predecessors[graph.label(prev_id)].push_back(prev_id);
                }
            }
            
            if (debug_determinize) {
                std::cerr << "predecessor sets:\n";
                
            }
            for (std::pair<const char, std::vector<uint64_t>>& pred_record : predecessors) {

                // sort and deduplicate the character group
                auto& pred_group = pred_record.second;
                std::sort(pred_group.begin(), pred_group.end());
                auto last_unique = std::unique(pred_group.begin(), pred_group.end());
                pred_group.resize(last_unique - pred_group.begin());
                
                if (debug_determinize) {
                    std::cerr << pred_record.first << ':';
                    for (auto nid : pred_record.second) {
                        std::cerr << ' ' << nid;
                    }
                    std::cerr << '\n';
                }
                
                // record the predecessor node set in the queue
                size_t max_index = 0;
                for (uint64_t node_id : pred_group) {
                    max_index = std::max<size_t>(top_index[node_id], max_index);
                }
                queue[max_index][std::move(pred_group)].push_back(new_node);
            }
        }
        // clean up the data as we go to keep memory use down
        queue[i].clear();
    }
    
    if (debug_determinize) {
        std::cerr << "finished determinizing:\n";
        print_graph(determinized, std::cerr);
    }
    
    return determinized;
}


template<class BGraph>
void rewalk_paths(BaseGraph& determinized,
                  const SentinelTableau& tableau,
                  const BGraph& graph) {
    
    static const bool debug = false;
    
    for (uint64_t path_id = 0; path_id < graph.path_size(); ++path_id) {
        
        if (debug) {
            std::cerr << "rewalk path " << path_id << '\n';
        }
        
        std::vector<uint64_t> translated_path;
        
        // walk backward from the component's sink (takes advantage of reverse determinism)
        uint64_t here = tableau.snk_id;
        for (uint64_t step_id : ReverseForEachAdapter<std::vector<uint64_t>>(graph.path(path_id))) {
            char base = graph.label(step_id);
            if (debug) {
                std::cerr << "look for predecessor of " << here << " with label " << base << '\n';
            }
            for (uint64_t prev_id : determinized.previous(here)) {
                if (determinized.label(prev_id) == base) {
                    if (debug) {
                        std::cerr << "found predecessor at " << prev_id << '\n';
                    }
                    translated_path.push_back(prev_id);
                    here = prev_id;
                    break;
                }
            }
        }
        
        // add the translated path into the determinized graph
        uint64_t new_path_id = determinized.add_path(graph.path_name(path_id));
        for (uint64_t node_id : ReverseForEachAdapter<std::vector<uint64_t>>(translated_path)) {
            determinized.extend_path(new_path_id, node_id);
        }
    }
    
}

}

#endif /* centrolign_determinize_hpp */
