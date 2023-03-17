#ifndef centrolign_determinize_hpp
#define centrolign_determinize_hpp

#include <vector>
#include <map>
#include <algorithm>
#include <utility>

#include "centrolign/graph.hpp"
#include "centrolign/topological_order.hpp"

namespace centrolign {


// create and return an equivalent base graph that is reverse deterministic
// TODO: in theory this can have exponential complexity, although it doesn't seem likely
// that it will crop up in practice
template<class BGraph>
BaseGraph determinize(const BGraph& graph) {
    
    // make a map of node ID to topological index
    std::vector<size_t> top_index(graph.node_size());
    {
        std::vector<uint64_t> top_order = topological_order(graph);
        for (size_t i = 0; i < top_index.size(); ++i) {
            top_index[top_order[i]] = i;
        }
    }
    
    // queue consists of unique sets of node IDs, bucketed by the highest topological
    // index in the set. each set is mapped to the IDs (in the new graph) that it must have an
    // edge to
    // note: the innermost vectors (the keys) are sorted and deduplicated before being inserted
    std::vector<std::map<std::vector<uint64_t>, std::vector<uint64_t>>> queue(top_index.size());
    
    // add sinks into the queue
    for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
        if (graph.next_size(node_id) == 0) {
            size_t index = top_index[node_id];
            queue[index].emplace(1, node_id);
        }
    }
    
    BaseGraph determinized;
    // work back to front along topological order
    for (int64_t i = queue.size() - 1; i >= 0; --i) {
        for (const auto& record : queue[i]) {
            const std::vector<uint64_t>& node_set = record.first;
            const std::vector<uint64_t>& successors = record.second;
            
            // make the new node for this set and add its successors
            uint64_t new_node = determinized.add_node(graph.base(node_set.front()));
            for (uint64_t succ : successors) {
                determinized.add_edge(new_node, succ);
            }
            
            // find all predecessors with equivalent characters
            std::map<char, std::vector<uint64_t>> predecessors;
            for (uint64_t node_id : node_set) {
                for (uint64_t prev_id : graph.previous(node_id)) {
                    predecessors[graph.base(prev_id)].push_back(prev_id);
                }
            }
            
            for (std::pair<const char, std::vector<uint64_t>>& pred_record : predecessors) {
                // sort and deduplicate the character group
                auto& pred_group = pred_record.second;
                std::sort(pred_group.begin(), pred_group.end());
                auto last_unique = std::unique(pred_group.begin(), pred_group.end());
                pred_group.resize(last_unique - pred_group.begin());
                
                // record the predecessor node set in the queue
                size_t max_index = 0;
                for (uint64_t node_id : pred_group) {
                    max_index = std::max(top_index[node_id], max_index);
                }
                queue[max_index][std::move(pred_group)].push_back(new_node);
            }
        }
        // clean up the data as we go to keep memory use down
        queue[i].clear();
    }
    
    return determinized;
}

}

#endif /* centrolign_determinize_hpp */
