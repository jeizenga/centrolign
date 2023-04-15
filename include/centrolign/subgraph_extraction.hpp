#ifndef centrolign_subgraph_extraction_hpp
#define centrolign_subgraph_extraction_hpp

#include <vector>
#include <iostream>
#include <unordered_map>
#include <utility>

#include "centrolign/graph.hpp"
#include "centrolign/chain_merge.hpp"

namespace centrolign {

struct ConnectingGraphInfo {
    // the actual subgraph
    BaseGraph subgraph;
    // the translation of subgraph node IDs to parent graph node IDs
    std::vector<uint64_t> back_translation;
    // the node IDs in the subgraph that were adjacent to the source ID
    std::vector<uint64_t> sources;
    // the node IDs in the subgraph that were adjacent to the sink ID
    std::vector<uint64_t> sinks;
    
    ConnectingGraphInfo() = default;
    ~ConnectingGraphInfo() = default;
};
// Extract the subgraph between two graph nodes, not including those nodes themselves.
template<class BGraph>
ConnectingGraphInfo extract_connecting_graph(const BGraph& graph,
                                             uint64_t from_id, uint64_t to_id,
                                             const ChainMerge& chain_merge) {
    
    static const bool debug = false;
    
    if (debug) {
        std::cerr << "extracting subgraph from " << from_id << " to " << to_id << '\n';
    }
    
    auto to_return = ConnectingGraphInfo();
    std::unordered_map<uint64_t, uint64_t> forward_translation;
    
    // DFS from one position
    std::vector<uint64_t> stack(1, from_id);
    while (!stack.empty()) {
        uint64_t node_id = stack.back();
        stack.pop_back();
        if (debug) {
            std::cerr << "destack " << node_id << '\n';
        }
        for (uint64_t next_id : graph.next(node_id)) {
            if (debug) {
                std::cerr << "follow edge to " << next_id << '\n';
            }
            if (next_id == to_id && node_id != from_id) {
                if (debug) {
                    std::cerr << "stop at the sink\n";
                }
                to_return.sinks.push_back(forward_translation[node_id]);
                continue;
            }
            if (!chain_merge.reachable(next_id, to_id)) {
                // this can't reach the end
                if (debug) {
                    std::cerr << "stop because not reachable\n";
                }
                continue;
            }
            auto it = forward_translation.find(next_id);
            if (it == forward_translation.end()) {
                // this is the first time we encountered this node, add it
                // to the subgraph and translations
                if (debug) {
                    std::cerr << "add as new node\n";
                }
                forward_translation[next_id] = to_return.subgraph.node_size();
                to_return.back_translation.push_back(next_id);
                to_return.subgraph.add_node(graph.label(next_id));
                // set up the iterator
                it = forward_translation.find(next_id);
                // and queue the node up
                stack.push_back(next_id);
            }
            
            if (node_id != from_id) {
                if (debug) {
                    std::cerr << "add edge\n";
                }
                // we need to add the edge from the predecessor
                to_return.subgraph.add_edge(forward_translation.at(node_id), it->second);
            }
            else {
                if (debug) {
                    std::cerr << "identify source\n";
                }
                to_return.sources.push_back(it->second);
            }
        }
    }
    
    return to_return;
}

}

#endif /* centrolign_subgraph_extraction_hpp */
