#include "centrolign/modify_graph.hpp"

#include <unordered_set>
#include <vector>

namespace centrolign {

using namespace std;

SentinelTableau add_sentinels(BaseGraph& graph,
                              char src_sentinel, char snk_sentinel) {
    
    std::unordered_set<uint64_t> path_endings, path_beginnings;
    for (uint64_t path_id = 0; path_id < graph.path_size(); ++path_id) {
        path_beginnings.insert(graph.path(path_id).front());
        path_endings.insert(graph.path(path_id).back());
    }
    
    SentinelTableau tableau;
    tableau.src_id = graph.add_node(src_sentinel);
    tableau.snk_id = graph.add_node(snk_sentinel);
    tableau.src_sentinel = src_sentinel;
    tableau.snk_sentinel = snk_sentinel;
    
    if (graph.node_size() == 2) {
        graph.add_edge(tableau.src_id, tableau.snk_id);
    }
    else {
        for (uint64_t node_id = 0; node_id < tableau.src_id; ++node_id) {
            if (graph.next_size(node_id) == 0 || path_endings.count(node_id)) {
                graph.add_edge(node_id, tableau.snk_id);
            }
            if (graph.previous_size(node_id) == 0 || path_beginnings.count(node_id)) {
                graph.add_edge(tableau.src_id, node_id);
            }
        }
    }
    
    return tableau;
}

void reassign_sentinels(BaseGraph& graph, SentinelTableau& tableau,
                        char src_sentinel, char snk_sentinel) {
    
    tableau.src_sentinel = src_sentinel;
    tableau.snk_sentinel = snk_sentinel;
    graph.relabel(tableau.src_id, src_sentinel);
    graph.relabel(tableau.snk_id, snk_sentinel);
}


}
