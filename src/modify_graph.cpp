#include "centrolign/modify_graph.hpp"

#include <unordered_set>
#include <vector>
#include <cassert>

#include "centrolign/utility.hpp"


namespace centrolign {

using namespace std;

SequenceGraph make_sequence_graph(const std::string& name,
                                  const std::string& sequence) {
    assert(!sequence.empty());
    assert(!name.empty());
    SequenceGraph graph;
    uint64_t node_id = graph.add_node(encode_seq(sequence));
    uint64_t path_id = graph.add_path(name);
    graph.extend_path(path_id, node_id);
    return graph;
}

BaseGraph make_base_graph(const std::string& name,
                          const std::string& sequence) {
    assert(!sequence.empty());
    assert(!name.empty());
    BaseGraph graph;
    uint64_t prev_id = graph.add_node(encode_base(sequence[0]));
    uint64_t path_id = graph.add_path(name);
    graph.extend_path(path_id, prev_id);
    for (size_t i = 1; i < sequence.size(); ++i) {
        uint64_t node_id = graph.add_node(encode_base(sequence[i]));
        graph.add_edge(prev_id, node_id);
        graph.extend_path(path_id, node_id);
        prev_id = node_id;
    }
    return graph;
}

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
        
        // add a dummy path for the sentinels
        auto path_id = graph.add_path("__dummy__");
        graph.extend_path(path_id, tableau.src_id);
        graph.extend_path(path_id, tableau.snk_id);
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
        
        // add the sentinels to an arbitrary path
        if (graph.path_size() != 0) {
            graph.pre_extend_path(0, tableau.src_id);
            graph.extend_path(0, tableau.snk_id);
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
