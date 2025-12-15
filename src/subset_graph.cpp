#include "centrolign/subset_graph.hpp"

#include <unordered_set>
#include <cassert>

namespace centrolign {


void subset_graph_internal(BaseGraph& subsetted, SentinelTableau* subset_tableau,
                           const BaseGraph& graph, const std::vector<uint64_t>& path_subset, const SentinelTableau* tableau) {
    
    // copy over the relevant nodes
    std::vector<uint64_t> translation(graph.node_size(), -1);
    for (uint64_t path_id : path_subset) {
        for (uint64_t node_id : graph.path(path_id)) {
            if (translation[node_id] == -1) {
                translation[node_id] = subsetted.add_node(graph.label(node_id));
            }
        }
    }
    
    if (tableau) {
        // add the source and sink nodes
        assert(subset_tableau);
        subset_tableau->src_id = subsetted.add_node(tableau->src_sentinel);
        subset_tableau->snk_id = subsetted.add_node(tableau->snk_sentinel);
        subset_tableau->src_sentinel = tableau->src_sentinel;
        subset_tableau->snk_sentinel = tableau->snk_sentinel;
    }
    
    // which edges have been created, for O(1) membership testing
    std::vector<std::unordered_set<uint64_t>> edges(subsetted.node_size());
    
    // copy over the edges and paths
    for (uint64_t path_id : path_subset) {
        uint64_t subset_path_id = subsetted.add_path(graph.path_name(path_id));
        
        uint64_t subset_prev_id = -1;
        for (uint64_t node_id : graph.path(path_id)) {
            subsetted.extend_path(subset_path_id, translation[node_id]);
            
            uint64_t subset_node_id = translation[node_id];
            if (subset_prev_id != -1 && !edges[subset_prev_id].count(subset_node_id)) {
                subsetted.add_edge(subset_prev_id, subset_node_id);
                edges[subset_prev_id].insert(subset_node_id);
            }
            subset_prev_id = subset_node_id;
        }
    }
    
    if (tableau) {
        // add edges to/from the source and sink nodes
        for (uint64_t path_id = 0; path_id < subsetted.path_size(); ++path_id) {
            const auto& path = subsetted.path(path_id);
            if (!edges[subset_tableau->src_id].count(path.front())) {
                subsetted.add_edge(subset_tableau->src_id, path.front());
                edges[subset_tableau->src_id].insert(path.front());
            }
            if (!edges[path.back()].count(subset_tableau->snk_id)) {
                subsetted.add_edge(path.back(), subset_tableau->snk_id);
                edges[path.back()].insert(subset_tableau->snk_id);
            }
        }
    }
}

BaseGraph subset_graph(const BaseGraph& graph, const std::vector<uint64_t>& path_subset) {
    
    BaseGraph subsetted;
    subset_graph_internal(subsetted, nullptr, graph, path_subset, nullptr);
    return subsetted;
}

std::pair<BaseGraph, SentinelTableau> subset_graph(const BaseGraph& graph, const SentinelTableau& tableau,
                                                   const std::vector<uint64_t>& path_subset) {
    
    std::pair<BaseGraph, SentinelTableau> subsetted;
    subset_graph_internal(subsetted.first, &subsetted.second, graph, path_subset, &tableau);
    return subsetted;
}


}
