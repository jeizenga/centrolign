#ifndef centrolign_subset_graph_hpp
#define centrolign_subset_graph_hpp

#include "centrolign/graph.hpp"
#include "centrolign/modify_graph.hpp"

namespace centrolign {


BaseGraph subset_graph(const BaseGraph& graph, const std::vector<uint64_t>& path_subset);

std::pair<BaseGraph, SentinelTableau> subset_graph(const BaseGraph& graph, const SentinelTableau& tableau,
                                                   const std::vector<uint64_t>& path_subset);

}

#endif /* centrolign_subset_graph_hpp */
