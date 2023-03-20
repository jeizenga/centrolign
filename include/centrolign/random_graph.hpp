#ifndef centrolign_random_graph_hpp
#define centrolign_random_graph_hpp

#include <random>
#include <algorithm>
#include <vector>
#include <utility>
#include <string>

#include "centrolign/graph.hpp"

namespace centrolign {

template <class Generator>
BaseGraph random_graph(size_t num_nodes, size_t num_edges,
                       Generator& gen) {
    
    BaseGraph graph;
    
    std::string alphabet = "ACGT";
    std::uniform_int_distribution<size_t> base_distr(0, alphabet.size() - 1);
    
    for (size_t i = 0; i < num_nodes; ++i) {
        graph.add_node(alphabet[base_distr(gen)]);
    }
    
    
    std::vector<std::pair<uint64_t, uint64_t>> edges;
    for (uint64_t node_from = 0; node_from < graph.node_size(); ++node_from) {
        for (uint64_t node_to = node_from + 1; node_to < graph.node_size(); ++node_to) {
            edges.emplace_back(node_from, node_to);
        }
    }
    
    std::shuffle(edges.begin(), edges.end(), gen);
    
    for (size_t i = 0; i < num_edges && i < edges.size(); ++i) {
        graph.add_edge(edges[i].first, edges[i].second);
    }
    
    return graph;
}

}

#endif /* centrolign_random_graph_hpp */
