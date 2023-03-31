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


template <class Generator>
std::string random_sequence(size_t length, Generator& gen) {
    
    std::string alphabet = "ACGT";
    std::uniform_int_distribution<size_t> base_distr(0, alphabet.size() - 1);
    
    std::string seq;
    for (size_t i = 0; i < length; ++i) {
        seq.push_back(alphabet[base_distr(gen)]);
    }
    
    return seq;
}

template <class Generator>
std::string mutate_sequence(const std::string& seq, double sub_rate, double indel_rate,
                            Generator& gen) {
    
    std::string alphabet = "ACGT";
    std::uniform_int_distribution<size_t> base_distr(0, alphabet.size() - 1);
    std::uniform_real_distribution<size_t> prob_distr(0.0, 1.0);
    
    std::string mutated;
    for (size_t i = 0; i < seq.size(); ++i) {
        double p = prob_distr(gen);
        if (p < sub_rate) {
            mutated.push_back(alphabet[base_distr(gen)]);
        }
        else if (p < sub_rate + indel_rate) {
            if (prob_distr(gen) < 0.5) {
                ++i;
            }
            else {
                mutated.push_back(alphabet[base_distr(gen)]);
                --i;
            }
        }
        else {
            mutated.push_back(seq[i]);
        }
    }
    
    return mutated;
}


}

#endif /* centrolign_random_graph_hpp */
