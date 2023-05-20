#ifndef centrolign_test_util_hpp
#define centrolign_test_util_hpp

#include <random>
#include <algorithm>
#include <vector>
#include <utility>
#include <string>
#include <list>

#include "centrolign/graph.hpp"
#include "centrolign/alignment.hpp"

namespace centrolign {

template <class Generator>
BaseGraph random_graph(size_t num_nodes, size_t num_edges,
                       Generator& gen);


template <class PGraph, class Generator>
void add_random_path_cover(PGraph& graph, Generator& gen);

template <class Generator>
std::string random_sequence(size_t length, Generator& gen);

template <class Generator>
std::string mutate_sequence(const std::string& seq, double sub_rate, double indel_rate,
                            Generator& gen);

bool is_reachable(const BaseGraph& graph, uint64_t id_from, uint64_t id_to);

std::vector<std::vector<uint64_t>> all_paths(const BaseGraph& graph,
                                             uint64_t id_from, uint64_t id_to);

std::string path_to_string(const BaseGraph& graph, const std::vector<uint64_t>& path);

// identical node IDs
bool graphs_are_identical(const BaseGraph& graph1, const BaseGraph& graph2);

// identical back-translated node IDs
bool translated_graphs_are_identical(const BaseGraph& subgraph1, const BaseGraph& subgraph2,
                                     const std::vector<uint64_t>& back_translation1,
                                     const std::vector<uint64_t>& back_translation2);

// isomorphism is hard, but we can check many simpler conditions that rule it out
bool possibly_isomorphic(const BaseGraph& graph1, const BaseGraph& graph2);

// equivalent in the sense of accepting the same string set
template<class Generator>
bool is_probably_equivalent(const BaseGraph& graph1, const BaseGraph& graph2,
                            Generator& gen);

void check_alignment(const Alignment& got, const Alignment& expected);

/*
 * Template implementations
 */

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

template <class PGraph, class Generator>
void add_random_path_cover(PGraph& graph, Generator& gen) {
    
    uint64_t num_covered = 0;
    std::vector<bool> covered(graph.node_size(), false);
    
    while (num_covered < graph.node_size()) {
        size_t rank = std::uniform_int_distribution<size_t>(0, graph.node_size() - num_covered - 1)(gen);
        uint64_t seed = -1;
        for (uint64_t node_id = 0, r = 0; node_id < graph.node_size(); ++node_id) {
            if (!covered[node_id]) {
                if (r == rank) {
                    seed = node_id;
                    break;
                }
                ++r;
            }
        }
        assert(seed != -1);
        
        std::list<uint64_t> path{seed};
        while (graph.previous_size(path.front()) != 0) {
            std::vector<uint64_t> uncovered_prevs;
            for (auto p : graph.previous(path.front())) {
                if (!covered[p]) {
                    uncovered_prevs.push_back(p);
                }
            }
            uint64_t nid;
            if (uncovered_prevs.empty()) {
                size_t idx = std::uniform_int_distribution<size_t>(0, graph.previous_size(path.front()) - 1)(gen);
                nid = graph.previous(path.front())[idx];
            }
            else {
                size_t idx = std::uniform_int_distribution<size_t>(0, uncovered_prevs.size() - 1)(gen);
                nid = uncovered_prevs[idx];
            }
            path.push_front(nid);
        }
        while (graph.next_size(path.back()) != 0) {
            std::vector<uint64_t> uncovered_nexts;
            for (auto n : graph.next(path.back())) {
                if (!covered[n]) {
                    uncovered_nexts.push_back(n);
                }
            }
            uint64_t nid;
            if (uncovered_nexts.empty()) {
                size_t idx = std::uniform_int_distribution<size_t>(0, graph.next_size(path.back()) - 1)(gen);
                nid = graph.next(path.back())[idx];
            }
            else {
                size_t idx = std::uniform_int_distribution<size_t>(0, uncovered_nexts.size() - 1)(gen);
                nid = uncovered_nexts[idx];
            }
            path.push_back(nid);
        }
        
        uint64_t path_id = graph.add_path(std::to_string(num_covered));
        for (auto nid : path) {
            if (!covered[nid]) {
                covered[nid] = true;
                ++num_covered;
            }
            graph.extend_path(path_id, nid);
        }
    }
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

template<class Generator>
bool is_probably_equivalent(const BaseGraph& graph,
                            const BaseGraph& determinized,
                            Generator& gen) {
    
    auto select_random_path = [&](const BaseGraph& graph) {
        
        std::vector<uint64_t> sources;
        for (uint64_t nid = 0; nid < graph.node_size(); ++nid) {
            if (graph.previous_size(nid) == 0) {
                sources.push_back(nid);
            }
        }
        
        uint64_t here = sources[std::uniform_int_distribution<size_t>(0, sources.size() - 1)(gen)];
        std::string sequence(1, graph.label(here));
        while (graph.next_size(here) != 0) {
            std::vector<uint64_t> nexts = graph.next(here);
            size_t i = std::uniform_int_distribution<size_t>(0, nexts.size() - 1)(gen);
            here = nexts[i];
            sequence.push_back(graph.label(here));
        }
        return sequence;
    };
    
    
    
    auto find_path = [&](const BaseGraph& graph, std::string& seq) {
        std::vector<uint64_t> sources;
        for (uint64_t nid = 0; nid < graph.node_size(); ++nid) {
            if (graph.previous_size(nid) == 0 && graph.label(nid) == seq.front()) {
                sources.push_back(nid);
            }
        }
        
        bool found_path = false;
        for (auto start_id : sources) {
            
            std::vector<std::pair<uint64_t, size_t>> stack(1, std::make_pair(start_id, (size_t) 0));
            while (!stack.empty()) {
                size_t i;
                uint64_t nid;
                std::tie(nid, i) = stack.back();
                stack.pop_back();
                if (i + 1 == seq.size()) {
                    if (graph.next_size(nid) == 0) {
                        found_path = true;
                        break;
                    }
                }
                else {
                    for (uint64_t next_id : graph.next(nid)) {
                        if (graph.label(next_id) == seq[i + 1]) {
                            stack.emplace_back(next_id, i + 1);
                        }
                    }
                }
            }
            if (found_path) {
                break;
            }
        }
        return found_path;
    };
    
    size_t num_tests = 5;
    for (size_t i = 0; i < num_tests; ++i) {
        std::string seq = select_random_path(graph);
        if (!find_path(determinized, seq)) {
            std::cerr << "failed to find sequence " << seq << " from original graph:\n";
            print_graph(graph, std::cerr);
            std::cerr << "in determinized graph:\n";
            print_graph(determinized, std::cerr);
            return false;
        }
    }
    
    return true;
}

}

#endif /* centrolign_test_util_hpp */
