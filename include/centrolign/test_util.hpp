#ifndef centrolign_test_util_hpp
#define centrolign_test_util_hpp

#include <random>
#include <algorithm>
#include <vector>
#include <utility>
#include <string>
#include <list>

#include "centrolign/alignment.hpp"
#include "centrolign/graph.hpp"
#include "centrolign/modify_graph.hpp"

namespace centrolign {

template <class Generator>
BaseGraph random_graph(size_t num_nodes, size_t num_edges,
                       Generator& gen);

// low complexity, bubble-rich graph with many indels
template <class Generator>
BaseGraph random_challenge_graph(size_t num_nodes, Generator& gen);

template <class PGraph, class Generator>
void add_random_path_cover(PGraph& graph, Generator& gen, const SentinelTableau* tableau = nullptr);

template <class Generator>
std::string random_sequence(size_t length, Generator& gen);

template <class Generator>
std::string random_low_entropy_sequence(size_t length, Generator& gen);

template <class Generator>
std::string mutate_sequence(const std::string& seq, double sub_rate, double indel_rate,
                            Generator& gen);

// check if only follows real edges
bool is_valid_path(const BaseGraph& graph, std::vector<uint64_t>& path);

bool is_reachable(const BaseGraph& graph, uint64_t id_from, uint64_t id_to);

// all paths to any sink node
std::vector<std::vector<uint64_t>> all_paths(const BaseGraph& graph, uint64_t id_from);

// all paths between two nodes
std::vector<std::vector<uint64_t>> all_paths(const BaseGraph& graph,
                                             uint64_t id_from, uint64_t id_to);

// identical node IDs
bool graphs_are_identical(const BaseGraph& graph1, const BaseGraph& graph2);

// identical back-translated node IDs (only valid for one-to-one translations)
bool translated_graphs_are_identical(const BaseGraph& subgraph1, const BaseGraph& subgraph2,
                                     const std::vector<uint64_t>& back_translation1,
                                     const std::vector<uint64_t>& back_translation2);

bool translations_possibly_consistent(const BaseGraph& subgraph1, const BaseGraph& subgraph2,
                                      const std::vector<uint64_t>& back_translation1,
                                      const std::vector<uint64_t>& back_translation2);

// isomorphism is hard, but we can check many simpler conditions that rule it out
bool possibly_isomorphic(const BaseGraph& graph1, const BaseGraph& graph2);

// equivalent in the sense of accepting the same string set
template<class Generator>
bool is_probably_equivalent(const BaseGraph& graph1, const BaseGraph& graph2,
                            Generator& gen);

void check_alignment(const Alignment& got, const Alignment& expected);

bool paths_match(const BaseGraph& graph1, const BaseGraph& graph2);

// c++ code for making the graph
std::string cpp_representation(const BaseGraph& graph, const std::string& name);

std::string pretty_alignment(const Alignment& aln, const std::string& seq1, const std::string& seq2);

std::string pretty_alignment(const Alignment& aln, const BaseGraph graph1, const BaseGraph& graph2);

/*
 * Template implementations
 */

template <class Generator>
BaseGraph random_graph(size_t num_nodes, size_t num_edges, bool acyclic,
                       Generator& gen) {
    
    BaseGraph graph;
    
    std::string alphabet = "ACGT";
    std::uniform_int_distribution<size_t> base_distr(0, alphabet.size() - 1);
    
    for (size_t i = 0; i < num_nodes; ++i) {
        graph.add_node(alphabet[base_distr(gen)]);
    }
    
    
    std::vector<std::pair<uint64_t, uint64_t>> edges;
    for (uint64_t node_from = 0; node_from < graph.node_size(); ++node_from) {
        for (uint64_t node_to = (acyclic ? node_from + 1 : 0); node_to < graph.node_size(); ++node_to) {
            if (node_from == node_to) {
                continue;
            }
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
BaseGraph random_challenge_graph(size_t num_nodes, Generator& gen) {
    
    size_t repeat_length = std::max<size_t>(num_nodes / 6, 1);
    
    std::uniform_real_distribution<double> prob_distr(0.0, 1.0);
    std::string alphabet = "ACGT";
    std::uniform_int_distribution<size_t> alphabet_distr(0, 3);
    
    std::string base_repeat;
    bool is_A = true;
    for (size_t i = 0; i < repeat_length; ++i) {
        if (prob_distr(gen) < .33) {
            is_A = !is_A;
        }
        base_repeat.push_back(is_A ? 'A' : 'C');
    }
    
    double mut_rate = 0.1;
    
    BaseGraph graph;
    uint64_t prev_id = -1;
    for (size_t i = 0; i < 5; ++i) {
        std::string repeat = base_repeat;
        for (size_t j = 0; j < repeat.size(); ++j) {
            if (prob_distr(gen) < mut_rate) {
                repeat[j] = alphabet[alphabet_distr(gen)];
            }
        }
        for (auto c : repeat) {
            uint64_t node_id = graph.add_node(c);
            if (prev_id != -1) {
                graph.add_edge(prev_id, node_id);
            }
            prev_id = node_id;
        }
    }
    
    size_t base_size = graph.node_size();
    
    std::uniform_int_distribution<uint64_t> base_node_distr(0, base_size - 1);
    
    while (graph.node_size() < num_nodes) {
        
        uint64_t node_id = base_node_distr(gen);
        
        uint64_t sub_id = graph.add_node(alphabet[alphabet_distr(gen)]);
        for (auto p : graph.previous(node_id)) {
            graph.add_edge(p, sub_id);
        }
        for (auto n : graph.next(node_id)) {
            graph.add_edge(sub_id, n);
        }
    }
    
    std::uniform_int_distribution<uint64_t> all_node_distr(0, graph.node_size() - 1);
    
    std::uniform_int_distribution<size_t> num_small_del_distr(1, 5);
    
    size_t num_small_del = num_small_del_distr(gen);
    for (size_t i = 0; i < num_small_del; ++i) {
        uint64_t node_id = all_node_distr(gen);
        for (auto p : graph.previous(node_id)) {
            for (auto n : graph.next(node_id)) {
                bool found = false;
                for (auto d : graph.next(p)) {
                    if (d == n) {
                        found = true;
                    }
                }
                if (!found) {
                    graph.add_edge(p, n);
                }
            }
        }
    }
    
    std::uniform_int_distribution<size_t> num_large_del_distr(2, 4);
    size_t num_large_del = num_large_del_distr(gen);
    
    std::uniform_int_distribution<size_t> large_del_len_distr(base_repeat.size() > 2 ? base_repeat.size() - 2 : 1,
                                                              std::min<size_t>(base_repeat.size() + 2, base_size - 1));
    
    size_t complete_dels = 0;
    while (complete_dels < num_large_del) {
        
        uint64_t node_id = all_node_distr(gen);
        size_t len = large_del_len_distr(gen);
        
        auto here = node_id;
        size_t i = 0;
        for (; i < len; ++i) {
            if (graph.next_size(here) == 0) {
                break;
            }
            std::uniform_int_distribution<size_t> edge_distr(0, graph.next_size(here) - 1);
            here = graph.next(here)[edge_distr(gen)];
        }
        
        if (i == len) {
            graph.add_edge(node_id, here);
            ++complete_dels;
        }
    }
    
    return graph;
}

template <class PGraph, class Generator>
void add_random_path_cover(PGraph& graph, Generator& gen, const SentinelTableau* tableau) {
    
    uint64_t num_covered = 0;
    std::vector<bool> covered(graph.node_size(), false);
    if (tableau) {
        covered[tableau->src_id] = true;
        covered[tableau->snk_id] = true;
        num_covered = 2;
    }
    
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
            if (tableau && nid == tableau->src_id) {
                break;
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
            if (tableau && nid == tableau->snk_id) {
                break;
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
std::string random_low_entropy_sequence(size_t length, Generator& gen) {
    
    std::string alphabet = "ACGT";
    
    std::vector<std::vector<double>> conditional_distrs;
    
    std::gamma_distribution<double> gam_distr(0.75, 1.0);
    
    for (size_t i = 0; i < 4; ++i) {
        
        // sample from a dirichlet
        conditional_distrs.emplace_back();
        auto& con_distr = conditional_distrs.back();
        double total = 0.0;
        for (size_t j = 0; j < 4; ++j) {
            con_distr.emplace_back(gam_distr(gen));
            total += con_distr.back();
        }
        for (size_t j = 0; j < 4; ++j) {
            con_distr[j] /= total;
            if (j > 0) {
                con_distr[j] += con_distr[j - 1];
            }
        }
    }
    
    int state = std::uniform_int_distribution<int>(0, 3)(gen);
    
    std::uniform_real_distribution<double> prob_distr(0.0, 1.0);
    
    std::string seq;
    while (seq.size() < length) {
        
        double p = prob_distr(gen);
        
        int next_state = 0;
        while (conditional_distrs[state][next_state] < p) {
            ++next_state;
        }
        
        seq.push_back(alphabet[next_state]);
    }
    
    return seq;
}

template <class Generator>
std::string mutate_sequence(const std::string& seq, double sub_rate, double indel_rate,
                            Generator& gen) {
    
    std::string alphabet = "ACGT";
    std::uniform_int_distribution<size_t> base_distr(0, alphabet.size() - 1);
    std::uniform_real_distribution<double> prob_distr(0.0, 1.0);
    
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
            std::cerr << "in putative equivalent graph:\n";
            print_graph(determinized, std::cerr);
            return false;
        }
    }
    
    return true;
}

}

#endif /* centrolign_test_util_hpp */
