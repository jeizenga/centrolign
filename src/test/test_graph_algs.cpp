#include <cstdio>
#include <cstdlib>
#include <cstdlib>
#include <unordered_set>
#include <unordered_map>
#include <set>

#include "centrolign/utility.hpp"
#include "centrolign/random_graph.hpp"
#include "centrolign/determinize.hpp"
#include "centrolign/topological_order.hpp"
#include "centrolign/antichain_partition.hpp"
#include "centrolign/reverse_graph.hpp"
#include "centrolign/subgraph_extraction.hpp"

using namespace std;
using namespace centrolign;

// DFS
bool is_reachable(const BaseGraph& graph, uint64_t id_from, uint64_t id_to) {
    
    vector<bool> traversed(graph.node_size(), false);
    vector<uint64_t> stack(1, id_from);
    traversed[id_from] = true;
    while (!stack.empty()) {
        auto here = stack.back();
        stack.pop_back();
        for (auto nid : graph.next(here)) {
            if (nid == id_to) {
                return true;
            }
            if (!traversed[nid]) {
                traversed[nid] = true;
                stack.push_back(nid);
            }
        }
    }
    return false;
}

ConnectingGraphInfo ugly_extract_connecting_graph(const BaseGraph& graph,
                                                  uint64_t from_id, uint64_t to_id) {
    
    ConnectingGraphInfo to_return;
    for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
        if (is_reachable(graph, from_id, node_id) && is_reachable(graph, node_id, to_id)) {
            to_return.subgraph.add_node(graph.label(node_id));
            to_return.back_translation.push_back(node_id);
        }
    }

    unordered_map<uint64_t, uint64_t> fwd_translation;
    for (uint64_t n = 0; n < to_return.back_translation.size(); ++n) {
        fwd_translation[to_return.back_translation[n]] = n;
    }
    
    for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
        if (!fwd_translation.count(node_id)) {
            continue;
        }
        for (uint64_t next_id : graph.next(node_id)) {
            if (!fwd_translation.count(next_id)) {
                continue;
            }
            
            to_return.subgraph.add_edge(fwd_translation[node_id],
                                        fwd_translation[next_id]);
        }
    }
    
    for (uint64_t next_id : graph.next(from_id)) {
        if (fwd_translation.count(next_id)) {
            to_return.sources.push_back(fwd_translation[next_id]);
        }
    }
    for (uint64_t prev_id : graph.previous(to_id)) {
        if (fwd_translation.count(prev_id)) {
            to_return.sinks.push_back(fwd_translation[prev_id]);
        }
    }
    
    return to_return;
}

bool subgraphs_are_equivalent(const BaseGraph& subgraph1, const BaseGraph& subgraph2,
                              const std::vector<uint64_t>& back_translation1,
                              const std::vector<uint64_t>& back_translation2) {
    // same number of nodes
    if (subgraph1.node_size() != subgraph2.node_size()) {
        return false;
    }
    
    size_t num_edges1 = 0, num_edges2 = 0;
    unordered_map<uint64_t, uint64_t> fwd_translation1, fwd_translation2;
    for (uint64_t n = 0; n < subgraph1.node_size(); ++n) {
        fwd_translation1[back_translation1[n]] = n;
        fwd_translation2[back_translation2[n]] = n;
        num_edges1 += subgraph1.next_size(n);
        num_edges2 += subgraph2.next_size(n);
    }
    
    if (num_edges1 != num_edges2) {
        return false;
    }
    
    // exact same nodes
    for (auto r : fwd_translation1) {
        if (!fwd_translation2.count(r.first)) {
            return false;
        }
        if (subgraph1.label(r.second) != subgraph2.label(fwd_translation2[r.first])) {
            return false;
        }
    }
    
    set<pair<uint64_t, uint64_t>> edges1, edges2;
    for (uint64_t n = 0; n < subgraph1.node_size(); ++n) {
        for (auto i : subgraph1.next(n)) {
            edges1.emplace(back_translation1[n], back_translation1[i]);
        }
        for (auto i : subgraph2.next(n)) {
            edges2.emplace(back_translation2[n], back_translation2[i]);
        }
    }
    if (edges1 != edges2) {
        return false;
    }
    
    return true;
}

bool node_sets_are_equivalent(const std::vector<uint64_t>& node_set1,
                              const std::vector<uint64_t>& node_set2,
                              const std::vector<uint64_t>& back_translation1,
                              const std::vector<uint64_t>& back_translation2) {
    set<uint64_t> back_set1, back_set2;
    for (auto n : node_set1) {
        back_set1.insert(back_translation1[n]);
    }
    for (auto n : node_set2) {
        back_set2.insert(back_translation2[n]);
    }
    return node_set1.size() == node_set2.size() && back_set1 == back_set2;
}

bool is_reverse_deterministic(const BaseGraph& graph) {
    for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
        unordered_set<char> prev_bases;
        for (uint64_t prev_id : graph.previous(node_id)) {
            char base = graph.label(prev_id);
            if (prev_bases.count(base)) {
                return false;
            }
            prev_bases.insert(base);
        }
    }
    return true;
}

bool is_probably_equivalent(const BaseGraph& graph,
                            const BaseGraph& determinized,
                            default_random_engine& gen) {
    
    auto select_random_path = [&](const BaseGraph& graph) {
        
        std::vector<uint64_t> sources;
        for (uint64_t nid = 0; nid < graph.node_size(); ++nid) {
            if (graph.previous_size(nid) == 0) {
                sources.push_back(nid);
            }
        }
        
        uint64_t here = sources[uniform_int_distribution<size_t>(0, sources.size() - 1)(gen)];
        std::string sequence(1, graph.label(here));
        while (graph.next_size(here) != 0) {
            vector<uint64_t> nexts = graph.next(here);
            size_t i = uniform_int_distribution<size_t>(0, nexts.size() - 1)(gen);
            here = nexts[i];
            sequence.push_back(graph.label(here));
        }
        return sequence;
    };
    
    
    
    auto find_path = [&](const BaseGraph& graph, string& seq) {
        std::vector<uint64_t> sources;
        for (uint64_t nid = 0; nid < graph.node_size(); ++nid) {
            if (graph.previous_size(nid) == 0 && graph.label(nid) == seq.front()) {
                sources.push_back(nid);
            }
        }
        
        bool found_path = false;
        for (auto start_id : sources) {
            
            vector<pair<uint64_t, size_t>> stack(1, make_pair(start_id, (size_t) 0));
            while (!stack.empty()) {
                size_t i;
                uint64_t nid;
                tie(nid, i) = stack.back();
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
        string seq = select_random_path(graph);
        if (!find_path(determinized, seq)) {
            cerr << "failed to find sequence " << seq << " from original graph:\n";
            print_graph(graph, cerr);
            cerr << "in determinized graph:\n";
            print_graph(determinized, cerr);
            exit(1);
        }
    }
    
    return true;
}


void test_topological_order(const BaseGraph& graph) {
    auto order = topological_order(graph);
    vector<size_t> index(order.size());
    for (size_t i = 0; i < order.size(); ++i) {
        index[order[i]] = i;
    }
    for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
        for (auto next_id : graph.next(node_id)) {
            if (index[node_id] >= index[next_id]) {
                std::cerr << "topological ordering failed on edge " << node_id << " -> " << next_id << " on graph:\n";
                print_graph(graph, cerr);
                exit(1);
            }
        }
        for (auto prev_id : graph.previous(node_id)) {
            if (index[prev_id] >= index[node_id]) {
                std::cerr << "topological ordering failed on edge " << prev_id << " -> " << node_id << " on graph:\n";
                print_graph(graph, cerr);
                exit(1);
            }
        }
    }
}

void test_subgraph_extraction(const BaseGraph& graph,
                              default_random_engine& gen) {
    
    ChainMerge chain_merge(graph);
    
    uniform_int_distribution<uint64_t> node_distr(0, graph.node_size() - 1);
    
    auto dump_connecting_info = [&](const ConnectingGraphInfo& info) {
        cerr << "subgraph:\n";
        print_graph(info.subgraph, cerr);
        cerr << "translation:\n";
        for (size_t i = 0; i < info.back_translation.size(); ++i) {
            cerr << " " << i << ":" << info.back_translation[i] << '\n';
        }
        cerr << "sources:";
        for (auto s : info.sources) {
            cerr << ' ' << s;
        }
        cerr << '\n';
        cerr << "sinks:";
        for (auto s : info.sinks) {
            cerr << ' ' << s;
        }
        cerr << '\n';
    };
    
    for (size_t rep = 0; rep < 10; ++rep) {
        uint64_t from = node_distr(gen);
        uint64_t to = node_distr(gen);
        
        auto extracted = extract_connecting_graph(graph, 2, 0, chain_merge);
        auto expected = ugly_extract_connecting_graph(graph, 2, 0);
        
        if (!subgraphs_are_equivalent(extracted.subgraph, expected.subgraph,
                                      extracted.back_translation, expected.back_translation) ||
            !node_sets_are_equivalent(extracted.sources, expected.sources,
                                      extracted.back_translation, expected.back_translation) ||
            !node_sets_are_equivalent(extracted.sinks, expected.sinks,
                                      extracted.back_translation, expected.back_translation)) {
            cerr << "connecting graph extraction failed betwen nodes " << from << " and " << to << " on graph:\n";
            print_graph(graph, cerr);
            cerr << "fast implementation:\n";
            dump_connecting_info(extracted);
            cerr << "slow implementation:\n";
            dump_connecting_info(expected);
            exit(1);
        }
    }
    
}

void test_antichain_partition(const BaseGraph& graph) {

    auto partition = antichain_partition(graph);
    std::vector<std::vector<uint64_t>> partition_sets;
    for (size_t node_id = 0; node_id < graph.node_size(); ++node_id) {
        size_t set = partition[node_id];
        while (partition_sets.size() <= set) {
            partition_sets.emplace_back();
        }
        partition_sets[set].push_back(node_id);
    }

    for (size_t i = 0; i < partition_sets.size(); ++i) {
        for (size_t j = i; j < partition_sets.size(); ++j) {
            for (uint64_t node_id : partition_sets[i]) {
                for (uint64_t other_id : partition_sets[j]) {
                    if (is_reachable(graph, other_id, node_id)) {
                        cerr << "unreachable antichains " << i << " and " << j << " contain comparable IDs " << node_id << " and " << other_id << " in graph:\n";
                        print_graph(graph, cerr);
                        exit(1);
                    }
                }
            }
        }
        
    }
}

void do_tests(const BaseGraph& graph, default_random_engine& gen) {
    
    BaseGraph determinized = determinize(graph);
    if (!is_reverse_deterministic(determinized)) {
        cerr << "determinization failed on graph:\n";
        print_graph(graph, cerr);
        cerr << "resulted in determinized graph:\n";
        print_graph(determinized, cerr);
        exit(1);
    }
    is_probably_equivalent(graph, determinized, gen);
    
    test_topological_order(graph);
    test_topological_order(determinized);
    
    test_antichain_partition(graph);
    test_antichain_partition(determinized);
    
    test_subgraph_extraction(graph, gen);
    //test_subgraph_extraction(determinized, gen); // it doesn't have a path cover
}

void add_sentinels(BaseGraph& graph) {
    assert(graph.node_size() != 0);
    uint64_t src_id = graph.add_node('#');
    uint64_t snk_id = graph.add_node('$');
    for (uint64_t node_id = 0; node_id + 2 < graph.node_size(); ++node_id) {
        if (graph.previous_size(node_id) == 0) {
            graph.add_edge(src_id, node_id);
        }
        if (graph.next_size(node_id) == 0) {
            graph.add_edge(node_id, snk_id);
        }
    }
}

int main(int argc, char* argv[]) {
     
    random_device rd;
    default_random_engine gen(rd());
    
    BaseGraph graph1;
    string graph1_labels = "GGGAT";
    for (size_t i = 0; i < 5; ++i) {
        graph1.add_node(graph1_labels[i]);
        for (size_t j = 0; j < i; ++j) {
            graph1.add_edge(j, i);
        }
    }
    // we have to have sentinels to make the determinize algorithm identify
    // the initial position of a sequence
    add_sentinels(graph1);
    add_random_path_cover(graph1, gen);
    do_tests(graph1, gen);
        
    size_t num_reps = 10;
    vector<pair<size_t, size_t>> graph_sizes;
    graph_sizes.emplace_back(8, 15);
    graph_sizes.emplace_back(10, 12);
    graph_sizes.emplace_back(20, 35);
    for (auto& sizes : graph_sizes) {
        size_t num_nodes = sizes.first;
        size_t num_edges = sizes.second;
        for (size_t i = 0; i < num_reps; ++i) {
            BaseGraph graph = random_graph(num_nodes, num_edges, gen);
            add_sentinels(graph);
            add_random_path_cover(graph, gen);
            do_tests(graph, gen);
        }
    }
    
    cerr << "passed all tests!" << endl;
}
