#include <cstdio>
#include <cstdlib>
#include <cstdlib>
#include <unordered_set>
#include <unordered_map>
#include <set>

#include "centrolign/utility.hpp"
#include "centrolign/test_util.hpp"
#include "centrolign/determinize.hpp"
#include "centrolign/topological_order.hpp"
#include "centrolign/antichain_partition.hpp"
#include "centrolign/reverse_graph.hpp"
#include "centrolign/subgraph_extraction.hpp"
#include "centrolign/count_walks.hpp"
#include "centrolign/minmax_distance.hpp"
#include "centrolign/fuse.hpp"
#include "centrolign/target_reachability.hpp"
#include "centrolign/shortest_path.hpp"
#include "centrolign/union_find.hpp"
#include "centrolign/connected_components.hpp"

using namespace std;
using namespace centrolign;

bool is_simple(const BaseGraph& graph) {
    
    for (uint64_t n = 0; n < graph.node_size(); ++n) {
        std::unordered_set<uint64_t> nexts;
        for (auto m : graph.next(n)) {
            if (nexts.count(m)) {
                return false;
            }
            nexts.insert(m);
        }
    }
    return true;
}

std::vector<std::vector<uint64_t>> dumb_connected_components(const BaseGraph& graph) {
    
    UnionFind uf(graph.node_size());
    
    for (uint64_t n1 = 0; n1 < graph.node_size(); ++n1) {
        for (uint64_t n2 = 0; n2 < graph.node_size(); ++n2) {
            if (n1 != n2) {
                if (is_reachable(graph, n1, n2)) {
                    uf.union_groups(n1, n2);
                }
            }
        }
    }
    
    // shim from size_t to uint64_t
    std::vector<std::vector<uint64_t>> components;
    for (auto& group : uf.get_groups()) {
        components.emplace_back(group.begin(), group.end());
    }
    return components;
}

std::vector<uint64_t> shortest_path_brute_force(const BaseGraph& graph,
                                                uint64_t from, uint64_t to) {
    auto paths = all_paths(graph, from, to);
    
    std::vector<uint64_t> shortest;
    for (auto& p : paths) {
        if (shortest.empty() || p.size() < shortest.size()) {
            shortest = move(p);
        }
    }
    return shortest;
}

uint64_t count_walks_brute_force(const BaseGraph& graph) {
    
    std::vector<uint64_t> sources;
    for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
        if (graph.previous_size(node_id) == 0) {
            sources.push_back(node_id);
        }
    }
    
    uint64_t total = 0;
    for (auto src_id : sources) {
        
        vector<pair<vector<uint64_t>, size_t>> stack;
        stack.emplace_back(graph.next(src_id), 0);
        
        while (!stack.empty()) {
            auto& top = stack.back();
            if (top.first.empty() && total != -1) {
                ++total;
            }
            if (top.second < top.first.size()) {
                stack.emplace_back(graph.next(top.first[top.second++]), 0);
            }
            else {
                stack.pop_back();
            }
        }
    }
    
    return total;
}

pair<int64_t, int64_t> minmax_distance_brute_force(const BaseGraph& graph, uint64_t id_from, uint64_t id_to) {
    
    pair<int64_t, int64_t> minmax(numeric_limits<int64_t>::max(), -1);
    auto paths = all_paths(graph, id_from, id_to);
    for (const auto& path : paths) {
        // remove 1 so we don't count the source node
        minmax.first = min<int64_t>(minmax.first, path.size() - 1);
        minmax.second = max<int64_t>(minmax.second, path.size() - 1);
    }
    
    return minmax;
}

void test_minmax_distance(const BaseGraph& graph, uint64_t src_id) {
        
    vector<uint64_t> src{src_id};
    auto got = minmax_distance(graph, &src);
    
    for (uint64_t n = 0; n < graph.node_size(); ++n) {
        if (n == src_id) {
            continue;
        }
        auto expected = minmax_distance_brute_force(graph, src_id, n);
        if (got[n] != expected) {
            cerr << "failed minimax distance test between " << src_id << " and " << n << '\n';
            cerr << "got: " << got[n].first << ", " << got[n].second << '\n';
            cerr << "expected: " << expected.first << ", " << expected.second << '\n';
            cerr << "full results:\n";
            for (auto g : got) {
                cerr << g.first << '\t' << g.second << '\n';
            }
            cerr << "graph:\n";
            cerr << cpp_representation(graph, "graph") << '\n';
            exit(1);
        }
    }
}


void test_count_walks(const BaseGraph& graph) {
    
    uint64_t expected = count_walks_brute_force(graph);
    uint64_t got = count_walks(graph);
    
    if (got != expected) {
        cerr << "count walks test failed with " << got << " rather than " << expected << " in graph:\n";
        print_graph(graph, cerr);
        exit(1);
    }
}

void fill_in_extracted_subgraph(const BaseGraph& graph, SubGraphInfo& info,
                                uint64_t from_id, uint64_t to_id) {
    
    unordered_map<uint64_t, uint64_t> fwd_translation;
    for (uint64_t n = 0; n < info.back_translation.size(); ++n) {
        fwd_translation[info.back_translation[n]] = n;
    }
    
    for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
        if (!fwd_translation.count(node_id)) {
            continue;
        }
        for (uint64_t next_id : graph.next(node_id)) {
            if (!fwd_translation.count(next_id)) {
                continue;
            }
            info.subgraph.add_edge(fwd_translation[node_id],
                                   fwd_translation[next_id]);
        }
    }
    
    if (from_id != -1) {
        for (uint64_t next_id : graph.next(from_id)) {
            if (fwd_translation.count(next_id)) {
                info.sources.push_back(fwd_translation[next_id]);
            }
        }
    }
    if (to_id != -1) {
        for (uint64_t prev_id : graph.previous(to_id)) {
            if (fwd_translation.count(prev_id)) {
                info.sinks.push_back(fwd_translation[prev_id]);
            }
        }
    }
}

SubGraphInfo ugly_extract_connecting_graph(const BaseGraph& graph,
                                           uint64_t from_id, uint64_t to_id) {
    
    SubGraphInfo to_return;
    for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
        if (is_reachable(graph, from_id, node_id) && is_reachable(graph, node_id, to_id)) {
            to_return.subgraph.add_node(graph.label(node_id));
            to_return.back_translation.push_back(node_id);
        }
    }

    fill_in_extracted_subgraph(graph, to_return, from_id, to_id);
    
    return to_return;
}

SubGraphInfo ugly_extract_extending_graph(const BaseGraph& graph,
                                          uint64_t from_id, bool forward) {

    SubGraphInfo to_return;
    for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
        if ((is_reachable(graph, from_id, node_id) && forward) ||
            (is_reachable(graph, node_id, from_id) && !forward)) {
            
            to_return.subgraph.add_node(graph.label(node_id));
            to_return.back_translation.push_back(node_id);
        }
    }
    
    fill_in_extracted_subgraph(graph, to_return, forward ? from_id : -1, forward ? -1 : from_id);
    
    return to_return;
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



void test_topological_order(const BaseGraph& graph) {
    auto order = topological_order(graph);
    auto index = invert(order);
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

void test_fuse(const BaseGraph& graph1, const BaseGraph& graph2,
               const SentinelTableau& tableau1, const SentinelTableau& tableau2,
               const Alignment& alignment, const BaseGraph& expected) {
    
    
    // copy to make them mutable
    BaseGraph destination = graph1;
    BaseGraph source = graph2;
    
    fuse(destination, source, tableau1, tableau2, alignment);
    
    if (!is_simple(destination)) {
        cerr << "fused graph is not simple\n";
        cerr << cpp_representation(graph1, "destination") << '\n';
        cerr << cpp_representation(graph2, "source") << '\n';
        exit(1);
    }
    if (!possibly_isomorphic(destination, expected)) {
        cerr << "failed fuse test\n";
        cerr << "graph1:\n";
        print_graph(graph1, cerr);
        cerr << "graph2:\n";
        print_graph(graph2, cerr);
        cerr << "alignment:\n";
        for (const auto& ap : alignment) {
            cerr << '\t' << (ap.node_id1 == AlignedPair::gap ? string("-") : to_string(ap.node_id1)) << '\t' << (ap.node_id2 == AlignedPair::gap ? string("-") : to_string(ap.node_id2)) << '\n';
        }
        cerr << "expected:\n";
        print_graph(expected, cerr);
        cerr << "got:\n";
        print_graph(destination, cerr);
        
        exit(1);
    }
}

void test_internal_fuse(const BaseGraph& graph, const vector<Alignment>& alignments,
                        const BaseGraph& expected) {
    
    BaseGraph got = internal_fuse(graph, alignments);
    
    if (!possibly_isomorphic(got, expected)) {
        
        cerr << "failed internal fuse test\n";
        cerr << "graph:\n";
        print_graph(graph, cerr);
        cerr << "alignment:\n";
        for (size_t i = 0; i < alignments.size(); ++i) {
            cerr << "alignment " << i << ":\n";
            for (const auto& ap : alignments[i]) {
                cerr << '\t' << (ap.node_id1 == AlignedPair::gap ? string("-") : to_string(ap.node_id1)) << '\t' << (ap.node_id2 == AlignedPair::gap ? string("-") : to_string(ap.node_id2)) << '\n';
            }
        }
        cerr << "expected:\n";
        print_graph(expected, cerr);
        cerr << "got:\n";
        print_graph(got, cerr);
        exit(1);
    }
    
}

void test_subgraph_extraction(const BaseGraph& graph, const SentinelTableau& tableau,
                              default_random_engine& gen) {
    
    ChainMerge chain_merge(graph, tableau);
    
    uniform_int_distribution<uint64_t> node_distr(0, graph.node_size() - 1);
    
    auto dump_connecting_info = [&](const SubGraphInfo& info) {
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
        {
            auto extracted = extract_connecting_graph(graph, from, to, chain_merge);
            auto expected = ugly_extract_connecting_graph(graph, from, to);
            
            if (!translated_graphs_are_identical(extracted.subgraph, expected.subgraph,
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
        
        for (auto node_id : {from, to}) {
            for (bool fwd : {true, false}) {
                
                auto extracted = extract_extending_graph(graph, node_id, fwd);
                auto expected = ugly_extract_extending_graph(graph, node_id, fwd);
                
                if (!translated_graphs_are_identical(extracted.subgraph, expected.subgraph,
                                                     extracted.back_translation,
                                                     expected.back_translation) ||
                    !node_sets_are_equivalent(extracted.sources, expected.sources,
                                              extracted.back_translation, expected.back_translation) ||
                    !node_sets_are_equivalent(extracted.sinks, expected.sinks,
                                              extracted.back_translation, expected.back_translation)) {
                    cerr << "extending graph extraction failed from node " << node_id << " forward? " << fwd << " on graph:\n";
                    print_graph(graph, cerr);
                    cerr << "fast implementation:\n";
                    dump_connecting_info(extracted);
                    cerr << "slow implementation:\n";
                    dump_connecting_info(expected);
                    exit(1);
                }
            }
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

void test_target_reachability(const BaseGraph& graph, default_random_engine& gen) {
    
    for (size_t rep = 0; rep < 10; ++rep) {
        
        size_t num_targets = uniform_int_distribution<size_t>(0, 1)(gen);
        
        std::vector<uint64_t> targets;
        for (size_t i = 0; i < num_targets; ++i) {
            targets.push_back(uniform_int_distribution<uint64_t>(0, graph.node_size() - 1)(gen));
        }
        
        std::vector<bool> expected(graph.node_size(), false);
        for (size_t i = 0; i < graph.node_size(); ++i) {
            for (auto t : targets) {
                expected[i] = expected[i] || (i == t) || is_reachable(graph, i, t);
            }
        }
        
        auto got = target_reachability(graph, targets);
        
        if (got != expected) {
            cerr << "failed target reachability tests\n";
            exit(1);
        }
    }
}

void test_connected_components(const BaseGraph& graph) {
    
    auto expected = dumb_connected_components(graph);
    auto got = connected_components(graph);
    
    for (auto& comp : expected) {
        std::sort(comp.begin(), comp.end());
    }
    for (auto& comp : got) {
        std::sort(comp.begin(), comp.end());
    }
    std::sort(expected.begin(), expected.end());
    std::sort(got.begin(), got.end());
    
    if (got != expected) {
        cerr << "connected components fail on graph:";
        cerr << cpp_representation(graph, "graph") << '\n';
        cerr << "got:\n";
        for (auto comp : got) {
            for (auto n : comp) {
                cerr << n << ' ';
            }
            cerr << '\n';
        }
        cerr << "expected:\n";
        for (auto comp : expected) {
            for (auto n : comp) {
                cerr << n << ' ';
            }
            cerr << '\n';
        }
        exit(1);
    }
}


void test_shortest_path(const BaseGraph& graph, default_random_engine& gen) {
    
    uniform_int_distribution<uint64_t> node_distr(0, graph.node_size() - 1);
    
    for (size_t rep = 0; rep < 5; ++rep) {
        
        auto from = node_distr(gen);
        auto to = node_distr(gen);
        
        auto got = shortest_path(graph, from, to);
        auto expected = shortest_path_brute_force(graph, from, to);
        
        if (got.size() != expected.size() || !is_valid_path(graph, got)) {
            cerr << "shortest path failed between " << from << " and " << to << " on graph:";
            cerr << cpp_representation(graph, "graph") << '\n';
            cerr << "got:\n";
            for (auto n : got) {
                std::cerr << '\t' << n << '\n';
            }
            cerr << "expected:\n";
            for (auto n : expected) {
                std::cerr << '\t' << n << '\n';
            }
            exit(1);
        }
    }
}

void do_tests(const BaseGraph& graph, const SentinelTableau& tableau, default_random_engine& gen) {
    
    BaseGraph determinized = determinize(graph);
    if (!is_reverse_deterministic(determinized)) {
        cerr << "determinization failed on graph:\n";
        print_graph(graph, cerr);
        cerr << "resulted in determinized graph:\n";
        print_graph(determinized, cerr);
        exit(1);
    }
    assert(is_probably_equivalent(graph, determinized, gen));
    
    auto determinized_tableau = translate_tableau(determinized, tableau);
    rewalk_paths(determinized, determinized_tableau, graph);
    if (!paths_match(graph, determinized)) {
        cerr << "rewalk paths failed on graph:\n";
        print_graph(graph, cerr);
        cerr << "resulted in determinized graph:\n";
        print_graph(determinized, cerr);
        exit(1);
    }
    
    test_count_walks(graph);
    test_count_walks(determinized);
    
    test_topological_order(graph);
    test_topological_order(determinized);
    
    test_antichain_partition(graph);
    test_antichain_partition(determinized);
    
    test_target_reachability(graph, gen);
    test_target_reachability(determinized, gen);
    
    test_shortest_path(graph, gen);
    test_shortest_path(determinized, gen);
    
    test_subgraph_extraction(graph, tableau, gen);
    // the paths aren't necessarily a cover anymore...
    //test_subgraph_extraction(determinized, determinized_tableau, gen);
    
    for (size_t i = 0; i < 4; ++i) {
        uint64_t n = uniform_int_distribution<uint64_t>(0, graph.node_size() - 1)(gen);
        test_minmax_distance(graph, n);
    }
}

int main(int argc, char* argv[]) {
     
    random_device rd;
    default_random_engine gen(rd());

    
    {
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
        auto tableau = add_sentinels(graph1, '^', '$');
        vector<vector<int>> paths{
            {0, 3, 4},
            {0, 1, 2, 3, 4}
        };
        for (size_t i = 0; i < paths.size(); ++i) {
            auto p = graph1.add_path(to_string(i));
            for (auto n : paths[i]) {
                graph1.extend_path(p, n);
            }
        }
        do_tests(graph1, tableau, gen);
    }

    size_t num_reps = 10;
    vector<pair<size_t, size_t>> graph_sizes;
    graph_sizes.emplace_back(5, 3);
    graph_sizes.emplace_back(8, 15);
    graph_sizes.emplace_back(10, 12);
    graph_sizes.emplace_back(20, 35);
    for (auto& sizes : graph_sizes) {
        size_t num_nodes = sizes.first;
        size_t num_edges = sizes.second;
        for (size_t i = 0; i < num_reps; ++i) {
            BaseGraph graph = random_graph(num_nodes, num_edges, true, gen);
            auto tableau = add_sentinels(graph, '^', '$');
            add_random_path_cover(graph, gen, &tableau);
            do_tests(graph, tableau, gen);
        }
    }

    {
        BaseGraph graph1, graph2, expected;
        string seq1 = "ACGCTAC";
        string seq2 = "ACTGTAG";
        string seq_exp = "ACTGCTACG";
        for (auto c : seq1) {
            graph1.add_node(c);
        }
        for (auto c : seq2) {
            graph2.add_node(c);
        }
        for (auto c : seq_exp) {
            expected.add_node(c);
        }
        vector<pair<int, int>> edges1{
            {0, 1},
            {1, 2},
            {2, 3},
            {2, 4},
            {3, 5},
            {4, 5},
            {5, 6}
        };
        vector<pair<int, int>> edges2{
            {0, 1},
            {0, 2},
            {1, 3},
            {2, 3},
            {3, 4},
            {4, 5},
            {5, 6}
        };
        vector<pair<int, int>> edges_exp{
            {0, 1},
            {0, 2},
            {1, 3},
            {2, 3},
            {3, 4},
            {3, 5},
            {4, 6},
            {5, 6},
            {6, 7},
            {6, 8}
        };
        for (auto e : edges1) {
            graph1.add_edge(e.first, e.second);
        }
        for (auto e : edges2) {
            graph2.add_edge(e.first, e.second);
        }
        for (auto e : edges_exp) {
            expected.add_edge(e.first, e.second);
        }

        auto tableau1 = add_sentinels(graph1, '^', '$');
        auto tableau2 = add_sentinels(graph2, '^', '$');
        auto tableau_exp = add_sentinels(expected, '^', '$');

        Alignment alignment{
            {0, 0},
            {1, 1},
            {2, 3},
            {4, 4},
            {5, 5},
            {6, 6}
        };

        test_fuse(graph1, graph2, tableau1, tableau2, alignment, expected);
    }
    
    {
        BaseGraph graph1, graph2, expected;
        string seq1 = "AGTACT";
        string seq2 = "ACGAGT";
        string seq_exp = "ACGTACGT";
        for (auto c : seq1) {
            graph1.add_node(c);
        }
        for (auto c : seq2) {
            graph2.add_node(c);
        }
        for (auto c : seq_exp) {
            expected.add_node(c);
        }
        for (int i = 1; i < seq1.size(); ++i) {
            graph1.add_edge(i - 1, i);
            graph2.add_edge(i - 1, i);
        }
        vector<pair<int, int>> edges_exp{
            {0, 1},
            {0, 2},
            {1, 2},
            {2, 3},
            {2, 4},
            {3, 4},
            {4, 5},
            {4, 6},
            {5, 7},
            {6, 7}
        };
        for (auto e : edges_exp) {
            expected.add_edge(e.first, e.second);
        }
        
        auto tableau1 = add_sentinels(graph1, '^', '$');
        auto tableau2 = add_sentinels(graph2, '^', '$');
        auto tableau_exp = add_sentinels(expected, '^', '$');
        
        Alignment alignment{
            {0, 0},
            {AlignedPair::gap, 1},
            {1, 2},
            {2, AlignedPair::gap},
            {3, 3},
            {4, 4},
            {5, 5}
        };
        
        test_fuse(graph1, graph2, tableau1, tableau2, alignment, expected);
    }
    
    {
        BaseGraph graph, expected;
        string seq = "ACGTACTG";
        string seq_exp = "ACGTTG";
        for (auto c : seq) {
            graph.add_node(c);
        }
        for (auto c : seq_exp) {
            expected.add_node(c);
        }
        vector<pair<int, int>> edges{
            {0, 1},
            {0, 2},
            {1, 3},
            {2, 3},
            {3, 4},
            {4, 5},
            {4, 6},
            {5, 7},
            {6, 7}
        };
        for (auto e : edges) {
            graph.add_edge(e.first, e.second);
        }
        vector<pair<int, int>> edges_exp{
            {0, 1},
            {0, 2},
            {0, 4},
            {1, 3},
            {1, 5},
            {2, 3},
            {3, 0},
            {4, 5}
        };
        for (auto e : edges_exp) {
            expected.add_edge(e.first, e.second);
        }
        
        vector<Alignment> alignments{{
            {0, 4},
            {1, 5}
        }};
        
        test_internal_fuse(graph, alignments, expected);
    }
    
    {
        BaseGraph graph, expected;
        string seq = "ACGTACTGTCTA";
        string seq_exp = "ACGTTGTTA";
        for (auto c : seq) {
            graph.add_node(c);
        }
        for (auto c : seq_exp) {
            expected.add_node(c);
        }
        vector<pair<int, int>> edges{
            {0, 1},
            {0, 2},
            {1, 3},
            {2, 3},
            {3, 4},
            {4, 5},
            {4, 6},
            {5, 7},
            {6, 7},
            {7, 8},
            {8, 9},
            {8, 10},
            {9, 11},
            {10, 11}
        };
        for (auto e : edges) {
            graph.add_edge(e.first, e.second);
        }
        vector<pair<int, int>> edges_exp{
            {0, 1},
            {0, 2},
            {0, 4},
            {1, 3},
            {1, 5},
            {1, 8},
            {2, 3},
            {3, 0},
            {4, 5},
            {5, 6},
            {6, 1},
            {6, 7},
            {7, 8}
        };
        for (auto e : edges_exp) {
            expected.add_edge(e.first, e.second);
        }
        
        vector<Alignment> alignments{
            {
                {0, 4},
                {1, 5}
            },
            {
                {0, 8},
                {1, 9}
            },
        };
        
        test_internal_fuse(graph, alignments, expected);
    }
    
    {
        BaseGraph graph;
        for (char c : string("AGCT")) {
            graph.add_node(c);
        }
        
        vector<pair<int, int>> edges{
            {0, 1},
            {0, 2},
            {1, 3},
            {2, 3}
        };
        for (auto e : edges) {
            graph.add_edge(e.first, e.second);
        }
        
        vector<vector<int>> paths{
            {0, 1, 3}
        };
        
        for (int i = 0; i < paths.size(); ++i) {
            auto p = graph.add_path(to_string(i));
            for (auto n : paths[i]) {
                graph.extend_path(p, n);
            }
        }
        
        auto tableau = add_sentinels(graph, '^', '$');
        
        purge_uncovered_nodes(graph, tableau);
        
        assert(graph.node_size() == 3 + 2);
        assert(graph.label(tableau.src_id) == tableau.src_sentinel);
        assert(graph.label(tableau.snk_id) == tableau.snk_sentinel);
    }
    
    {
        BaseGraph graph;
        for (auto c : string("AAACCAAG")) {
            graph.add_node(c);
        }
        graph.add_edge(0, 1);
        graph.add_edge(0, 2);
        graph.add_edge(0, 3);
        graph.add_edge(1, 4);
        graph.add_edge(2, 4);
        graph.add_edge(3, 4);
        graph.add_edge(4, 5);
        graph.add_edge(4, 6);
        graph.add_edge(4, 7);
        graph.add_edge(5, 7);
        graph.add_edge(6, 7);
        std::vector<std::vector<int>> paths{
            {0, 1, 4, 5, 7},
            {0, 2, 4, 6, 7},
            {0, 3, 4, 7}
        };
        for (size_t i = 0; i < paths.size(); ++i) {
            auto p = graph.add_path(to_string(i));
            auto path = paths[i];
            for (auto n : path) {
                graph.extend_path(p, n);
            }
        }
        
        BaseGraph expected;
        for (auto c : string("AACCAAG")) {
            expected.add_node(c);
        }
        expected.add_edge(0, 1);
        expected.add_edge(0, 2);
        expected.add_edge(1, 3);
        expected.add_edge(2, 3);
        expected.add_edge(3, 4);
        expected.add_edge(3, 5);
        expected.add_edge(3, 6);
        expected.add_edge(4, 6);
        expected.add_edge(5, 6);
        std::vector<std::vector<int>> exp_paths{
            {0, 1, 3, 4, 6},
            {0, 1, 3, 5, 6},
            {0, 2, 3, 6}
        };
        for (size_t i = 0; i < exp_paths.size(); ++i) {
            auto p = expected.add_path(to_string(i));
            auto path = exp_paths[i];
            for (auto n : path) {
                expected.extend_path(p, n);
            }
        }
        
        auto tableau = add_sentinels(graph, '^', '$');
        auto dummy = add_sentinels(expected, '^', '$');
        
        simplify_bubbles(graph, tableau);
        
        assert(possibly_isomorphic(graph, expected));
    }
    
    cerr << "passed all tests!" << endl;
}
