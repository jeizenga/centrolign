#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <unordered_set>
#include <vector>
#include <list>
#include <string>
#include <random>
#include <iostream>

#include "centrolign/graph.hpp"
#include "centrolign/modify_graph.hpp"
#include "centrolign/simplifier.hpp"
#include "centrolign/test_util.hpp"

using namespace std;
using namespace centrolign;

void check_paths(const BaseGraph& graph, const BaseGraph& expanded_graph) {
    for (uint64_t path_id = 0; path_id < graph.path_size(); ++path_id) {
        
        auto simp_path_id = expanded_graph.path_id(graph.path_name(path_id));
        
        auto seq1 = path_to_string(graph, graph.path(path_id));
        auto seq2 = path_to_string(expanded_graph, expanded_graph.path(simp_path_id));
        
        if (seq1 != seq2) {
            cerr << "did not find matching path sequence for path " << graph.path_name(path_id) << " in simplified graph.\n";
            cerr << "original graph:\n";
            print_graph(graph, cerr);
            cerr << "simplified graph:\n";
            print_graph(expanded_graph, cerr);
            exit(1);
        }
        
    }
}


int main(int argc, char* argv[]) {
     
    random_device rd;
    default_random_engine gen(rd());
    
    {
        BaseGraph graph;
        for (char c : string("AAAAAG")) {
            graph.add_node(c);
        }
        vector<pair<int, int>> edges{
            {0, 1},
            {1, 2},
            {1, 5},
            {1, 3},
            {2, 3},
            {2, 4},
            {3, 4},
            {5, 3}
        };
        for (auto e : edges) {
            graph.add_edge(e.first, e.second);
        }
        vector<string> pnames{"1", "2"};
        vector<vector<int>> paths{
            {0, 1, 2, 3, 4},
            {0, 1, 5, 3, 4}
        };
        for (int i = 0; i < paths.size(); ++i) {
            auto pid = graph.add_path(pnames[i]);
            for (auto n : paths[i]) {
                graph.extend_path(pid, n);
            }
        }

        auto tableau = add_sentinels(graph, '^', '$');

        vector<uint64_t> nodes{3, 4};
        size_t dist = 3;

        Simplifier simplifier;

        auto expanded = simplifier.targeted_simplify(graph, tableau, nodes, dist);

        check_paths(graph, expanded.graph);
    }

    {
        Simplifier simplifier;

        uniform_int_distribution<size_t> size_distr(5, 8);
        uniform_int_distribution<size_t> simplify_dist_distr(0, 10);
        uniform_int_distribution<size_t> num_simplify_nodes_distr(1, 5);

        size_t num_reps = 100;
        for (size_t rep = 0; rep < num_reps; ++rep) {

            size_t num_nodes = size_distr(gen);

            BaseGraph graph = random_challenge_graph(num_nodes, gen);
            add_random_path_cover(graph, gen);

            uniform_int_distribution<uint64_t> node_distr(0, graph.node_size() - 1);

            auto tableau = add_sentinels(graph, '^', '$');

            size_t num_simplify_nodes = num_simplify_nodes_distr(gen);
            vector<uint64_t> nodes;
            for (size_t i = 0; i < num_simplify_nodes; ++i) {
                nodes.push_back(node_distr(gen));
            }
            sort(nodes.begin(), nodes.end());
            auto e = unique(nodes.begin(), nodes.end());
            nodes.resize(e - nodes.begin());

            size_t dist = simplify_dist_distr(gen);


//            cerr << "next graph:\n";
//            print_graph(graph, cerr);
//            cerr << "simplify dist " << dist << " from nodes:\n";
//            for (auto n : nodes) {
//                cerr << '\t' << n << '\n';
//            }

            auto expanded = simplifier.targeted_simplify(graph, tableau, nodes, dist);

            check_paths(graph, expanded.graph);

        }
    }

    {
        BaseGraph graph;
        for (char c : string("ACAACGTACGTA")) {
            graph.add_node(c);
        }
        vector<pair<int, int>> edges{
            {0, 1},
            {0, 2},
            {1, 3},
            {2, 3},
            {3, 4},
            {3, 5},
            {3, 6},
            {4, 7},
            {5, 7},
            {6, 7},
            {7, 8},
            {7, 9},
            {7, 10},
            {8, 11},
            {9, 11},
            {10, 11}
        };
        for (auto e : edges) {
            graph.add_edge(e.first, e.second);
        }
        vector<string> pnames{"1", "2", "3", "4"};
        vector<vector<int>> paths{
            {0, 1, 3, 4, 7, 8, 11},
            {0, 2, 3, 5, 7, 9, 11},
            {0, 1, 3, 6, 7, 10, 11},
            {0, 2, 3, 6, 7, 8, 11}
        };
        for (int i = 0; i < paths.size(); ++i) {
            auto pid = graph.add_path(pnames[i]);
            for (auto n : paths[i]) {
                graph.extend_path(pid, n);
            }
        }

        auto tableau = add_sentinels(graph, '^', '$');

        Simplifier simp;
        simp.min_dist_window = 5;
        simp.max_walks = 8;

        auto expanded_graph = simp.simplify(graph, tableau);

        BaseGraph expected;
        for (auto c : string("^ACAACGTAAACGTA$")) {
            expected.add_node(c);
        }

        vector<uint64_t> back_trans{
            tableau.src_id,
            0, 1, 2, 3, 4, 5, 6, 7, 7, 7, 8, 9, 10, 11,
            tableau.snk_id
        };
        vector<pair<int, int>> expected_edges{
            {0, 1},
            {1, 2},
            {1, 3},
            {2, 4},
            {3, 4},
            {4, 5},
            {4, 6},
            {4, 7},
            {5, 8},
            {6, 9},
            {7, 8},
            {7, 10},
            {8, 11},
            {9, 12},
            {10, 13},
            {11, 14},
            {12, 14},
            {13, 14},
            {14, 15}
        };
        for (auto e : expected_edges) {
            expected.add_edge(e.first, e.second);
        }
        vector<string> expected_pnames{"1", "2", "3", "4"};
        vector<vector<int>> expected_paths{
            {1, 2, 4, 5, 8, 11, 14},
            {1, 3, 4, 6, 9, 12, 14},
            {1, 2, 4, 7, 10, 13, 14},
            {1, 3, 4, 7, 8, 11, 14}
        };

        for (int i = 0; i < expected_paths.size(); ++i) {
            auto pid = expected.add_path(expected_pnames[i]);
            for (auto n : expected_paths[i]) {
                expected.extend_path(pid, n);
            }
        }

        if (!possibly_isomorphic(expanded_graph.graph, expected) ||
            !translations_possibly_consistent(expanded_graph.graph, expected,
                                              expanded_graph.back_translation, back_trans)) {

            cerr << "failed graph simplification test 1\n";
            cerr << "got:\n";
            print_graph(expanded_graph.graph, cerr);
            cerr << "expected:\n";
            print_graph(expected, cerr);

            exit(1);
        }
    }

    // respects min distance window
    {
        BaseGraph graph;
        for (char c : string("ACAAAA")) {
            graph.add_node(c);
        }
        vector<pair<int, int>> edges{
            {0, 1},
            {0, 2},
            {1, 3},
            {2, 3},
            {3, 4},
            {3, 5},
            {4, 5}
        };
        for (auto e : edges) {
            graph.add_edge(e.first, e.second);
        }
        vector<string> pnames{"1", "2"};
        vector<vector<int>> paths{
            {0, 1, 3, 4, 5},
            {0, 2, 3, 5}
        };
        for (int i = 0; i < paths.size(); ++i) {
            auto pid = graph.add_path(pnames[i]);
            for (auto n : paths[i]) {
                graph.extend_path(pid, n);
            }
        }

        auto tableau = add_sentinels(graph, '^', '$');

        Simplifier simp;
        simp.min_dist_window = 4;
        simp.max_walks = 3;

        auto expanded_graph = simp.simplify(graph, tableau);

        BaseGraph expected;
        for (auto c : string("^ACAAAAA$")) {
            expected.add_node(c);
        }
        vector<uint64_t> back_trans{tableau.src_id, 0, 1, 3, 4, 2, 3, 5, tableau.snk_id
        };
        vector<pair<int, int>> expected_edges{
            {0, 1},
            {1, 2},
            {2, 3},
            {3, 4},
            {4, 7},
            {1, 5},
            {5, 6},
            {6, 7},
            {7, 8}
        };
        for (auto e : expected_edges) {
            expected.add_edge(e.first, e.second);
        }
        vector<vector<int>> expected_paths{
            {1, 2, 3, 4, 7},
            {1, 5, 6, 7}
        };

        for (int i = 0; i < expected_paths.size(); ++i) {
            auto pid = expected.add_path(pnames[i]);
            for (auto n : expected_paths[i]) {
                expected.extend_path(pid, n);
            }
        }

        if (!possibly_isomorphic(expanded_graph.graph, expected) ||
            !translations_possibly_consistent(expanded_graph.graph, expected,
                                              expanded_graph.back_translation, back_trans)) {

            cerr << "failed graph simplification test 2\n";
            cerr << "got:\n";
            print_graph(expanded_graph.graph, cerr);
            cerr << "expected:\n";
            print_graph(expected, cerr);

            exit(1);
        }
    }

    // can count hierarchical walks
    {
        BaseGraph graph;
        for (char c : string("ACCAAACAAA")) {
            graph.add_node(c);
        }
        vector<pair<int, int>> edges{
            {0, 1},
            {0, 5},
            {1, 2},
            {1, 3},
            {2, 4},
            {3, 4},
            {4, 9},
            {5, 6},
            {5, 7},
            {6, 8},
            {7, 8},
            {8, 9}
        };
        for (auto e : edges) {
            graph.add_edge(e.first, e.second);
        }
        vector<string> pnames{"1", "2", "3", "4"};
        vector<vector<int>> paths{
            {0, 1, 2, 4, 9},
            {0, 1, 3, 4, 9},
            {0, 5, 6, 8, 9},
            {0, 5, 7, 8, 9}
        };
        for (int i = 0; i < paths.size(); ++i) {
            auto pid = graph.add_path(pnames[i]);
            for (auto n : paths[i]) {
                graph.extend_path(pid, n);
            }
        }

        auto tableau = add_sentinels(graph, '^', '$');

        Simplifier simp;
        simp.min_dist_window = 5;
        simp.max_walks = 3;

        auto expanded_graph = simp.simplify(graph, tableau);

        vector<uint64_t> back_trans{
            0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
            tableau.src_id, tableau.snk_id
        };

        // this one actually gets remerged back into the same graph
        if (!possibly_isomorphic(expanded_graph.graph, graph) ||
            !translations_possibly_consistent(expanded_graph.graph, graph,
                                              expanded_graph.back_translation, back_trans)) {

            cerr << "failed graph simplification test 3\n";
            cerr << "got:\n";
            print_graph(expanded_graph.graph, cerr);
            cerr << "expected:\n";
            print_graph(graph, cerr);

            exit(1);
        }
    }

    // preserves long alleles
    {
        BaseGraph graph;
        for (char c : string("AACAAA")) {
            graph.add_node(c);
        }
        vector<pair<int, int>> edges{
            {0, 1},
            {0, 2},
            {1, 5},
            {2, 3},
            {3, 4},
            {4, 5}
        };
        for (auto e : edges) {
            graph.add_edge(e.first, e.second);
        }
        vector<string> pnames{"1", "2"};
        vector<vector<int>> paths{
            {0, 1, 5},
            {0, 2, 3, 4, 5}
        };
        for (int i = 0; i < paths.size(); ++i) {
            auto pid = graph.add_path(pnames[i]);
            for (auto n : paths[i]) {
                graph.extend_path(pid, n);
            }
        }

        auto tableau = add_sentinels(graph, '^', '$');

        Simplifier simp;
        simp.min_dist_window = 5;
        simp.preserve_bubble_size = 4;
        simp.max_walks = 1;

        auto expanded_graph = simp.simplify(graph, tableau);

        BaseGraph expected;
        for (auto c : string("^AACAAA$")) {
            expected.add_node(c);
        }
        vector<uint64_t> back_trans{
            tableau.src_id,
            0, 1, 2, 3, 4, 5,
            tableau.snk_id
        };
        vector<pair<int, int>> expected_edges{
            {0, 1},
            {1, 2},
            {1, 3},
            {2, 6},
            {3, 4},
            {4, 5},
            {5, 6},
            {6, 7}
        };
        for (auto e : expected_edges) {
            expected.add_edge(e.first, e.second);
        }
        vector<vector<int>> expected_paths{
            {1, 2, 6},
            {1, 3, 4, 5, 6}
        };

        for (int i = 0; i < expected_paths.size(); ++i) {
            auto pid = expected.add_path(pnames[i]);
            for (auto n : expected_paths[i]) {
                expected.extend_path(pid, n);
            }
        }

        if (!possibly_isomorphic(expanded_graph.graph, expected) ||
            !translations_possibly_consistent(expanded_graph.graph, expected,
                                              expanded_graph.back_translation, back_trans)) {

            cerr << "failed graph simplification test 4\n";
            cerr << "got:\n";
            print_graph(expanded_graph.graph, cerr);
            cerr << "expected:\n";
            print_graph(expected, cerr);

            exit(1);
        }
    }

    {
        BaseGraph graph;
        for (char c : string("ACAAACACA")) {
            graph.add_node(c);
        }
        vector<pair<int, int>> edges{
            {0, 1},
            {0, 2},
            {1, 3},
            {2, 3},
            {3, 4},
            {4, 5},
            {4, 6},
            {5, 6},
            {6, 7},
            {6, 8},
            {7, 8}
        };
        for (auto e : edges) {
            graph.add_edge(e.first, e.second);
        }
        vector<string> pnames{"1", "2"};
        vector<vector<int>> paths{
            {0, 1, 3, 4, 5, 6, 7, 8},
            {0, 2, 3, 4,    6,    8}
        };
        for (int i = 0; i < paths.size(); ++i) {
            auto pid = graph.add_path(pnames[i]);
            for (auto n : paths[i]) {
                graph.extend_path(pid, n);
            }
        }

        auto tableau = add_sentinels(graph, '^', '$');

        Simplifier simp;

        BaseGraph expected1;
        vector<uint64_t> back_trans1;

        {
            vector<uint64_t> simp_ids{1, 3, 5};
            size_t simp_dist = 0;

            auto expanded_graph = simp.targeted_simplify(graph, tableau, simp_ids, simp_dist);

            BaseGraph& expected = expected1;

            for (auto c : string("^ACAAAAACACA$")) {
                expected.add_node(c);
            }
            back_trans1 = vector<uint64_t>{
                tableau.src_id,
                0, 1, 2, 3, 3, 4, 4, 5, 6, 7, 8,
                tableau.snk_id
            };
            auto& back_trans = back_trans1;
            vector<pair<int, int>> expected_edges{
                {0, 1},
                {1, 2},
                {1, 3},
                {2, 4},
                {3, 5},
                {4, 6},
                {5, 7},
                {6, 8},
                {7, 9},
                {8, 9},
                {9, 10},
                {9, 11},
                {10, 11},
                {11, 12}
            };
            for (auto e : expected_edges) {
                expected.add_edge(e.first, e.second);
            }
            vector<vector<int>> expected_paths{
                {1, 2, 4, 6, 8, 9, 10, 11},
                {1, 3, 5, 7,    9,     11}
            };

            for (int i = 0; i < expected_paths.size(); ++i) {
                auto pid = expected.add_path(pnames[i]);
                for (auto n : expected_paths[i]) {
                    expected.extend_path(pid, n);
                }
            }

            if (!possibly_isomorphic(expanded_graph.graph, expected) ||
                !translations_possibly_consistent(expanded_graph.graph, expected,
                                                  expanded_graph.back_translation, back_trans)) {

                cerr << "failed graph simplification test\n";
                cerr << "got:\n";
                print_graph(expanded_graph.graph, cerr);
                cerr << "expected:\n";
                print_graph(expected, cerr);

                exit(1);
            }
        }

        {
            vector<uint64_t> simp_ids{2};
            size_t simp_dist = 2;

            auto expanded_graph = simp.targeted_simplify(graph, tableau, simp_ids, simp_dist);

            if (!possibly_isomorphic(expanded_graph.graph, expected1) ||
                !translations_possibly_consistent(expanded_graph.graph, expected1,
                                                  expanded_graph.back_translation, back_trans1)) {

                cerr << "failed graph simplification test\n";
                cerr << "got:\n";
                print_graph(expanded_graph.graph, cerr);
                cerr << "expected:\n";
                print_graph(expected1, cerr);

                exit(1);
            }
        }

        BaseGraph expected2;
        vector<uint64_t> back_trans2;

        {
            vector<uint64_t> simp_ids{4, 5, 8};
            size_t simp_dist = 2;

            auto expanded_graph = simp.targeted_simplify(graph, tableau, simp_ids, simp_dist);

            BaseGraph& expected = expected2;

            for (auto c : string("^ACAAACAACA$")) {
                expected.add_node(c);
            }
            back_trans2 = vector<uint64_t>{
                tableau.src_id,
                0, 1, 2, 3, 4, 5, 6, 6, 7, 8,
                tableau.snk_id
            };
            auto& back_trans = back_trans2;
            vector<pair<int, int>> expected_edges{
                {0, 1},
                {1, 2},
                {1, 3},
                {2, 4},
                {3, 4},
                {4, 5},
                {5, 6},
                {5, 8},
                {6, 7},
                {7, 9},
                {8, 10},
                {9, 10},
                {10, 11}
            };
            for (auto e : expected_edges) {
                expected.add_edge(e.first, e.second);
            }
            vector<vector<int>> expected_paths{
                {1, 2, 4, 5, 6, 7, 9, 10},
                {1, 3, 4, 5,    8,    10}
            };

            for (int i = 0; i < expected_paths.size(); ++i) {
                auto pid = expected.add_path(pnames[i]);
                for (auto n : expected_paths[i]) {
                    expected.extend_path(pid, n);
                }
            }

            if (!possibly_isomorphic(expanded_graph.graph, expected) ||
                !translations_possibly_consistent(expanded_graph.graph, expected,
                                                  expanded_graph.back_translation, back_trans)) {

                cerr << "failed graph simplification test\n";
                cerr << "got:\n";
                print_graph(expanded_graph.graph, cerr);
                cerr << "expected:\n";
                print_graph(expected, cerr);

                exit(1);
            }
        }

        {
            vector<uint64_t> simp_ids{5};
            size_t simp_dist = 1;

            auto expanded_graph = simp.targeted_simplify(graph, tableau, simp_ids, simp_dist);

            if (!possibly_isomorphic(expanded_graph.graph, expected2) ||
                !translations_possibly_consistent(expanded_graph.graph, expected2,
                                                  expanded_graph.back_translation, back_trans2)) {

                cerr << "failed graph simplification test\n";
                cerr << "got:\n";
                print_graph(expanded_graph.graph, cerr);
                cerr << "expected:\n";
                print_graph(expected1, cerr);

                exit(1);
            }
        }
    }
    

    // nested bubbles
    {
        BaseGraph graph;
        for (char c : string("AACAACAACACCACCA")) {
            graph.add_node(c);
        }
        vector<pair<int, int>> edges{
            {0, 1},
            {0, 8},
            {1, 2},
            {1, 3},
            {2, 4},
            {3, 4},
            {4, 5},
            {4, 6},
            {5, 7},
            {6, 7},
            {7, 15},
            {8, 9},
            {8, 10},
            {9, 11},
            {10, 11},
            {11, 12},
            {11, 13},
            {12, 14},
            {13, 14},
            {14, 15}
        };
        for (auto e : edges) {
            graph.add_edge(e.first, e.second);
        }
        vector<string> pnames{"1", "2", "3", "4"};
        vector<vector<int>> paths{
            {0, 1, 2, 4, 5, 7, 15},
            {0, 1, 3, 4, 6, 7, 15},
            {0, 8, 9, 11, 12, 14, 15},
            {0, 8, 10, 11, 13, 14, 15}
        };
        for (int i = 0; i < paths.size(); ++i) {
            auto pid = graph.add_path(pnames[i]);
            for (auto n : paths[i]) {
                graph.extend_path(pid, n);
            }
        }
        
        auto tableau = add_sentinels(graph, '^', '$');
        
        BaseGraph expected1, expected2;
        vector<uint64_t> back_trans1, back_trans2;
        
        {
            auto& expected = expected1;
            auto& back_trans = back_trans1;
            
            for (auto c : string("^AACAAACAACACCACCA$")) {
                expected.add_node(c);
            }
            back_trans = vector<uint64_t>{
                tableau.src_id,
                0, 1, 2, 3, 4, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                tableau.snk_id
            };
            vector<pair<int, int>> expected_edges{
                {0, 1},
                {1, 2},
                {1, 10},
                {2, 3},
                {2, 4},
                {3, 5},
                {4, 6},
                {5, 7},
                {6, 8},
                {7, 9},
                {8, 9},
                {9, 17},
                {10, 11},
                {10, 12},
                {11, 13},
                {12, 13},
                {13, 14},
                {13, 15},
                {14, 16},
                {15, 16},
                {16, 17},
                {17, 18}
            };
            for (auto e : expected_edges) {
                expected.add_edge(e.first, e.second);
            }
            vector<vector<int>> expected_paths{
                {1, 2, 3, 5, 7, 9, 17},
                {1, 2, 4, 6, 8, 9, 17},
                {1, 10, 11, 13, 14, 16, 17},
                {1, 10, 12, 13, 15, 16, 17}
            };
            
            for (int i = 0; i < expected_paths.size(); ++i) {
                auto pid = expected.add_path(pnames[i]);
                for (auto n : expected_paths[i]) {
                    expected.extend_path(pid, n);
                }
            }
        }
        
        {
            auto& expected = expected2;
            auto& back_trans = back_trans2;
            
            for (auto c : string("^AACAAACAACACCCACCA$")) {
                expected.add_node(c);
            }
            back_trans = vector<uint64_t>{
                tableau.src_id,
                0, 1, 2, 3, 4, 4, 5, 6, 7, 8, 9, 10, 11, 11, 12, 13, 14, 15,
                tableau.snk_id
            };
            vector<pair<int, int>> expected_edges{
                {0, 1},
                {1, 2},
                {1, 10},
                {2, 3},
                {2, 4},
                {3, 5},
                {4, 6},
                {5, 7},
                {6, 8},
                {7, 9},
                {8, 9},
                {9, 18},
                {10, 11},
                {10, 12},
                {11, 13},
                {12, 14},
                {13, 15},
                {14, 16},
                {15, 17},
                {16, 17},
                {17, 18},
                {18, 19}
            };
            for (auto e : expected_edges) {
                expected.add_edge(e.first, e.second);
            }
            vector<vector<int>> expected_paths{
                {1, 2, 3, 5, 7, 9, 18},
                {1, 2, 4, 6, 8, 9, 18},
                {1, 10, 11, 13, 15, 17, 18},
                {1, 10, 12, 14, 16, 17, 18}
            };
            
            for (int i = 0; i < expected_paths.size(); ++i) {
                auto pid = expected.add_path(pnames[i]);
                for (auto n : expected_paths[i]) {
                    expected.extend_path(pid, n);
                }
            }
        }
        
        Simplifier simp;
        
        {

            vector<uint64_t> simp_ids{1};
            size_t simp_dist = 2;
            
            auto expanded_graph = simp.targeted_simplify(graph, tableau, simp_ids, simp_dist);
                        
            if (!possibly_isomorphic(expanded_graph.graph, expected1) ||
                !translations_possibly_consistent(expanded_graph.graph, expected1,
                                                  expanded_graph.back_translation, back_trans1)) {
                
                cerr << "failed graph simplification test\n";
                cerr << "got:\n";
                print_graph(expanded_graph.graph, cerr);
                cerr << "expected:\n";
                print_graph(graph, cerr);
                
                exit(1);
            }
        }
        
        {
            
            vector<uint64_t> simp_ids{1, 4, 9, 11};
            size_t simp_dist = 1;
            
            auto expanded_graph = simp.targeted_simplify(graph, tableau, simp_ids, simp_dist);
            
            if (!possibly_isomorphic(expanded_graph.graph, expected2) ||
                !translations_possibly_consistent(expanded_graph.graph, expected2,
                                                  expanded_graph.back_translation, back_trans2)) {
                
                cerr << "failed graph simplification test\n";
                cerr << "got:\n";
                print_graph(expanded_graph.graph, cerr);
                cerr << "expected:\n";
                print_graph(graph, cerr);
                
                exit(1);
            }
        }
        
        {
            
            vector<uint64_t> simp_ids{0};
            size_t simp_dist = 1;
            
            auto expanded_graph = simp.targeted_simplify(graph, tableau, simp_ids, simp_dist);
            
            if (!possibly_isomorphic(expanded_graph.graph, expected2) ||
                !translations_possibly_consistent(expanded_graph.graph, expected2,
                                                  expanded_graph.back_translation, back_trans2)) {
                
                cerr << "failed graph simplification test\n";
                cerr << "got:\n";
                print_graph(expanded_graph.graph, cerr);
                cerr << "expected:\n";
                print_graph(graph, cerr);
                
                exit(1);
            }
        }
    }
        
    cerr << "passed all tests!" << endl;
}
