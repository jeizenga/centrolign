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


int main(int argc, char* argv[]) {
     
    random_device rd;
    default_random_engine gen(rd());
    
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
        
//        BaseGraph expected;
//        for (auto c : string("^ACCCAAACAAAA$")) {
//            expected.add_node(c);
//        }
        vector<uint64_t> back_trans{
            0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
            tableau.src_id, tableau.snk_id
        };
//        vector<pair<int, int>> expected_edges{
//            {0, 1},
//            {1, 2},
//            {1, 4},
//            {1, 7},
//            {1, 9},
//            {2, 3},
//            {3, 6},
//            {4, 5},
//            {5, 6},
//            {6, 12},
//            {7, 8},
//            {8, 11},
//            {9, 10},
//            {10, 11},
//            {11, 12},
//            {12, 13}
//        };
//        for (auto e : expected_edges) {
//            expected.add_edge(e.first, e.second);
//        }
//        vector<vector<int>> expected_paths{
//            {1, 2, 3, 6, 12},
//            {1, 4, 5, 6, 12},
//            {1, 7, 8, 11, 12},
//            {1, 9, 10, 11, 12}
//        };
//
//        for (int i = 0; i < expected_paths.size(); ++i) {
//            auto pid = expected.add_path(pnames[i]);
//            for (auto n : expected_paths[i]) {
//                expected.extend_path(pid, n);
//            }
//        }
        
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
    
    cerr << "passed all tests!" << endl;
}
