#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <utility>
#include <string>

#include "centrolign/stitcher.hpp"
#include "centrolign/graph.hpp"
#include "centrolign/modify_graph.hpp"
#include "centrolign/chain_merge.hpp"
#include "centrolign/anchorer.hpp"
#include "centrolign/match_finder.hpp"
#include "centrolign/test_util.hpp"

using namespace std;
using namespace centrolign;

int main(int argc, char* argv[]) {
    
    Stitcher stitcher;
    
    // test the partitioning algorithm
    {
        BaseGraph graph1;
        for (auto c : std::string("GATCGAT")) {
            graph1.add_node(c);
        }
        
        std::vector<std::pair<int, int>> graph1_edges{
            {0, 1},
            {1, 2},
            {2, 3},
            {3, 4},
            {4, 5},
            {5, 6}
        };
        
        std::vector<std::vector<int>> graph1_paths{
            {0, 1, 2, 3, 4, 5, 6}
        };
        
        for (auto e : graph1_edges) {
            graph1.add_edge(e.first, e.second);
        }
        
        for (size_t i = 0; i < graph1_paths.size(); ++i) {
            auto p = graph1.add_path(std::to_string(i));
            for (auto n : graph1_paths[i]) {
                graph1.extend_path(p, n);
            }
        }
        
        BaseGraph graph2;
        for (auto c : std::string("TAGCTAG")) {
            graph2.add_node(c);
        }
        
        std::vector<std::pair<int, int>> graph2_edges{
            {0, 1},
            {1, 2},
            {2, 3},
            {3, 4},
            {4, 5},
            {5, 6}
        };
        
        std::vector<std::vector<int>> graph2_paths{
            {0, 1, 2, 3, 4, 5, 6}
        };
        
        for (auto e : graph2_edges) {
            graph2.add_edge(e.first, e.second);
        }
        
        for (size_t i = 0; i < graph2_paths.size(); ++i) {
            auto p = graph2.add_path(std::to_string(i));
            for (auto n : graph2_paths[i]) {
                graph2.extend_path(p, n);
            }
        }
                
        BaseGraph sub1;
        for (auto c : std::string("GAT")) {
            sub1.add_node(c);
        }
        
        std::vector<std::pair<int, int>> sub1_edges{
            {0, 1},
            {1, 2}
        };
        
        std::vector<std::vector<int>> sub1_paths{
            {0, 1, 2}
        };
        
        for (auto e : sub1_edges) {
            sub1.add_edge(e.first, e.second);
        }
        
        for (size_t i = 0; i < sub1_paths.size(); ++i) {
            auto p = sub1.add_path(std::to_string(i));
            for (auto n : sub1_paths[i]) {
                sub1.extend_path(p, n);
            }
        }
        
        BaseGraph sub2;
        for (auto c : std::string("TAG")) {
            sub2.add_node(c);
        }
        
        std::vector<std::pair<int, int>> sub2_edges{
            {0, 1},
            {1, 2}
        };
        
        std::vector<std::vector<int>> sub2_paths{
            {0, 1, 2}
        };
        
        for (auto e : sub2_edges) {
            sub2.add_edge(e.first, e.second);
        }
        
        for (size_t i = 0; i < sub2_paths.size(); ++i) {
            auto p = sub2.add_path(std::to_string(i));
            for (auto n : sub2_paths[i]) {
                sub2.extend_path(p, n);
            }
        }
        
        auto tableau1 = add_sentinels(graph1, '^', '$');
        auto tableau2 = add_sentinels(graph2, '^', '$');
        
        std::vector<match_set_t> match_sets(4);
        // A
        match_sets[0].walks1.push_back({1});
        match_sets[0].walks1.push_back({5});
        match_sets[0].walks2.push_back({1});
        match_sets[0].walks2.push_back({5});
        // C
        match_sets[1].walks1.push_back({3});
        match_sets[1].walks2.push_back({3});
        // G
        match_sets[2].walks1.push_back({0});
        match_sets[2].walks1.push_back({4});
        match_sets[2].walks2.push_back({2});
        match_sets[2].walks2.push_back({6});
        // T
        match_sets[3].walks1.push_back({2});
        match_sets[3].walks1.push_back({6});
        match_sets[3].walks2.push_back({0});
        match_sets[3].walks2.push_back({4});
        
        for (auto& m : match_sets) {
            m.count1 = m.walks1.size();
            m.count2 = m.walks2.size();
        }
        
        std::vector<anchor_t> chain(1);
        chain.front().walk1 = match_sets[1].walks1.front();
        chain.front().walk2 = match_sets[1].walks2.front();
        chain.front().count1 = 1;
        chain.front().count2 = 1;
        
        ChainMerge merge1(graph1, tableau1);
        ChainMerge merge2(graph2, tableau2);
        
        auto stitch_graphs = stitcher.extract_stitch_graphs(chain, graph1, graph2,
                                                            tableau1, tableau2,
                                                            merge1, merge2);
        
        stitcher.project_paths(graph1, graph2, stitch_graphs);
        
        auto divvied = stitcher.divvy_matches(match_sets, graph1, graph2, stitch_graphs);
        
        assert(stitch_graphs.size() == 2);
        assert(divvied.size() == 2);
        
        for (size_t i = 0; i < stitch_graphs.size(); ++i) {
            
            auto& stitch_graph = stitch_graphs[i];
            
            assert(stitch_graph.first.back_translation.size() == stitch_graph.first.subgraph.node_size());
            assert(stitch_graph.second.back_translation.size() == stitch_graph.second.subgraph.node_size());
            
            if (!possibly_isomorphic(sub1, stitch_graph.first.subgraph)) {
                cerr << "isomorphism 1 failure\n";
                exit(1);
            }
            if (!possibly_isomorphic(sub2, stitch_graph.second.subgraph)) {
                cerr << "isomorphism 2 failure\n";
                exit(1);
            }
            
            assert(divvied[i].size() == 3);
            for (auto& set : divvied[i]) {
                assert(set.walks1.size() == 1);
                assert(set.walks2.size() == 1);
            }
        }
    }
    
    BaseGraph graph1, graph2;
    for (auto c : string("ACCAGTCGTTGA")) {
        graph1.add_node(c);
    }
    for (auto c : string("GATCGTGAACTATGC")) {
        graph2.add_node(c);
    }
    vector<pair<int, int>> edges1{
        {0, 1},
        {0, 2},
        {1, 3},
        {2, 3},
        {3, 4},
        {4, 5},
        {5, 6},
        {5, 7},
        {6, 8},
        {7, 8},
        {8, 9},
        {8, 10},
        {9, 10},
        {10, 11}
    };
    vector<pair<int, int>> edges2{
        {0, 1},
        {1, 2},
        {1, 3},
        {2, 3},
        {3, 4},
        {4, 5},
        {4, 6},
        {5, 7},
        {6, 7},
        {7, 8},
        {7, 9},
        {8, 10},
        {9, 10},
        {10, 11},
        {10, 12},
        {11, 13},
        {12, 13},
        {13, 14}
    };
    for (auto e : edges1) {
        graph1.add_edge(e.first, e.second);
    }
    for (auto e : edges2) {
        graph2.add_edge(e.first, e.second);
    }
    vector<vector<int>> paths1{
        {0, 1, 3, 4, 5, 6, 8,    10, 11},
        {0, 2, 3, 4, 5, 7, 8, 9, 10, 11}
    };
    vector<vector<int>> paths2{
        {0, 1,    3, 4, 5, 7, 8, 10, 11, 13, 14},
        {0, 1, 2, 3, 4, 6, 7, 9, 10, 12, 13, 14}
    };
    int n = 0;
    for (auto p : paths1) {
        auto i = graph1.add_path(to_string(n++));
        for (auto j : p) {
            graph1.extend_path(i, j);
        }
    }
    for (auto p : paths2) {
        auto i = graph2.add_path(to_string(n++));
        for (auto j : p) {
            graph2.extend_path(i, j);
        }
    }
    
    auto tableau1 = add_sentinels(graph1, '^', '$');
    auto tableau2 = add_sentinels(graph2, '^', '$');
    
    ChainMerge chain_merge1(graph1, tableau1);
    ChainMerge chain_merge2(graph2, tableau2);
    
    vector<vector<pair<int, int>>> anchors{
        {{0, 1}, {1, 3}},
        {{4, 4}, {5, 5}},
        {{9, 12}, {10, 13}}
    };
    
    vector<anchor_t> anchor_chain;
    for (auto a : anchors) {
        anchor_chain.emplace_back();
        anchor_chain.back().count1 = 1;
        anchor_chain.back().count2 = 1;
        for (auto b : a) {
            anchor_chain.back().walk1.push_back(b.first);
            anchor_chain.back().walk2.push_back(b.second);
        }
    }
    
    Alignment expected;
    expected.emplace_back(AlignedPair::gap, 0);
    expected.emplace_back(0, 1);
    expected.emplace_back(1, 3);
    expected.emplace_back(3, AlignedPair::gap);
    expected.emplace_back(4, 4);
    expected.emplace_back(5, 5);
    expected.emplace_back(AlignedPair::gap, 7);
    expected.emplace_back(6, 9);
    expected.emplace_back(8, 10);
    expected.emplace_back(9, 12);
    expected.emplace_back(10, 13);
    expected.emplace_back(11, 14);
    
    auto got = stitcher.stitch(anchor_chain,
                               graph1, graph2,
                               tableau1, tableau2,
                               chain_merge1, chain_merge2);
    
    check_alignment(got, expected);
    
    
    cerr << "passed all tests!" << endl;
    return 0;
}
