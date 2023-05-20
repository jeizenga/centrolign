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
#include "centrolign/test_util.hpp"

using namespace std;
using namespace centrolign;

int main(int argc, char* argv[]) {
    
    Stitcher stitcher;
    
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
