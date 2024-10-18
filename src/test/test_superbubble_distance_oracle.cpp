#include <cstdio>
#include <cstdlib>
#include <cstdlib>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <utility>
#include <cassert>
#include <limits>

#include "centrolign/utility.hpp"
#include "centrolign/test_util.hpp"
#include "centrolign/graph.hpp"
#include "centrolign/minmax_distance.hpp"
#include "centrolign/superbubble_distance_oracle.hpp"

using namespace std;
using namespace centrolign;



void do_tests(const BaseGraph& graph) {
    
    SuperbubbleDistanceOracle oracle(graph);
    
    for (uint64_t n = 0; n < graph.node_size(); ++n) {
        
        vector<uint64_t> source{n};
        auto minmax = minmax_distance(graph, &source);
        
        for (uint64_t m = 0; m < graph.node_size(); ++m) {
            
            size_t d = oracle.min_distance(n, m);
            bool match = d == -1 ? (minmax[m].first == numeric_limits<int64_t>::max()) : (minmax[m].first == d);
            if (!match) {
                
                cerr << "failed distance test between " << n << " and " << m << ", expected " << minmax[m].first << ", got " << d << '\n';
                cerr << "graph:\n";
                cerr << cpp_representation(graph, "graph") << '\n';
                exit(1);
            }
        }
    }
}

int main(int argc, char* argv[]) {
     
    random_device rd;
    default_random_engine gen(rd());

    
    {
        BaseGraph graph;
        
        for (auto c : string("ACGTACGTACGTACG")) {
            graph.add_node(c);
        }
        vector<pair<int, int>> edges{
            {0, 1},
            {0, 7},
            {0, 14},
            {1, 2},
            {2, 3},
            {2, 4},
            {3, 5},
            {4, 5},
            {5, 6},
            {6, 9},
            {7, 8},
            {8, 9},
            {9, 10},
            {9, 11},
            {10, 11},
            {11, 12},
            {11, 13},
            {14, 9}
        };
        for (auto e : edges) {
            graph.add_edge(e.first, e.second);
        }
        
        do_tests(graph);
    }
    

    size_t num_reps = 50;
    vector<pair<size_t, size_t>> graph_sizes{
        {10, 12},
        {20, 30}
    };
    for (auto& sizes : graph_sizes) {
        size_t num_nodes = sizes.first;
        size_t num_edges = sizes.second;
        for (size_t i = 0; i < num_reps; ++i) {
            BaseGraph graph = random_graph(num_nodes, num_edges, true, gen);
            do_tests(graph);
            BaseGraph hard_graph = random_challenge_graph(num_nodes, gen);
            do_tests(hard_graph);
        }
    }
    
    
    cerr << "passed all tests!" << endl;
}
