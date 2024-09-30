#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <utility>
#include <string>

#include "centrolign/inconsistency_identifier.hpp"
#include "centrolign/graph.hpp"
#include "centrolign/modify_graph.hpp"
#include "centrolign/test_util.hpp"

using namespace std;
using namespace centrolign;




int main(int argc, char* argv[]) {
    
    random_device rd;
    default_random_engine gen(rd());
    
    
    {
        BaseGraph graph;
        for (int i = 0; i < 8; ++i) {
            graph.add_node('A');
        }
        graph.add_edge(0, 1);
        graph.add_edge(1, 2);
        graph.add_edge(1, 3);
        graph.add_edge(2, 4);
        graph.add_edge(3, 4);
        graph.add_edge(4, 5);
        graph.add_edge(5, 5);
        graph.add_edge(5, 6);
        graph.add_edge(6, 1);
        graph.add_edge(6, 7);
        
        auto p = graph.add_path("1");
        vector<int> path{0, 1, 2, 4, 5, 6, 1, 3, 4, 5, 5, 6, 7};
        for (auto n : path) {
            graph.extend_path(p, n);
        }
        
        auto tableau = add_sentinels(graph, '^', '$');
        
        InconsistencyIdentifier inconsistency_identifier;
        
        {
            inconsistency_identifier.max_tight_cycle_size = 4;
            auto got = inconsistency_identifier.identify_tight_cycles(graph, tableau);
            std::vector<std::pair<uint64_t, uint64_t>> expected{{4, 6}};
            
            assert(got == expected);
        }
        {
            inconsistency_identifier.max_tight_cycle_size = 20;
            auto got = inconsistency_identifier.identify_tight_cycles(graph, tableau);
            std::vector<std::pair<uint64_t, uint64_t>> expected{{0, 7}};
            
            assert(got == expected);
        }
    }
    
    
    
    cerr << "passed all tests!" << endl;
    return 0;
}
