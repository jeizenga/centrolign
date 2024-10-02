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

class TestInconsistencyIdentifier : public InconsistencyIdentifier {
public:
    TestInconsistencyIdentifier() : InconsistencyIdentifier() {}
    using InconsistencyIdentifier::identify_tight_cycles;
    using InconsistencyIdentifier::identify_inconsistent_bonds;
};


int main(int argc, char* argv[]) {
    
    random_device rd;
    default_random_engine gen(rd());
    
    {
        BaseGraph graph;
        for (int i = 0; i < 14; ++i) {
            graph.add_node('A');
            if (i > 0) {
                graph.add_edge(i - 1, i);
            }
        }
        graph.add_edge(2, 4);
        graph.add_edge(4, 9);
        graph.add_edge(6, 11);
        graph.add_edge(12, 1);
        
        
        vector<vector<int>> paths{
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13},
            {0, 1, 2, 3, 4, 9, 10, 11, 12, 1, 2, 4, 5, 6, 11, 12, 13}
        };
        for (int i = 0; i < paths.size(); ++i) {
            auto p = graph.add_path(to_string(i));
            for (auto n : paths[i]) {
                graph.extend_path(p, n);
            }
        }
        
        auto tableau = add_sentinels(graph, '^', '$');
        
        TestInconsistencyIdentifier inconsistency_identifier;
        
        SnarlTree snarls(graph, tableau);
        StepIndex step_index(graph);
        
        // label the nodes that can be boundaries of non-trivial snarls
        std::vector<bool> nontrivial_left(graph.node_size(), false);
        {
            CompactedGraph compacted_graph(graph);
            for (uint64_t node_id = 0; node_id < compacted_graph.node_size(); ++node_id) {
                nontrivial_left[compacted_graph.back(node_id)] = true;
            }
        }
        
        {
            inconsistency_identifier.max_bond_inconsistency_window = 3;
            inconsistency_identifier.min_inconsistency_disjoint_length = 2;
            inconsistency_identifier.min_inconsistency_total_length = 2;
            
            auto got = inconsistency_identifier.identify_inconsistent_bonds(snarls, step_index, nontrivial_left);
            
            std::cerr << "got:\n";
            for (auto g : got) {
                std::cerr << g.first << ' ' << g.second << '\n';
            }
        }
    }
    
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
        
        TestInconsistencyIdentifier inconsistency_identifier;
        
        SnarlTree snarls(graph, tableau);
        StepIndex step_index(graph);
        std::vector<bool> nontrivial_left(graph.node_size(), true);
        
        {
            inconsistency_identifier.max_tight_cycle_size = 4;
            auto got = inconsistency_identifier.identify_tight_cycles(snarls, step_index, nontrivial_left);
            std::vector<std::pair<uint64_t, uint64_t>> expected{{4, 6}};
            
            assert(got == expected);
        }
        {
            inconsistency_identifier.max_tight_cycle_size = 20;
            auto got = inconsistency_identifier.identify_tight_cycles(snarls, step_index, nontrivial_left);
            std::vector<std::pair<uint64_t, uint64_t>> expected{{0, 7}};
            
            assert(got == expected);
        }
    }
    
    
    
    cerr << "passed all tests!" << endl;
    return 0;
}
