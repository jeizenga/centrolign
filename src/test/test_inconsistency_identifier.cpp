#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <utility>
#include <string>
#include <sstream>

#include "centrolign/inconsistency_identifier.hpp"
#include "centrolign/graph.hpp"
#include "centrolign/modify_graph.hpp"
#include "centrolign/test_util.hpp"
#include "centrolign/gfa.hpp"


using namespace std;
using namespace centrolign;

class TestInconsistencyIdentifier : public InconsistencyIdentifier {
public:
    TestInconsistencyIdentifier() : InconsistencyIdentifier() {}
    using InconsistencyIdentifier::identify_tight_cycles;
    using InconsistencyIdentifier::identify_inconsistent_bonds;
    using InconsistencyIdentifier::expand_inconsistencies;
};


int main(int argc, char* argv[]) {
    
    random_device rd;
    default_random_engine gen(rd());
    
    {
        std::string gfa =
        "H\tVN:Z:1.1\n"
        "S\t8294\tG\n"
        "S\t20136\tA\n"
        "S\t20137\tC\n"
        "S\t8295\tC\n"
        "S\t13753\tC\n"
        "S\t13755\tA\n"
        "S\t8309\tTGTTACACAGAAGAATTCTCAGTAACTTATTTGTGGTGTGTGTATTCAACTCACAGAGATGAACCTTCCTTCAGAAAGAGCAGATTTGAAACACTCTTTTTGTGGAGTTTCCATGTGGAGATTTCAATCGCTTTGAGACCAAAGGTAGAAAAGGAAACATCTTCGTATAACAACTAGACAGAATCATTCACAGAAACTACTTTGTGATGTGTGTGTTCAACTCAAGGAGTTTAACCTTTCTTTTGATGGAGCAGTTTGGAAACACTCTGTCTGTAAAGTCTGCAAGCAGATATTTGGACCTCTTTGAGGCCTTCGTTGGAAACGGGATTTCTTCATATAATGTTTGATAGGAGAAGTCTCAGTAACTTCTTTGTGCTGTGTGTATTCAACTCATAGAGTTGAACTTTCCTTTAGAAGAGCAGATGTTAAACACCCTTTTTGTGGAATTTGCAGCTGGAGATTTCAAGCGCTTTGAGACCTATGGTAGAAAAGGAAACATCTTCTTATAAAATCTAGACAGAATCATTCACAGAAACTTCTTTTCGATGTGTGTGTTCAGCTCACAGAGTTTAACCTTTCTTTTGATGGAGCAGTTTGGAAACACTCTGTTTGTAATGTCTGCAAGTGGATATTTGGACCTCTTTGAGGCCTTCGTTGGAAACGGGATTTCTTCAAGTAATGTTCGACAGAAGAATTCTCAGTAACTT\n"
        "S\t8305\tGAAACTACTTTGTGATGTGTGCGTTCAACTCAAGGAGTTTAAGCTTTCTTTTCATAGAGTAGTTTGGAAACACTCTGTCTGTAAAGTCTGCAAGCAGATATTTGGACCTCTTTGGGGCCTTCGTTGGAAACGGGATTTCTTCATAGAACGCTAGAAAGAAGAATACTGAGTAAGTTCTTTGTGTTGCCTCTATTCAACTCACAGAGGTGAACTGTCCTTTAGACAGAGCAGATGTGAAACCCTCTTTTTGTGATATTTGCAGGTGGAGATTTCAAGCGCTTTTAGGCCAAATGTAGAAAAGGAAATATCTTCGTATAAAAACTAGACAGAATCATTCTCAGAAACTACTTTGTGATGTGTGC\n"
        "S\t13754\tT\n"
        "S\t15938\tA\n"
        "S\t8303\tGATGTGTGTGTTCAGCTCACAGAGTTTAACCTTTCTTTTGATGGAGCAGTTTGGAAACACTCTGTTTGTAATGTCTGCAAGTGGATATTTGGACCTCTTTGAGGCCTTCGTTGGAAACGGGATTTCTTCAAGTAATGTTCGACAGAAGAATTCTCAGTAACTTATTTGTGGTGTGTGTATTCAACTCACAGAGTTGAACCTTCCTTTAGACAGAGCAGATTTGAAACACCCTATTTGTGCAGTTTCCAGTTGGAGATTTCAATCGCTTTGAGACCAAATGTAGAAAAGGAAACATCTTCGTATAAAAACTAGACAGAATCATTCTC\n"
        "S\t8297\tAGAGGTGAACTGTCCTTTAGACAGAGCAGATGTGAAACCCTCTTTTTGTGATATTTGCAGGTGGAGATTTCAAGCGCTTTTAGGCCAAATGTAGAAAAGGAAATATCTTCGTATAAAAACTAGACAGAATCATTCTCAGAAACTACTTTGTGATGAGTGCGTTCAATTCACAGAGTATAACCTTTCTTTTGATGGAGGAGTTTGGAGACACTGTCTTTGTAAAGTCTGCAAGTGGATATTTGGACCTCTTT\n"
        "S\t8306\tA\n"
        "S\t8307\tTTCAATTCACAGAGTATAACCTTTCTTTTGATGGAGGAGTTTGGAGACACTGTCTTTGTAAAGTCTGCAAGTGGATATTTGGACCTCTTTGAGGCCTTCGTTGGAAACGGGATTTCCTCATATA\n"
        "S\t8302\tC\n"
        "S\t8308\tG\n"
        "S\t15937\tT\n"
        "S\t8298\tG\n"
        "S\t15939\tG\n"
        "S\t8304\tC\n"
        "S\t8296\tG\n"
        "S\t8300\tA\n"
        "S\t8301\tAGACAGAATCATTCACAGAAACTACTTTGTGATGTGTGTGTTCAACTCAAGGAGTTTAACCTTTCTTTTGATGGAGCAGTTTGGAAACACTCTGTCTGTAAAGTCTGCAAGCAGATATTTGGACCTCTTTGAGGCCTTCGTTGGAAACGGGATTTCTTCATATAATGTTTGATAGGAGAAGTCTCAGTAACTTCTTTGTGCTGTGTGTATTCAACTCATAGAGTTGAACTTTCCTTTAGAAGAGCAGATGTTAAACACCCTTTTTGTGGAATTTGCAGCTGGAGATTTCAAGCGCTTTGAGGCCTACGGTAGAAAAGGAAACATCTTCTTATAAAATCTAGACAGAATCATTCACAGAAACTTCTTTT\n"
        "S\t8311\tTTTGTGGTGTGTGTATTCAACTCACAGAGTTGAACCTTCCTTTAGACAGAGCAGATTTGAAACTCCCTATTTGTGCAGTTTCCAGTTGGAGATTTCAATCGCTTTGAGACCAAATGTAGAAAAGGAAACATCTTCGTATAAAAACTAGACAGAATCATTCTCCGAAACTACTTTGTGATGTGTGCGTTCAACTCAAGGAGTTTAAGCTTTCTTTTCATAGAGTAGTTTGGAAACACTCTGTCTGTAAAGTCTGCAAGCAGATATTTGGACCTCTTTGGGGCCTTCGTTGGAAACGGGATTTCTTCATAGAACGCTAGAAAGAAGAATACTGAGTAAGTTCTTTGTGTTGCCTCTATTCAACTCA\n"
        "S\t8310\tA\n"
        "S\t8299\tAGGCCTTCGTTGGAAACGGGATTTCCTCATATAATGTTACACAGAAGAATTCTCAGTAACTTATTTGTGGTGTGTGTATTGAACTCACAGAGATGAACCTTCCTTCAGAAAGAGCAGATTTGAAACACTCTTTTTGTGGAGTTTCCATGTGGAGATTTCAATCGCTTTGAGACCAAAGGTAGAAAAGGAAACATCTTCGTATAACAAC\n"
        "S\t8312\tC\n"
        "P\tHG00438.1\t8294+,8295+,8297+,20136+,8299+,8300+,8301+,15937+,8303+,15938+,8305+,8306+,8307+,13755+,8309+,20137+,8311+,8296+,8297+,8298+,8299+,8300+,8301+,8302+,8303+,8304+,8305+,8306+,8307+,8308+,8309+,8310+,8311+,8312+\t*\n"
        "P\tHG00438.2_2\t8294+,8295+,8297+,8298+,8299+,8300+,8301+,8302+,8303+,8304+,8305+,8306+,8307+,13755+,8309+,8310+,8311+,8312+\t*\n"
        "P\tHG01928.2_1\t8294+,13753+,8297+,8298+,8299+,13754+,8301+,8302+,8303+,8304+,8305+,8306+,8307+,13755+,8309+,8310+,8311+,8312+\t*\n"
        "P\tHG01928.2_2\t8294+,8295+,8297+,8298+,8299+,8300+,8301+,8302+,8303+,8304+,8305+,8306+,8307+,13755+,8309+,8310+,8311+,8312+\t*\n"
        "P\tHG01978.1\t8294+,13753+,8297+,8298+,8299+,8300+,8301+,15937+,8303+,15938+,8305+,8306+,8307+,13755+,8309+,8310+,8311+,8295+,8297+,8298+,8299+,8300+,8301+,15937+,8303+,15938+,8305+,15939+,8307+,13755+,8309+,8310+,8311+,8312+\t*\n"
        "L\t8295\t+\t8297\t+\t0M\n"
        "L\t8311\t+\t8295\t+\t0M\n"
        "L\t8309\t+\t8310\t+\t0M\n"
        "L\t8309\t+\t20137\t+\t0M\n"
        "L\t13755\t+\t8309\t+\t0M\n"
        "L\t8305\t+\t8306\t+\t0M\n"
        "L\t8305\t+\t15939\t+\t0M\n"
        "L\t15938\t+\t8305\t+\t0M\n"
        "L\t8303\t+\t8304\t+\t0M\n"
        "L\t8303\t+\t15938\t+\t0M\n"
        "L\t15937\t+\t8303\t+\t0M\n"
        "L\t8294\t+\t8295\t+\t0M\n"
        "L\t8294\t+\t8296\t+\t0M\n"
        "L\t8294\t+\t13753\t+\t0M\n"
        "L\t8297\t+\t8298\t+\t0M\n"
        "L\t8297\t+\t20136\t+\t0M\n"
        "L\t13753\t+\t8297\t+\t0M\n"
        "L\t8306\t+\t8307\t+\t0M\n"
        "L\t8307\t+\t8308\t+\t0M\n"
        "L\t8307\t+\t13755\t+\t0M\n"
        "L\t15939\t+\t8307\t+\t0M\n"
        "L\t8302\t+\t8303\t+\t0M\n"
        "L\t8308\t+\t8309\t+\t0M\n"
        "L\t8298\t+\t8299\t+\t0M\n"
        "L\t8304\t+\t8305\t+\t0M\n"
        "L\t8296\t+\t8297\t+\t0M\n"
        "L\t8300\t+\t8301\t+\t0M\n"
        "L\t8301\t+\t8302\t+\t0M\n"
        "L\t8301\t+\t15937\t+\t0M\n"
        "L\t13754\t+\t8301\t+\t0M\n"
        "L\t8311\t+\t8312\t+\t0M\n"
        "L\t20137\t+\t8311\t+\t0M\n"
        "L\t8310\t+\t8311\t+\t0M\n"
        "L\t8299\t+\t8300\t+\t0M\n"
        "L\t8299\t+\t13754\t+\t0M\n"
        "L\t20136\t+\t8299\t+\t0M\n";
        
        stringstream strm(gfa);
        
        auto graph = read_gfa(strm, false);
                
        auto tableau = add_sentinels(graph, '^', '$');
        
        TestInconsistencyIdentifier inconsistency_identifier;
        
        SnarlTree snarls(graph, tableau);
        StepIndex step_index(graph);
        std::vector<bool> nontrivial_left(graph.node_size(), false);
        {
            CompactedGraph compacted_graph(graph);
            for (uint64_t node_id = 0; node_id < compacted_graph.node_size(); ++node_id) {
                nontrivial_left[compacted_graph.back(node_id)] = true;
            }
        }
        
        inconsistency_identifier.max_tight_cycle_size = 10000;
        
        auto got = inconsistency_identifier.identify_tight_cycles(snarls, step_index, nontrivial_left);
        assert(got.size() == 1);
        assert(got.front().first == graph.path(0).front());
        assert(got.front().second == graph.path(0).back());
    }
    
    
    {
        BaseGraph graph;
        for (int i = 0; i < 9; ++i) {
            graph.add_node('A');
            if (i > 0) {
                graph.add_edge(i - 1, i);
            }
        }
        graph.add_edge(1, 7);
        graph.add_edge(3, 3);
        graph.add_edge(4, 6);
        
        SentinelTableau tableau;
        tableau.src_id = 0;
        tableau.snk_id = 8;
        
        TestInconsistencyIdentifier inconsistency_identifier;
        
        SnarlTree snarls(graph, tableau);
        StepIndex step_index(graph);
        
        {
            inconsistency_identifier.padding_target_min_length = 20;
            inconsistency_identifier.padding_max_length_limit = 100;
            
            std::vector<std::pair<uint64_t, uint64_t>> inconsistencies{{4, 6}};
            
            inconsistency_identifier.expand_inconsistencies(inconsistencies, graph, snarls);
            
            assert(inconsistencies.size() == 1);
            assert(inconsistencies.front().first == 4);
            assert(inconsistencies.front().second == 6);
        }
    }
    
    {
        BaseGraph graph;
        for (int i = 0; i < 12; ++i) {
            graph.add_node('A');
            if (i > 0) {
                graph.add_edge(i - 1, i);
            }
        }
        graph.add_edge(1, 3);
        graph.add_edge(3, 5);
        graph.add_edge(6, 8);
        graph.add_edge(8, 10);
        
        SentinelTableau tableau;
        tableau.src_id = 0;
        tableau.snk_id = 11;
        
        TestInconsistencyIdentifier inconsistency_identifier;
        
        SnarlTree snarls(graph, tableau);
        StepIndex step_index(graph);
        
        {
            inconsistency_identifier.padding_target_min_length = 5;
            inconsistency_identifier.padding_max_length_limit = 10;
            
            std::vector<std::pair<uint64_t, uint64_t>> inconsistencies{{1, 3}, {8, 10}};
            
            inconsistency_identifier.expand_inconsistencies(inconsistencies, graph, snarls);
            
            assert(inconsistencies.size() == 2);
            std::sort(inconsistencies.begin(), inconsistencies.end());
            assert(inconsistencies.front().first == 1);
            assert(inconsistencies.front().second == 5);
            assert(inconsistencies.back().first == 6);
            assert(inconsistencies.back().second == 10);
        }
    }
    
    {
        BaseGraph graph;
        for (int i = 0; i < 17; ++i) {
            graph.add_node('A');
            if (i > 0) {
                graph.add_edge(i - 1, i);
            }
        }
        graph.add_edge(2, 4);
        graph.add_edge(5, 7);
        graph.add_edge(8, 15);
        
        SentinelTableau tableau;
        tableau.src_id = 0;
        tableau.snk_id = 16;
        
        TestInconsistencyIdentifier inconsistency_identifier;
        
        SnarlTree snarls(graph, tableau);
        StepIndex step_index(graph);
        
        {
            inconsistency_identifier.padding_target_min_length = 2;
            inconsistency_identifier.padding_max_length_limit = 5;
            
            std::vector<std::pair<uint64_t, uint64_t>> inconsistencies{{5, 7}};
            
            inconsistency_identifier.expand_inconsistencies(inconsistencies, graph, snarls);
            
            assert(inconsistencies.size() == 1);
            assert(inconsistencies.front().first == 2);
            assert(inconsistencies.front().second == 8);
        }
    }
    
    {
        BaseGraph graph;
        for (int i = 0; i < 12; ++i) {
            graph.add_node('A');
            if (i > 0) {
                graph.add_edge(i - 1, i);
            }
        }
        graph.add_edge(2, 6);
        graph.add_edge(4, 4);
        graph.add_edge(6, 9);
        graph.add_edge(10, 1);
        
        vector<vector<int>> paths{
            {0, 1, 2, 6, 9, 10, 11},
            {0, 1, 2, 3, 4, 4, 5, 6, 9, 10, 1, 2, 6, 7, 8, 9, 10, 11}
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
            inconsistency_identifier.max_tight_cycle_size = 10;
            inconsistency_identifier.max_bond_inconsistency_window = 4;
            inconsistency_identifier.min_inconsistency_disjoint_length = 2;
            inconsistency_identifier.min_inconsistency_total_length = 3;
            
            auto got = inconsistency_identifier.identify_inconsistent_bonds(snarls, step_index, nontrivial_left);
            
            std::vector<std::pair<uint64_t, uint64_t>> expected{{2, 9}};
            
            assert(got == expected);
        }
    }
    
    {
        BaseGraph graph;
        for (int i = 0; i < 14; ++i) {
            graph.add_node('A');
            if (i > 0) {
                graph.add_edge(i - 1, i);
            }
        }
        graph.add_edge(3, 7);
        graph.add_edge(8, 12);
        graph.add_edge(12, 1);
        
        vector<vector<int>> paths{
            {0, 1, 2, 3, 7, 8, 9, 12, 13},
            {0, 1, 2, 3, 7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13}
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
            inconsistency_identifier.max_bond_inconsistency_window = 4;
            inconsistency_identifier.min_inconsistency_disjoint_length = 2;
            inconsistency_identifier.min_inconsistency_total_length = 3;
            
            auto got = inconsistency_identifier.identify_inconsistent_bonds(snarls, step_index, nontrivial_left);
            
            std::vector<std::pair<uint64_t, uint64_t>> expected{{3, 12}};
            
            assert(got == expected);
        }
        
        {
            inconsistency_identifier.max_bond_inconsistency_window = 2;
            inconsistency_identifier.min_inconsistency_disjoint_length = 2;
            inconsistency_identifier.min_inconsistency_total_length = 3;
            
            auto got = inconsistency_identifier.identify_inconsistent_bonds(snarls, step_index, nontrivial_left);
            
            std::vector<std::pair<uint64_t, uint64_t>> expected{};
            
            assert(got == expected);
        }
    }
    
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
            inconsistency_identifier.min_inconsistency_total_length = 4;
            
            auto got = inconsistency_identifier.identify_inconsistent_bonds(snarls, step_index, nontrivial_left);
            
            std::vector<std::pair<uint64_t, uint64_t>> expected{{4, 11}};
            
            assert(got == expected);
        }
        
        {
            inconsistency_identifier.max_bond_inconsistency_window = 3;
            inconsistency_identifier.min_inconsistency_disjoint_length = 2;
            inconsistency_identifier.min_inconsistency_total_length = 6;
            
            auto got = inconsistency_identifier.identify_inconsistent_bonds(snarls, step_index, nontrivial_left);
            
            std::vector<std::pair<uint64_t, uint64_t>> expected{};
            
            assert(got == expected);
        }
        
        {
            inconsistency_identifier.max_bond_inconsistency_window = 3;
            inconsistency_identifier.min_inconsistency_disjoint_length = 3;
            inconsistency_identifier.min_inconsistency_total_length = 3;
            
            auto got = inconsistency_identifier.identify_inconsistent_bonds(snarls, step_index, nontrivial_left);
            
            std::vector<std::pair<uint64_t, uint64_t>> expected{};
            
            assert(got == expected);
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
