#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <unordered_set>
#include <vector>
#include <list>
#include <string>
#include <random>
#include <sstream>
#include <algorithm>

#include "centrolign/gfa.hpp"
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
        for (char c : string("ACGTACGTAC")) {
            graph.add_node(c);
        }
        vector<pair<int, int>> edges {
            {1, 0},
            {0, 2},
            {2, 3},
            {2, 4},
            {3, 6},
            {4, 5},
            {5, 6},
            {6, 7},
            {7, 8},
            {8, 9}
        };
        for (auto p : edges) {
            graph.add_edge(p.first, p.second);
        }
        vector<vector<int>> paths{
            {2, 3, 6, 7},
            {5, 6, 7}
        };
        
        int i = 0;
        for (auto p : paths) {
            auto pid = graph.add_path(to_string(i++));
            for (auto j : p) {
                graph.extend_path(pid, j);
            }
        }
        
        string gfa;
        {
            stringstream strm;
            write_gfa(graph, strm, false);
            
            gfa = strm.str();
            // just verify that it has the expected number of nodes and edges
            assert(count(gfa.begin(), gfa.end(), 'S') == 7);
            assert(count(gfa.begin(), gfa.end(), 'L') == 7);
            //std::cerr << gfa;
        }
        
        auto tableau = add_sentinels(graph, '#', '$');
        
        string gfa_sentinels;
        {
            stringstream strm;
            write_gfa(graph, tableau, strm, false);
            
            gfa_sentinels = strm.str();
            // just verify that it has the expected number of nodes and edges
            assert(count(gfa_sentinels.begin(), gfa_sentinels.end(), 'S') == 7);
            assert(count(gfa_sentinels.begin(), gfa_sentinels.end(), 'L') == 7);
            //std::cerr << gfa;
        }
        
        assert(gfa == gfa_sentinels);
    }
    
    int num_reps = 10;
    int size = 100;
    for (int rep = 0; rep < num_reps; ++rep) {
        
        auto graph = random_challenge_graph(size, gen);
        
        stringstream strm_out;
        write_gfa(graph, strm_out, false);
        
        auto gfa = strm_out.str();
        
        stringstream strm_in(gfa);
        
        auto loaded_graph = read_gfa(strm_in, false);
        
        if (!possibly_isomorphic(graph, loaded_graph)) {
            cerr << "failed re-load test on graph:\n";
            cerr << cpp_representation(graph, "graph") << '\n';
            
            cerr << "got this graph instead:\n";
            cerr << cpp_representation(loaded_graph, "loaded") << '\n';
            exit(1);
        }
    }
    
    
    cerr << "passed all tests!" << endl;
}
