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
        
        stringstream strm;
        write_gfa(graph, strm, false);
        
        string gfa = strm.str();
        // just verify that it has the expected number of nodes and edges
        assert(count(gfa.begin(), gfa.end(), 'S') == 7);
        assert(count(gfa.begin(), gfa.end(), 'L') == 7);
        std::cerr << gfa;
    }
    
    cerr << "passed all tests!" << endl;
}
