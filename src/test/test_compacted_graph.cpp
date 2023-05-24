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
#include <iostream>

#include "centrolign/compacted_graph.hpp"
#include "centrolign/graph.hpp"

using namespace std;
using namespace centrolign;

int main(int argc, char* argv[]) {
     
    random_device rd;
    default_random_engine gen(rd());
    
    {
        BaseGraph graph;
        for (char c : string("AAAAAAAA")) {
            graph.add_node(c);
        }
        vector<pair<int, int>> edges {
            {0, 1},
            {1, 2},
            {1, 4},
            {2, 3},
            {3, 5},
            {4, 5},
            {5, 6},
            {6, 7}
        };
        for (auto p : edges) {
            graph.add_edge(p.first, p.second);
        }
        
        CompactedGraph compacted(graph);
        
        assert(compacted.node_size() == 4);
        
        std::vector<uint64_t> cnode(4, -1);
        for (uint64_t n = 0; n < compacted.node_size(); ++n) {
            if (compacted.front(n) == 0 && compacted.back(n) == 1 && compacted.label_size(n) == 2) {
                cnode[0] = n;
            }
            else if (compacted.front(n) == 2 && compacted.back(n) == 3 && compacted.label_size(n) == 2) {
                cnode[1] = n;
            }
            else if (compacted.front(n) == 4 && compacted.back(n) == 4 && compacted.label_size(n) == 1) {
                cnode[2] = n;
            }
            else if (compacted.front(n) == 5 && compacted.back(n) == 7 && compacted.label_size(n) == 3) {
                cnode[3] = n;
            }
            else {
                assert(false);
            }
        }
        
        vector<vector<uint64_t>> nexts{
            {cnode[1], cnode[2]},
            {cnode[3]},
            {cnode[3]},
            {}
        };
        for (size_t i = 0; i < nexts.size(); ++i) {
            auto cnexts = compacted.next(cnode[i]);
            sort(nexts[i].begin(), nexts[i].end());
            sort(cnexts.begin(), cnexts.end());
            assert(nexts[i] == cnexts);
        }
    }
    
    cerr << "passed all tests!" << endl;
}
