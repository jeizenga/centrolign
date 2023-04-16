#include <cstdio>
#include <cstdlib>
#include <cstdlib>
#include <unordered_set>
#include <vector>
#include <list>
#include <string>

#include "centrolign/alignment.hpp"
#include "centrolign/utility.hpp"
#include "centrolign/graph.hpp"

using namespace std;
using namespace centrolign;

int main(int argc, char* argv[]) {
     
    BaseGraph graph1;
    for (char c : "ACGTGCA") {
        graph1.add_node(c);
    }
    
    BaseGraph graph2;
    for (char c : "AGTTTGA") {
        graph2.add_node(c);
    }
    
    // both graphs will have the same topology (two bubbles)
    for (auto gp : {&graph1, &graph2}) {
        auto& g = *gp;
        g.add_edge(0, 1);
        g.add_edge(0, 2);
        g.add_edge(1, 3);
        g.add_edge(2, 3);
        g.add_edge(3, 4);
        g.add_edge(3, 5);
        g.add_edge(4, 6);
        g.add_edge(5, 6);
    }
    
    AlignmentParameters<1> params;
    params.match = 1;
    params.mismatch = 1;
    params.gap_extend[0] = 1;
    params.gap_open[0] = 1;
    
    Alignment alignment = po_poa(graph1, graph2, {0}, {0}, {6}, {6}, params);
    Alignment expected;
    expected.emplace_back(0, 0);
    expected.emplace_back(2, 1);
    expected.emplace_back(3, 3);
    expected.emplace_back(4, 5);
    expected.emplace_back(6, 6);
    
    if (alignment != expected) {
        cerr << "got:\n";
        for (auto& p : alignment) {
            cerr << '\t' << p.node_id1 << ", " << p.node_id2 << '\n';
        }
        cerr << "expected:\n";
        for (auto& p : expected) {
            cerr << '\t' << p.node_id1 << ", " << p.node_id2 << '\n';
        }
        exit(1);
    }
    
    cerr << "passed all tests!" << endl;
}
