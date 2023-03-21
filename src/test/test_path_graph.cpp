#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <utility>

#include "centrolign/path_graph.hpp"
#include "centrolign/graph.hpp"

using namespace std;
using namespace centrolign;

int main(int argc, char* argv[]) {
    
    // an extended example that i worked out the answer to:
    
    BaseGraph graph;
    // component 1
    graph.add_node('!');
    graph.add_node('T');
    graph.add_node('A');
    graph.add_node('A');
    graph.add_node('$');
    graph.add_edge(0, 1);
    graph.add_edge(1, 2);
    graph.add_edge(1, 3);
    graph.add_edge(2, 3);
    graph.add_edge(3, 4);
    // component 2
    graph.add_node('#');
    graph.add_node('C');
    graph.add_node('T');
    graph.add_node('G');
    graph.add_node('A');
    graph.add_node('%');
    graph.add_edge(5, 6);
    graph.add_edge(6, 7);
    graph.add_edge(6, 8);
    graph.add_edge(7, 9);
    graph.add_edge(8, 9);
    graph.add_edge(9, 10);
    
    PathGraph path_graph(graph);
    assert(!path_graph.is_prefix_sorted());
    assert(path_graph.node_size() == 12);
    // first doubling
    path_graph = PathGraph(path_graph);
    assert(!path_graph.is_prefix_sorted());
    assert(path_graph.node_size() == 12);
    path_graph = PathGraph(path_graph);
    // second doubling
    assert(path_graph.is_prefix_sorted());
    assert(path_graph.node_size() == 12);
    
    path_graph.order_by_rank();
    path_graph.construct_edges(graph);
    
    uint64_t nid = 0;
    assert(path_graph.from(nid) == 0);
    assert(path_graph.next_size(nid) == 2);
    assert(path_graph.next(nid)[0] == 9 || path_graph.next(nid)[1] == 9);
    assert(path_graph.next(nid)[0] == 11 || path_graph.next(nid)[1] == 11);
    nid = 1;
    assert(path_graph.from(nid) == 5);
    assert(path_graph.next_size(nid) == 1);
    assert(path_graph.next(nid)[0] == 7);
    nid = 2;
    assert(path_graph.from(nid) == 4);
    assert(path_graph.next_size(nid) == 0);
    nid = 3;
    assert(path_graph.from(nid) == 10);
    assert(path_graph.next_size(nid) == 0);
    nid = 4;
    assert(path_graph.from(nid) == 3);
    assert(path_graph.next_size(nid) == 1);
    assert(path_graph.next(nid)[0] == 2);
    nid = 5;
    assert(path_graph.from(nid) == 9);
    assert(path_graph.next_size(nid) == 1);
    assert(path_graph.next(nid)[0] == 3);
    nid = 6;
    assert(path_graph.from(nid) == 2);
    assert(path_graph.next_size(nid) == 1);
    assert(path_graph.next(nid)[0] == 4);
    nid = 7;
    assert(path_graph.from(nid) == 6);
    assert(path_graph.next_size(nid) == 2);
    assert(path_graph.next(nid)[0] == 8 || path_graph.next(nid)[1] == 8);
    assert(path_graph.next(nid)[0] == 10 || path_graph.next(nid)[1] == 10);
    nid = 8;
    assert(path_graph.from(nid) == 8);
    assert(path_graph.next_size(nid) == 1);
    assert(path_graph.next(nid)[0] == 5);
    nid = 9;
    assert(path_graph.from(nid) == 1);
    assert(path_graph.next_size(nid) == 1);
    assert(path_graph.next(nid)[0] == 4);
    nid = 10;
    assert(path_graph.from(nid) == 7);
    assert(path_graph.next_size(nid) == 1);
    assert(path_graph.next(nid)[0] == 5);
    nid = 11;
    assert(path_graph.from(nid) == 1);
    assert(path_graph.next_size(nid) == 1);
    assert(path_graph.next(nid)[0] == 6);
    
    cerr << "passed all tests!" << endl;
    return 0;
}
