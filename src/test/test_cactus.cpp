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
#include <cmath>
#include <iostream>

#include "centrolign/adjacency_graph.hpp"

#include "centrolign/graph.hpp"
#include "centrolign/test_util.hpp"
#include "centrolign/bridges.hpp"


using namespace std;
using namespace centrolign;

std::vector<std::pair<uint64_t, uint64_t>> find_bridges_brute_force(const BaseGraph& graph) {
    
    std::vector<std::pair<uint64_t, uint64_t>> bridges_found;
    
    for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
        
        for (uint64_t next_id : graph.next(node_id)) {
            
            // create a copy with this edge removed
            BaseGraph copy;
            bool found_edge = false;
            for (uint64_t copy_id = 0; copy_id < graph.node_size(); ++copy_id) {
                copy.add_node(graph.label(copy_id));
            }
            for (uint64_t copy_id = 0; copy_id < graph.node_size(); ++copy_id) {
                for (uint64_t next_copy_id : graph.next(copy_id)) {
                    if (!found_edge && copy_id == node_id && next_copy_id == next_id) {
                        found_edge = true;
                        continue;
                    }
                    // add both direction edges so the reachability is equivalent to
                    // an undirected graph
                    copy.add_edge(copy_id, next_copy_id);
                    copy.add_edge(next_copy_id, copy_id);
                }
            }
            
            if (!is_reachable(copy, node_id, next_id)) {
                bridges_found.emplace_back(node_id, next_id);
            }
        }
    }
    
    return bridges_found;
}

void test_bridges(const BaseGraph& graph) {
    
    auto expected = find_bridges_brute_force(graph);
    auto got = bridges(graph);
    
    std::sort(got.begin(), got.end());
    std::sort(expected.begin(), expected.end());
    
    if (got != expected) {
        std::cerr << "failure getting expected bridges\n";
        
        std::cerr << cpp_representation(graph, "graph") << '\n';
        
        std::cerr << "expected:\n";
        for (auto e : expected) {
            std::cerr << '\t' << e.first << '\t' << e.second << '\n';
        }
        std::cerr << "got:\n";
        for (auto e : got) {
            std::cerr << '\t' << e.first << '\t' << e.second << '\n';
        }
    }
    
}


int main(int argc, char* argv[]) {
     
    random_device rd;
    default_random_engine gen(rd());
    
    {
        
        BaseGraph graph;
        
        graph.add_node('A');
        graph.add_node('A');
        graph.add_node('A');
        graph.add_node('A');
        graph.add_node('A');
        
        graph.add_edge(0, 1);
        graph.add_edge(0, 3);
        graph.add_edge(1, 2);
        graph.add_edge(1, 4);
        graph.add_edge(2, 3);
        graph.add_edge(3, 4);
        
        AdjacencyGraph adj_graph(graph);
        
        assert(adj_graph.node_size() == 4);
        
        uint64_t n0 = -1, n1 = -1, n2 = -1, n3 = -1;
        for (uint64_t n = 0; n < adj_graph.node_size(); ++n) {
            if (adj_graph.previous_size(n) == 0 && adj_graph.next_size(n) == 1) {
                n0 = n;
            }
            else if (adj_graph.next_size(n) == 0 && adj_graph.previous_size(n) == 1) {
                n3 = n;
            }
        }
        
        for (uint64_t n = 0; n < adj_graph.node_size(); ++n) {
            if (adj_graph.next_size(n) == 2 && adj_graph.previous_size(n) == 2) {
                bool found = false;
                for (const auto& e : adj_graph.next_edges(n)) {
                    if (e.target == n3) {
                        found = true;
                        break;
                    }
                }
                if (found) {
                    n2 = n;
                }
                else {
                    n1 = n;
                }
            }
        }
        
        assert(n0 != -1);
        assert(n1 != -1);
        assert(n2 != -1);
        assert(n3 != -1);
        
        {
            assert(adj_graph.previous_size(n0) == 0);
            assert(adj_graph.next_size(n0) == 1);
            assert(adj_graph.next_edges(n0).front().target == n1);
            assert(adj_graph.next_edges(n0).front().label == 0);
        }
        
        {
            assert(adj_graph.previous_size(n1) == 2);
            assert(adj_graph.next_size(n1) == 2);
            {
                bool found1 = false, found2 = false;
                for (auto e : adj_graph.previous_edges(n1)) {
                    if (e.label == 0 && e.target == n0) {
                        found1 = true;
                    }
                    else if (e.label == 2 && e.target == n2) {
                        found2 = true;
                    }
                    else {
                        assert(false);
                    }
                }
                assert(found1);
                assert(found2);
            }
            {
                bool found1 = false, found2 = false;
                for (auto e : adj_graph.next_edges(n1)) {
                    assert(e.target == n2);
                    if (e.label == 1) {
                        found1 = true;
                    }
                    else if (e.label == 3) {
                        found2 = true;
                    }
                    else {
                        assert(false);
                    }
                }
                assert(found1);
                assert(found2);
            }
        }
        
        {
            assert(adj_graph.previous_size(n2) == 2);
            assert(adj_graph.next_size(n2) == 2);
            {
                bool found1 = false, found2 = false;
                for (auto e : adj_graph.previous_edges(n2)) {
                    assert(e.target == n1);
                    if (e.label == 1) {
                        found1 = true;
                    }
                    else if (e.label == 3) {
                        found2 = true;
                    }
                    else {
                        assert(false);
                    }
                }
                assert(found1);
                assert(found2);
            }
            {
                bool found1 = false, found2 = false;
                for (auto e : adj_graph.next_edges(n2)) {
                    if (e.label == 2 && e.target == n1) {
                        found1 = true;
                    }
                    else if (e.label == 4 && e.target == n3) {
                        found2 = true;
                    }
                    else {
                        assert(false);
                    }
                }
                assert(found1);
                assert(found2);
            }
        }
        
        {
            assert(adj_graph.previous_size(n3) == 1);
            assert(adj_graph.next_size(n3) == 0);
            assert(adj_graph.previous_edges(n3).front().target == n2);
            assert(adj_graph.previous_edges(n3).front().label == 4);
        }
    }
    
    {
        
        BaseGraph graph;
        for (size_t i = 0; i < 8; ++i) {
            graph.add_node('A');
        }
        graph.add_edge(1, 0);
        graph.add_edge(1, 2);
        graph.add_edge(1, 3);
        graph.add_edge(2, 4);
        graph.add_edge(3, 4);
        graph.add_edge(4, 5);
        graph.add_edge(5, 6);
        graph.add_edge(5, 7);
        graph.add_edge(6, 7);
        
        test_bridges(graph);
    }
    
    vector<pair<int, int>> graph_sizes{{7, 15}, {10, 30}, {16, 70}};
    for (auto size : graph_sizes) {
        for (int rep = 0; rep < 5; ++rep) {
            BaseGraph graph = random_graph(size.first, size.second, false, gen);
            test_bridges(graph);
        }
    }
    vector<int> challenge_graph_sizes{10, 20, 30};
    for (auto size : challenge_graph_sizes) {
        for (int rep = 0; rep < 5; ++rep) {
            BaseGraph graph = random_challenge_graph(size, gen);
            test_bridges(graph);
        }
    }
    
    cerr << "passed all tests!" << endl;
}
