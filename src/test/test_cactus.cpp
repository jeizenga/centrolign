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
#include "centrolign/cactus.hpp"

#include "centrolign/graph.hpp"
#include "centrolign/test_util.hpp"
#include "centrolign/bridges.hpp"
#include "centrolign/connected_components.hpp"
#include "centrolign/three_edge_connected_components.hpp"


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

bool are_three_edge_connected(const BaseGraph& graph, uint64_t node_id1, uint64_t node_id2) {
    
    std::vector<std::pair<uint64_t, uint64_t>> edges;
    
    for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
        
        for (auto next_id : graph.next(node_id)) {
            edges.emplace_back(node_id, next_id);
        }
    }
    
    if (edges.size() < 2) {
        return false;
    }
    
    for (size_t i = 0; i < edges.size(); ++i) {
        for (size_t j = i + 1; j < edges.size(); ++j) {
            
            std::unordered_set<std::pair<uint64_t, uint64_t>> masked_edges;
            masked_edges.insert(edges[i]);
            masked_edges.insert(edges[j]);
            
            BaseGraph copy;
            for (uint64_t n = 0;  n < graph.node_size(); ++n) {
                copy.add_node(graph.label(n));
            }
            for (uint64_t n = 0;  n < graph.node_size(); ++n) {
                for (auto m : graph.next(n)) {
                    if (!masked_edges.count(std::make_pair(n, m))) {
                        copy.add_edge(n, m);
                    }
                }
            }
            
            // get weakly connected components
            auto comps = connected_components(copy);
            for (auto& comp : comps) {
                int count = 0;
                for (auto n : comp) {
                    if (n == node_id1 || n == node_id2) {
                        ++count;
                    }
                }
                if (count == 1) {
                    // they are in two separate components
                    return false;
                }
            }
        }
    }
    
    return true;
}

std::vector<std::vector<uint64_t>> three_connected_components_brute_force(const BaseGraph& graph) {
    
    std::vector<std::vector<uint64_t>> comps;
    
    for (uint64_t n = 0; n < graph.node_size(); ++n) {
        
        bool found = false;
        for (auto& comp : comps) {
            if (are_three_edge_connected(graph, n, comp.front())) {
                comp.push_back(n);
                found = true;
                break;
            }
        }
        
        if (!found) {
            comps.emplace_back();
            comps.back().push_back(n);
        }
    }
    
    return comps;
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
        
        exit(1);
    }
    
}

void test_three_connected_components(const BaseGraph& graph) {
    
    auto expected = three_connected_components_brute_force(graph);
    auto got = three_edge_connected_components(graph);
    
    for (auto& c : got) {
        std::sort(c.begin(), c.end());
    }
    for (auto& c : expected) {
        std::sort(c.begin(), c.end());
    }
    std::sort(got.begin(), got.end());
    std::sort(expected.begin(), expected.end());
    
    if (got != expected) {
        std::cerr << "failure getting expected three edge connected components\n";
        
        std::cerr << cpp_representation(graph, "graph") << '\n';
        
        std::cerr << "expected:\n";
        for (auto c : expected) {
            for (auto n : c) {
                std::cerr << n << ' ';
            }
            std::cerr << '\n';
        }
        std::cerr << "got:\n";
        for (auto c : got) {
            for (auto n : c) {
                std::cerr << n << ' ';
            }
            std::cerr << '\n';
        }
        
        exit(1);
    }
}

int main(int argc, char* argv[]) {
     
    random_device rd;
    default_random_engine gen(rd());
    
    {
        BaseGraph graph;
        for (int i = 0; i < 11; ++i) {
            graph.add_node('A');
        }
        graph.add_edge(0, 1);
        graph.add_edge(1, 2);
        graph.add_edge(2, 3);
        graph.add_edge(2, 6);
        graph.add_edge(3, 4);
        graph.add_edge(3, 8);
        graph.add_edge(4, 5);
        graph.add_edge(5, 6);
        graph.add_edge(6, 7);
        graph.add_edge(7, 8);
        graph.add_edge(8, 9);
        graph.add_edge(9, 10);
        
        CactusGraph<BaseGraph> cactus(graph);
        
        assert(cactus.node_size() == 3);
        
        std::vector<uint64_t>
        p0{0, 1, 2},
        p1{3},
        p2{4, 5},
        p3{6, 7},
        p4{8, 9, 10};
        
        std::pair<uint64_t, size_t> e0, e1, e2, e3, e4;
        bool f0, f1, f2, f3, f4;
        f0 = f1 = f2 = f3 = f4 = false;
        for (uint64_t n = 0; n < cactus.node_size(); ++n) {
            for (size_t i = 0; i < cactus.next_size(n); ++i) {
                auto p = cactus.next_edge_label(n, i);
                if (p == p0) {
                    assert(!f0);
                    e0 = make_pair(n, i);
                    f0 = true;
                }
                else if (p == p1) {
                    assert(!f1);
                    e1 = make_pair(n, i);
                    f1 = true;
                }
                else if (p == p2) {
                    assert(!f2);
                    e2 = make_pair(n, i);
                    f2 = true;
                }
                else if (p == p3) {
                    assert(!f3);
                    e3 = make_pair(n, i);
                    f3 = true;
                }
                else if (p == p4) {
                    assert(!f4);
                    e4 = make_pair(n, i);
                    f4 = true;
                }
            }
        }
        
        assert(f0);
        assert(f1);
        assert(f2);
        assert(f3);
        assert(f4);
        
        
        assert(e0.first != e1.first);
        assert(e1.first == e2.first);
        assert(e1.first == e3.first);
        assert(e1.first == e4.first);
        
        assert(cactus.next(e0.first)[e0.second] == e1.first);
        assert(cactus.next(e1.first)[e1.second] == e1.first);
        assert(cactus.next(e2.first)[e2.second] == e1.first);
        assert(cactus.next(e3.first)[e3.second] == e1.first);
        assert(cactus.next(e4.first)[e4.second] != e0.first);
        assert(cactus.next(e4.first)[e4.second] != e1.first);
    }
    
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
        test_three_connected_components(graph);
    }
    
    // randomized tests for bridges
    {
        vector<pair<int, int>> graph_sizes{{7, 15}, {10, 30}, {16, 70}};
        for (auto size : graph_sizes) {
            for (int rep = 0; rep < 20; ++rep) {
                BaseGraph graph = random_graph(size.first, size.second, false, gen);
                test_bridges(graph);
            }
        }
        vector<int> challenge_graph_sizes{10, 20, 30};
        for (auto size : challenge_graph_sizes) {
            for (int rep = 0; rep < 20; ++rep) {
                BaseGraph graph = random_challenge_graph(size, gen);
                test_bridges(graph);
            }
        }
    }
    // randomized tests for 3-connected components
    {
        vector<pair<int, int>> graph_sizes{{7, 15}, {10, 30}};
        for (auto size : graph_sizes) {
            for (int rep = 0; rep < 5; ++rep) {
                BaseGraph graph = random_graph(size.first, size.second, false, gen);
                test_three_connected_components(graph);
            }
        }
        vector<int> challenge_graph_sizes{10, 15, 20};
        for (auto size : challenge_graph_sizes) {
            for (int rep = 0; rep < 5; ++rep) {
                BaseGraph graph = random_challenge_graph(size, gen);
                test_three_connected_components(graph);
            }
        }
    }
    
    cerr << "passed all tests!" << endl;
}
