#include <cstdio>
#include <cstdlib>
#include <cstdlib>
#include <unordered_set>
#include <vector>
#include <list>
#include <string>

#include "centrolign/chain_merge.hpp"
#include "centrolign/random_graph.hpp"
#include "centrolign/utility.hpp"

using namespace std;
using namespace centrolign;

// DFS
bool is_reachable(const BaseGraph& graph, uint64_t id_from, uint64_t id_to) {
    
    vector<bool> traversed(graph.node_size(), false);
    vector<uint64_t> stack(1, id_from);
    traversed[id_from] = true;
    while (!stack.empty()) {
        auto here = stack.back();
        stack.pop_back();
        for (auto nid : graph.next(here)) {
            if (nid == id_to) {
                return true;
            }
            if (!traversed[nid]) {
                traversed[nid] = true;
                stack.push_back(nid);
            }
        }
    }
    return false;
}

void make_path_cover(BaseGraph& graph, default_random_engine& gen) {
    
    uint64_t num_covered = 0;
    vector<bool> covered(graph.node_size(), false);
    
    while (num_covered < graph.node_size()) {
        size_t rank = uniform_int_distribution<size_t>(0, graph.node_size() - num_covered - 1)(gen);
        uint64_t seed = -1;
        for (uint64_t node_id = 0, r = 0; node_id < graph.node_size(); ++node_id) {
            if (!covered[node_id]) {
                if (r == rank) {
                    seed = node_id;
                    break;
                }
                ++r;
            }
        }
        assert(seed != -1);
        
        list<uint64_t> path{seed};
        while (graph.previous_size(path.front()) != 0) {
            vector<uint64_t> uncovered_prevs;
            for (auto p : graph.previous(path.front())) {
                if (!covered[p]) {
                    uncovered_prevs.push_back(p);
                }
            }
            uint64_t nid;
            if (uncovered_prevs.empty()) {
                size_t idx = uniform_int_distribution<size_t>(0, graph.previous_size(path.front()) - 1)(gen);
                nid = graph.previous(path.front())[idx];
            }
            else {
                size_t idx = uniform_int_distribution<size_t>(0, uncovered_prevs.size() - 1)(gen);
                nid = uncovered_prevs[idx];
            }
            path.push_front(nid);
        }
        while (graph.next_size(path.back()) != 0) {
            vector<uint64_t> uncovered_nexts;
            for (auto n : graph.next(path.back())) {
                if (!covered[n]) {
                    uncovered_nexts.push_back(n);
                }
            }
            uint64_t nid;
            if (uncovered_nexts.empty()) {
                size_t idx = uniform_int_distribution<size_t>(0, graph.next_size(path.back()) - 1)(gen);
                nid = graph.next(path.back())[idx];
            }
            else {
                size_t idx = uniform_int_distribution<size_t>(0, uncovered_nexts.size() - 1)(gen);
                nid = uncovered_nexts[idx];
            }
            path.push_back(nid);
        }
        
        uint64_t path_id = graph.add_path(to_string(num_covered));
        for (auto nid : path) {
            if (!covered[nid]) {
                covered[nid] = true;
                ++num_covered;
            }
            graph.extend_path(path_id, nid);
        }
    }
}

void do_tests(BaseGraph& graph) {
    
    ChainMerge chain_merge(graph);
    
    for (uint64_t id1 = 0; id1 < graph.node_size(); ++id1) {
        for (uint64_t id2 = 0; id2 < graph.node_size(); ++id2) {
            
            bool expected = is_reachable(graph, id1, id2);
            bool actual = chain_merge.reachable(id1, id2);
            if (expected != actual) {
                cerr << "got reachability " << actual << " when expecting " << expected << " between " << id1 << " and " << id2 << " on graph\n";
                print_graph(graph, cerr);
                exit(1);
            }
            
        }
    }
}


int main(int argc, char* argv[]) {
     
    random_device rd;
    default_random_engine gen(rd());
    
    {
        BaseGraph graph;
        for (size_t i = 0; i < 5; ++i) {
            graph.add_node('C');
        }
        graph.add_edge(0, 1);
        graph.add_edge(0, 2);
        graph.add_edge(0, 4);
        graph.add_edge(1, 2);
        graph.add_edge(1, 3);
        graph.add_edge(2, 3);
        graph.add_edge(3, 4);
        
        graph.add_path("1");
        graph.add_path("2");
        graph.extend_path(0, 0);
        graph.extend_path(0, 4);
        graph.extend_path(1, 0);
        graph.extend_path(1, 1);
        graph.extend_path(1, 2);
        graph.extend_path(1, 3);
        graph.extend_path(1, 4);
        
        do_tests(graph);
    }
    
    size_t num_reps = 10;
    vector<pair<size_t, size_t>> graph_sizes;
    graph_sizes.emplace_back(5, 7);
    graph_sizes.emplace_back(10, 12);
    graph_sizes.emplace_back(20, 35);
    for (auto& sizes : graph_sizes) {
        size_t num_nodes = sizes.first;
        size_t num_edges = sizes.second;
        for (size_t i = 0; i < num_reps; ++i) {
            BaseGraph graph = random_graph(num_nodes, num_edges, gen);
            make_path_cover(graph, gen);
            do_tests(graph);
        }
    }
    
    cerr << "passed all tests!" << endl;
}
