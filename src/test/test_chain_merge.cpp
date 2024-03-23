#include <cstdio>
#include <cstdlib>
#include <cstdlib>
#include <unordered_set>
#include <vector>
#include <list>
#include <string>

#include "centrolign/chain_merge.hpp"
#include "centrolign/path_merge.hpp"
#include "centrolign/test_util.hpp"
#include "centrolign/modify_graph.hpp"
#include "centrolign/utility.hpp"
#include "centrolign/step_index.hpp"

using namespace std;
using namespace centrolign;

template<class XMerge>
void do_tests(BaseGraph& graph, SentinelTableau* tableau = nullptr) {
    
    XMerge chain_merge = tableau ? XMerge(graph, *tableau) : XMerge(graph);
    
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
    
    for (uint64_t n = 0; n < graph.node_size(); ++n) {
        size_t num = 0;
        for (auto c : chain_merge.chains_on(n)) {
            ++num;
            auto i = chain_merge.index_on(n, c);
            if (chain_merge.node_at(c, i) != n) {
                cerr << "inconsistent bookkeepping for node " << n << ", chain " << c << ", index " << i << " on graph\n";
                print_graph(graph, cerr);
                exit(1);
            }
            if (chain_merge.predecessor_indexes(n)[c] != -1) {
                assert(chain_merge.predecessor_indexes(n)[c] == i - 1);
            }
            else {
                assert(i == 0);
            }
        }
        assert(num != 0);
    }
    
    // test this here too, just for convenience
    StepIndex step_index(graph);
    for (uint64_t n = 0; n < graph.node_size(); ++n) {
        for (const auto& p : step_index.path_steps(n)) {
            if (graph.path(p.first)[p.second] != n) {
                cerr << "step index tests failed\n";
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
        graph.add_edge(0, 4);
        graph.add_edge(0, 1);
        graph.add_edge(0, 2);
        graph.add_edge(1, 4);
        graph.add_edge(1, 3);
        graph.add_edge(2, 4);
        graph.add_edge(3, 4);
        
        graph.add_path("1");
        graph.add_path("2");
        graph.extend_path(0, 0);
        graph.extend_path(0, 1);
        graph.extend_path(0, 3);
        graph.extend_path(0, 4);
        graph.extend_path(1, 0);
        graph.extend_path(1, 2);
        graph.extend_path(1, 4);
        
        do_tests<ChainMerge>(graph);
        do_tests<PathMerge<>>(graph);
        auto tableau = add_sentinels(graph, '^', '$');
        do_tests<ChainMerge>(graph, &tableau);
        do_tests<PathMerge<>>(graph, &tableau);
    }
    
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
        
        do_tests<ChainMerge>(graph);
        do_tests<PathMerge<>>(graph);
        auto tableau = add_sentinels(graph, '^', '$');
        do_tests<ChainMerge>(graph, &tableau);
        do_tests<PathMerge<>>(graph, &tableau);
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
            add_random_path_cover(graph, gen);
            do_tests<ChainMerge>(graph);
            do_tests<PathMerge<>>(graph);
            auto tableau = add_sentinels(graph, '^', '$');
            do_tests<ChainMerge>(graph, &tableau);
            do_tests<PathMerge<>>(graph, &tableau);
        }
    }
    
    cerr << "passed all tests!" << endl;
}
