#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <utility>

#include "centrolign/modify_graph.hpp"
#include "centrolign/path_graph.hpp"
#include "centrolign/graph.hpp"
#include "centrolign/test_util.hpp"
#include "centrolign/determinize.hpp"
#include "centrolign/count_walks.hpp"

using namespace std;
using namespace centrolign;


BaseGraph to_base_graph(const PathGraph& path_graph, const BaseGraph& graph) {
    
    BaseGraph converted;
    for (uint64_t node_id = 0; node_id < path_graph.node_size(); ++node_id) {
        converted.add_node(graph.label(path_graph.from(node_id)));
    }
    
    for (uint64_t node_id = 0; node_id < path_graph.node_size(); ++node_id) {
        for (uint64_t next_id : path_graph.next(node_id)) {
            converted.add_edge(node_id, next_id);
        }
    }
    
    return converted;
}

bool is_prefix_sorted(const PathGraph& path_graph, const BaseGraph& graph) {
    
    std::vector<std::string> prefixes(path_graph.node_size());
    
    BaseGraph base_graph = to_base_graph(path_graph, graph);
    
    for (uint64_t node_id = 0; node_id < path_graph.node_size(); ++node_id) {
        auto paths = all_paths(base_graph, node_id);
        
        size_t lcp = 0;
        if (node_id > 0) {
            lcp = max(lcp, path_graph.lcp(node_id - 1));
        }
        if (node_id + 1 < path_graph.node_size()) {
            lcp = max(lcp, path_graph.lcp(node_id));
        }
        size_t prefix_len = lcp + 1;
        
        prefixes[node_id] = path_to_string(base_graph, paths.front()).substr(0, prefix_len);
        if (prefixes[node_id].size() != prefix_len) {
            cerr << "too short prefix, expected " << prefix_len << " got prefix " << prefixes[node_id] << '\n';
            return false;
        }
        
        for (size_t i = 1; i < paths.size(); ++i) {
            auto pref = path_to_string(base_graph, paths[i]).substr(0, prefix_len);
            if (pref != prefixes[node_id]) {
                cerr << "non matching prefixes for node " << node_id << '\n';
                cerr << prefixes[node_id] << '\n';
                cerr << pref << '\n';
                return false;
            }
        }
    }
    
    for (size_t i = 1; i < prefixes.size(); ++i) {
        if (prefixes[i] <= prefixes[i - 1]) {
            cerr << "prefix sorting failed for prefixes of " << (i - 1) << " and " << i << ":\n";
            cerr << prefixes[i - 1] << '\n';
            cerr << prefixes[i] << '\n';
            return false;
        }
    }
    
    return true;
}

void do_test(const BaseGraph& graph) {
    
    PathGraph path_graph(graph);
    while (!path_graph.is_prefix_sorted()) {
        path_graph = PathGraph(path_graph);
    }
    
    path_graph.finish(graph);
    
    // coarse way to measure path complexity
    uint64_t num_walks = count_walks(path_graph);
    
    // to prevent excessive run time
    if (num_walks < 3000) { // this still completes at 50k
        if (!is_prefix_sorted(path_graph, graph)) {
            cerr << "prefix sorted test failed on graph:\n";
            print_graph(graph, cerr);
            exit(1);
        }
    }
}


int main(int argc, char* argv[]) {
    
    random_device rd;
    default_random_engine gen(rd());
    
    {
        BaseGraph graph;
        for (auto c : string("$AAAGGAGCA#")) {
            graph.add_node(c);
        }
        graph.add_edge(1, 0);
        graph.add_edge(2, 1);
        graph.add_edge(3, 2);
        graph.add_edge(4, 1);
        graph.add_edge(5, 2);
        graph.add_edge(6, 3);
        graph.add_edge(6, 4);
        graph.add_edge(6, 5);
        graph.add_edge(7, 3);
        graph.add_edge(7, 4);
        graph.add_edge(7, 5);
        graph.add_edge(8, 2);
        graph.add_edge(8, 3);
        graph.add_edge(8, 4);
        graph.add_edge(8, 5);
        graph.add_edge(9, 6);
        graph.add_edge(9, 7);
        graph.add_edge(9, 8);
        graph.add_edge(10, 2);
        graph.add_edge(10, 3);
        graph.add_edge(10, 6);
        graph.add_edge(10, 9);
        
        do_test(graph);
    }
    
    {
        BaseGraph graph;
        for (auto c : string("$ATGAATACCAA#")) {
            graph.add_node(c);
        }
        graph.add_edge(1, 0);
        graph.add_edge(2, 0);
        graph.add_edge(3, 0);
        graph.add_edge(4, 1);
        graph.add_edge(5, 2);
        graph.add_edge(5, 3);
        graph.add_edge(6, 1);
        graph.add_edge(6, 2);
        graph.add_edge(6, 3);
        graph.add_edge(7, 4);
        graph.add_edge(7, 5);
        graph.add_edge(7, 6);
        graph.add_edge(8, 3);
        graph.add_edge(9, 4);
        graph.add_edge(9, 5);
        graph.add_edge(9, 6);
        graph.add_edge(10, 7);
        graph.add_edge(10, 8);
        graph.add_edge(10, 9);
        graph.add_edge(11, 10);
        graph.add_edge(12, 7);
        graph.add_edge(12, 10);
        graph.add_edge(12, 11);
        
        do_test(graph);
    }
    
    {
        BaseGraph graph;
        for (auto c : string("$AGAACAATAGTA#")) {
            graph.add_node(c);
        }
        graph.add_edge(1, 0);
        graph.add_edge(2, 0);
        graph.add_edge(3, 1);
        graph.add_edge(4, 2);
        graph.add_edge(5, 1);
        graph.add_edge(5, 2);
        graph.add_edge(6, 3);
        graph.add_edge(6, 4);
        graph.add_edge(7, 5);
        graph.add_edge(8, 1);
        graph.add_edge(8, 2);
        graph.add_edge(8, 3);
        graph.add_edge(8, 4);
        graph.add_edge(8, 5);
        graph.add_edge(9, 6);
        graph.add_edge(9, 7);
        graph.add_edge(9, 8);
        graph.add_edge(10, 1);
        graph.add_edge(10, 2);
        graph.add_edge(10, 3);
        graph.add_edge(10, 4);
        graph.add_edge(10, 5);
        graph.add_edge(10, 6);
        graph.add_edge(10, 7);
        graph.add_edge(10, 8);
        graph.add_edge(11, 6);
        graph.add_edge(11, 7);
        graph.add_edge(11, 8);
        graph.add_edge(12, 9);
        graph.add_edge(12, 10);
        graph.add_edge(12, 11);
        graph.add_edge(13, 3);
        graph.add_edge(13, 6);
        graph.add_edge(13, 9);
        graph.add_edge(13, 12);

        do_test(graph);
    }
    
    {
        BaseGraph graph;
        for (auto c : string("$CAAACAAAAAA#")) {
            graph.add_node(c);
        }
        graph.add_edge(1, 0);
        graph.add_edge(2, 0);
        graph.add_edge(3, 1);
        graph.add_edge(3, 2);
        graph.add_edge(4, 3);
        graph.add_edge(5, 4);
        graph.add_edge(6, 5);
        graph.add_edge(7, 6);
        graph.add_edge(8, 4);
        graph.add_edge(8, 7);
        graph.add_edge(9, 8);
        graph.add_edge(10, 9);
        graph.add_edge(11, 10);
        graph.add_edge(12, 10);
        graph.add_edge(12, 11);

        do_test(graph);
    }

    {
        BaseGraph graph;
        for (auto c : string("$CACAAGCCACTCAA#")) {
            graph.add_node(c);
        }
        graph.add_edge(1, 0);
        graph.add_edge(2, 1);
        graph.add_edge(3, 2);
        graph.add_edge(4, 2);
        graph.add_edge(5, 3);
        graph.add_edge(5, 4);
        graph.add_edge(6, 3);
        graph.add_edge(6, 4);
        graph.add_edge(7, 6);
        graph.add_edge(8, 3);
        graph.add_edge(8, 4);
        graph.add_edge(8, 5);
        graph.add_edge(9, 7);
        graph.add_edge(9, 8);
        graph.add_edge(10, 9);
        graph.add_edge(11, 7);
        graph.add_edge(11, 9);
        graph.add_edge(11, 10);
        graph.add_edge(12, 11);
        graph.add_edge(13, 11);
        graph.add_edge(14, 12);
        graph.add_edge(14, 13);
        graph.add_edge(15, 13);
        graph.add_edge(15, 14);

        do_test(graph);
    }

    for (int rep = 0; rep < 20; ++rep) {
        
        uniform_int_distribution<int> num_nodes_distr(10, 50);
        
        BaseGraph graph = random_challenge_graph(num_nodes_distr(gen), gen);
        auto tableau = add_sentinels(graph, '#', '$');
        
        BaseGraph determinized = determinize(graph);
        
//        cerr << "next graph:\n";
//        print_graph(determinized, cerr);
        
        do_test(determinized);
    }

    // an extended example that i worked out the answer to by hand:
    {
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

        path_graph.finish(graph);

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
    }

     // a debugging example designed to create distinct path nodes with the same
     // origin node that are lexicographically adjacent in the prefix sort
    {
        BaseGraph graph;
        // component 1
        graph.add_node('!');
        graph.add_node('G');
        graph.add_node('A');
        graph.add_node('C');
        graph.add_node('C');
        graph.add_node('A');
        graph.add_node('$');
        graph.add_edge(0, 1);
        graph.add_edge(1, 2);
        graph.add_edge(1, 4);
        graph.add_edge(2, 3);
        graph.add_edge(3, 6);
        graph.add_edge(4, 5);
        graph.add_edge(5, 6);
        // component 2
        graph.add_node('#');
        graph.add_node('G');
        graph.add_node('A');
        graph.add_node('A');
        graph.add_node('C');
        graph.add_node('T');
        graph.add_node('%');
        graph.add_edge(7, 8);
        graph.add_edge(8, 9);
        graph.add_edge(8, 11);
        graph.add_edge(9, 10);
        graph.add_edge(10, 13);
        graph.add_edge(11, 12);
        graph.add_edge(12, 13);

        PathGraph path_graph(graph);
        while (!path_graph.is_prefix_sorted()) {
            path_graph = PathGraph(path_graph);
        }

        path_graph.finish(graph);

        // the prefix should double to allow the 'G' nodes to have only 1
        // forward edge
        for (uint64_t node_id = 0; node_id < path_graph.node_size(); ++node_id) {
            if (graph.label(path_graph.from(node_id)) == 'G') {
                assert(path_graph.next_size(node_id) == 1);
            }
        }
    }

    {
        BaseGraph graph;
        // component 1
        graph.add_node('!');
        graph.add_node('G');
        graph.add_node('A');
        graph.add_node('A');
        graph.add_node('A');
        graph.add_node('A');
        graph.add_node('C');
        graph.add_node('A');
        graph.add_node('C');
        graph.add_node('A');
        graph.add_node('A');
        graph.add_node('A');
        graph.add_node('T');
        graph.add_node('%');
        graph.add_edge(0, 1);
        graph.add_edge(1, 2);
        graph.add_edge(1, 9);
        graph.add_edge(2, 3);
        graph.add_edge(2, 6);
        graph.add_edge(3, 4);
        graph.add_edge(4, 5);
        graph.add_edge(5, 13);
        graph.add_edge(6, 7);
        graph.add_edge(7, 8);
        graph.add_edge(8, 13);
        graph.add_edge(9, 10);
        graph.add_edge(10, 11);
        graph.add_edge(11, 12);
        graph.add_edge(12, 13);
        // component 2
        graph.add_node('#');
        graph.add_node('G');
        graph.add_node('A');
        graph.add_node('A');
        graph.add_node('A');
        graph.add_node('A');
        graph.add_node('C');
        graph.add_node('A');
        graph.add_node('C');
        graph.add_node('$');
        graph.add_node('&');
        graph.add_edge(14, 15);
        graph.add_edge(15, 16);
        graph.add_edge(16, 17);
        graph.add_edge(16, 20);
        graph.add_edge(17, 18);
        graph.add_edge(18, 19);
        graph.add_edge(19, 23);
        graph.add_edge(20, 21);
        graph.add_edge(21, 22);
        graph.add_edge(22, 24);

        PathGraph path_graph(graph);
        while (!path_graph.is_prefix_sorted()) {
            path_graph = PathGraph(path_graph);
        }

        path_graph.finish(graph);

        // the prefix should double to allow the 'G' nodes to have only 1
        // forward edge
        for (uint64_t node_id = 0; node_id < path_graph.node_size(); ++node_id) {
            if (graph.label(path_graph.from(node_id)) == 'G') {
                assert(path_graph.next_size(node_id) == 1);
            }
        }
    }
    
    cerr << "passed all tests!" << endl;
    return 0;
}
