#include <cstdio>
#include <cstdlib>
#include <cstdlib>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <cassert>

#include "centrolign/utility.hpp"
#include "centrolign/random_graph.hpp"
#include "centrolign/graph.hpp"
#include "centrolign/superbubbles.hpp"
#include "centrolign/topological_order.hpp"
#include "centrolign/modify_graph.hpp"

using namespace std;
using namespace centrolign;

vector<pair<uint64_t, uint64_t>> brute_force_superbubbles(BaseGraph& graph) {
    
    vector<tuple<uint64_t, uint64_t, set<uint64_t>>> candidates;
    
    auto order = topological_order(graph);
    
    for (size_t i = 0; i < order.size(); ++i) {
        for (size_t j = i + 1; j < order.size(); ++j) {
            
            uint64_t id_from = order[i];
            uint64_t id_to = order[j];
            
            set<uint64_t> reachable_forward;
            set<uint64_t> reachable_backward;
            
            // find reachable without crossing end from start
            vector<bool> traversed(graph.node_size(), false);
            vector<uint64_t> stack(1, id_from);
            traversed[id_from] = true;
            bool reachable = false;
            while (!stack.empty()) {
                auto here = stack.back();
                stack.pop_back();
                for (auto nid : graph.next(here)) {
                    if (nid == id_to) {
                        reachable = true;
                    }
                    else if (!traversed[nid]) {
                        traversed[nid] = true;
                        stack.push_back(nid);
                        reachable_forward.insert(nid);
                    }
                }
            }
            
            // find reachable without crossing start from end
            
            traversed.clear();
            traversed.resize(graph.node_size(), false);
            stack.emplace_back(id_to);
            traversed[id_to] = true;
            
            while (!stack.empty()) {
                auto here = stack.back();
                stack.pop_back();
                for (auto nid : graph.previous(here)) {
                    if (nid != id_from && !traversed[nid]) {
                        traversed[nid] = true;
                        stack.push_back(nid);
                        reachable_backward.insert(nid);
                    }
                }
            }
            
//            cerr << "from " << id_from << ", to " << id_to << ", reachable? " << reachable << '\n';
//            cerr << "reachable fwd:\n";
//            for (auto n : reachable_forward) {
//                cerr <<  ' ' << n;
//            }
//            cerr << '\n';
//            cerr << "reachable bwd:\n";
//            for (auto n : reachable_backward) {
//                cerr <<  ' ' << n;
//            }
//            cerr << '\n';
            
            if (reachable && reachable_forward == reachable_backward) {
//                cerr << "recording candidate\n";
                candidates.emplace_back(id_from, id_to, move(reachable_forward));
            }
            
        }
    }
    
    // verify minimality
    vector<pair<uint64_t, uint64_t>> superbubbles;
    for (size_t i = 0; i < candidates.size(); ++i) {
        auto s = get<0>(candidates[i]);
        auto t = get<1>(candidates[i]);
        auto& contents = get<2>(candidates[i]);
        bool keep = true;
        for (size_t j = 0; j < candidates.size(); ++j) {
            if (j == i) {
                continue;
            }
            
            if (s == get<0>(candidates[j]) && contents.count(get<1>(candidates[j]))) {
//                cerr << "(" << s << ", " << t << ") split by (" << s << ", " << get<1>(candidates[j]) << ")\n";
                keep = false;
                break;
            }
        }
        
        if (keep) {
            superbubbles.emplace_back(s, t);
        }
    }
    
    return superbubbles;
}


void do_test(BaseGraph& graph) {
    
    auto expected = brute_force_superbubbles(graph);
    auto got = find_superbubbles(graph);
    
    sort(expected.begin(), expected.end());
    sort(got.begin(), got.end());
    
    if (got != expected) {
        cerr << "failed test on graph:\n";
        print_graph(graph, cerr);
        cerr << "expected:\n";
        for (auto p : expected) {
            cerr << '\t' << p.first << ", " << p.second << '\n';
        }
        cerr << "got:\n";
        for (auto p : got) {
            cerr << '\t' << p.first << ", " << p.second << '\n';
        }
        exit(1);
    }
}

int main(int argc, char* argv[]) {
     
    random_device rd;
    default_random_engine gen(rd());

    {
        BaseGraph graph;
        for (int i = 0; i < 10; ++i) {
            graph.add_node('A');
        }
        graph.add_edge(0, 7);
        graph.add_edge(0, 2);
        graph.add_edge(0, 1);
        graph.add_edge(0, 4);
        graph.add_edge(1, 3);
        graph.add_edge(1, 6);
        graph.add_edge(1, 7);
        graph.add_edge(2, 5);
        graph.add_edge(2, 6);
        graph.add_edge(2, 4);
        graph.add_edge(3, 9);
        graph.add_edge(4, 8);
        graph.add_edge(4, 6);
        graph.add_edge(5, 9);
        graph.add_edge(7, 8);
        
        // make it single-source, single-sink
        add_sentinels(graph, '^', '$');
        
        do_test(graph);
    }
    
    {
        BaseGraph graph;
        for (int i = 0; i < 10; ++i) {
            graph.add_node('A');
        }
        graph.add_edge(0, 1);
        graph.add_edge(0, 3);
        graph.add_edge(0, 7);
        graph.add_edge(0, 9);
        graph.add_edge(1, 4);
        graph.add_edge(2, 3);
        graph.add_edge(2, 8);
        graph.add_edge(2, 9);
        graph.add_edge(3, 7);
        graph.add_edge(3, 8);
        graph.add_edge(3, 9);
        graph.add_edge(4, 6);
        graph.add_edge(4, 8);
        graph.add_edge(5, 8);
        graph.add_edge(7, 8);
        
        // make it single-source, single-sink
        add_sentinels(graph, '^', '$');
        
        do_test(graph);
    }
    
    {
        BaseGraph graph;
        for (int i = 0; i < 10; ++i) {
            graph.add_node('A');
        }
        graph.add_edge(0, 3);
        graph.add_edge(0, 7);
        graph.add_edge(0, 9);
        graph.add_edge(0, 1);
        graph.add_edge(1, 4);
        graph.add_edge(2, 8);
        graph.add_edge(2, 9);
        graph.add_edge(2, 3);
        graph.add_edge(3, 9);
        graph.add_edge(3, 8);
        graph.add_edge(3, 7);
        graph.add_edge(4, 8);
        graph.add_edge(4, 6);
        graph.add_edge(5, 8);
        graph.add_edge(7, 8);
        
        // make it single-source, single-sink
        add_sentinels(graph, '^', '$');
        
        do_test(graph);
    }
    
    {
        BaseGraph graph;
        for (int i = 0; i < 4; ++i) {
            graph.add_node('A');
        }
        graph.add_edge(0, 1);
        graph.add_edge(0, 2);
        graph.add_edge(1, 3);
        graph.add_edge(2, 3);

        do_test(graph);
    }

    {
        BaseGraph graph;
        for (int i = 0; i < 7; ++i) {
            graph.add_node('A');
        }
        graph.add_edge(0, 1);
        graph.add_edge(0, 2);
        graph.add_edge(1, 6);
        graph.add_edge(2, 3);
        graph.add_edge(2, 4);
        graph.add_edge(3, 5);
        graph.add_edge(4, 5);
        graph.add_edge(5, 6);

        do_test(graph);
    }
    
    {
        BaseGraph graph;
        for (int i = 0; i < 8; ++i) {
            graph.add_node('A');
        }
        graph.add_edge(0, 1);
        graph.add_edge(0, 2);
        graph.add_edge(1, 3);
        graph.add_edge(2, 3);
        graph.add_edge(3, 4);
        graph.add_edge(4, 5);
        graph.add_edge(4, 6);
        graph.add_edge(5, 7);
        graph.add_edge(6, 7);
        
        do_test(graph);
    }
    
    size_t num_reps = 50;
    vector<pair<size_t, size_t>> graph_sizes{
        {10, 12},
        {20, 30}
    };
    for (auto& sizes : graph_sizes) {
        size_t num_nodes = sizes.first;
        size_t num_edges = sizes.second;
        for (size_t i = 0; i < num_reps; ++i) {
            BaseGraph graph = random_graph(num_nodes, num_edges, gen);
            // make it single-source, single-sink
            add_sentinels(graph, '^', '$');
            do_test(graph);
        }
    }
    
    cerr << "passed all tests!" << endl;
}
