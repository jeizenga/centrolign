#include <cstdio>
#include <cstdlib>
#include <cstdlib>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <utility>
#include <cassert>

#include "centrolign/utility.hpp"
#include "centrolign/random_graph.hpp"
#include "centrolign/graph.hpp"
#include "centrolign/superbubbles.hpp"
#include "centrolign/topological_order.hpp"
#include "centrolign/modify_graph.hpp"

using namespace std;
using namespace centrolign;

class TestSuperbubbleTree : public SuperbubbleTree {
public:
    TestSuperbubbleTree(const BaseGraph& graph) : SuperbubbleTree(graph) {}
    using SuperbubbleTree::find_superbubbles;
};

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
            
            if (reachable && reachable_forward == reachable_backward) {
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
    auto got = TestSuperbubbleTree::find_superbubbles(graph);
    
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
        for (int i = 0; i < 12; ++i) {
            graph.add_node('A');
        }
        graph.add_edge(0, 1);
        graph.add_edge(0, 7);
        graph.add_edge(0, 8);
        graph.add_edge(1, 2);
        graph.add_edge(1, 3);
        graph.add_edge(2, 4);
        graph.add_edge(3, 4);
        graph.add_edge(4, 5);
        graph.add_edge(4, 6);
        graph.add_edge(5, 6);
        graph.add_edge(6, 11);
        graph.add_edge(7, 8);
        graph.add_edge(8, 9);
        graph.add_edge(8, 10);
        graph.add_edge(9, 10);
        graph.add_edge(10, 11);
        
        SuperbubbleTree tree(graph);
        
        assert(tree.superbubble_size() == 4);
        assert(tree.chain_size() == 3);
        
        for (int n : {2, 3, 5, 6, 7, 9, 10, 11}) {
            assert(tree.superbubble_beginning_at(n) == -1);
        }
        for (int n : {0, 1, 2, 3, 5, 7, 8, 9}) {
            assert(tree.superbubble_ending_at(n) == -1);
        }
        
        uint64_t bub1 = tree.superbubble_beginning_at(0);
        uint64_t bub2 = tree.superbubble_beginning_at(1);
        uint64_t bub3 = tree.superbubble_beginning_at(4);
        uint64_t bub4 = tree.superbubble_beginning_at(8);
        assert(bub1 != -1);
        assert(bub2 != -1);
        assert(bub3 != -1);
        assert(bub4 != -1);
        assert(bub1 != bub2);
        assert(bub1 != bub3);
        assert(bub1 != bub4);
        assert(bub2 != bub3);
        assert(bub2 != bub4);
        assert(bub3 != bub4);
        assert(bub1 == tree.superbubble_ending_at(11));
        assert(bub2 == tree.superbubble_ending_at(4));
        assert(bub3 == tree.superbubble_ending_at(6));
        assert(bub4 == tree.superbubble_ending_at(10));
        pair<uint64_t, uint64_t> boundary1(0, 11);
        pair<uint64_t, uint64_t> boundary2(1, 4);
        pair<uint64_t, uint64_t> boundary3(4, 6);
        pair<uint64_t, uint64_t> boundary4(8, 10);
        assert(tree.superbubble_boundaries(bub1) == boundary1);
        assert(tree.superbubble_boundaries(bub2) == boundary2);
        assert(tree.superbubble_boundaries(bub3) == boundary3);
        assert(tree.superbubble_boundaries(bub4) == boundary4);
        
        uint64_t chain1 = tree.chain_containing(bub1);
        uint64_t chain2 = tree.chain_containing(bub2);
        uint64_t chain3 = tree.chain_containing(bub4);
        assert(chain1 != -1);
        assert(chain2 != -1);
        assert(chain3 != -1);
        assert(chain1 != chain2);
        assert(chain1 != chain3);
        assert(chain2 != chain3);
        assert(chain2 == tree.chain_containing(bub3));
        
        assert(tree.chains_inside(bub1).size() == 2);
        assert(find(tree.chains_inside(bub1).begin(),
                    tree.chains_inside(bub1).end(),
                    chain2) != tree.chains_inside(bub1).end());
        assert(find(tree.chains_inside(bub1).begin(),
                    tree.chains_inside(bub1).end(),
                    chain3) != tree.chains_inside(bub1).end());
        assert(tree.chains_inside(bub2).empty());
        assert(tree.chains_inside(bub3).empty());
        assert(tree.chains_inside(bub4).empty());
        
        vector<uint64_t> chain_bubs1{bub1};
        vector<uint64_t> chain_bubs2{bub2, bub3};
        vector<uint64_t> chain_bubs3{bub4};
        assert(tree.superbubbles_inside(chain1) == chain_bubs1);
        assert(tree.superbubbles_inside(chain2) == chain_bubs2);
        assert(tree.superbubbles_inside(chain3) == chain_bubs3);
        
        assert(tree.superbubble_containing(chain1) == -1);
        assert(tree.superbubble_containing(chain2) == bub1);
        assert(tree.superbubble_containing(chain3) == bub1);
        
        NetGraph ng1(graph, tree, bub1);
        NetGraph ng2(graph, tree, bub2);
        NetGraph ng3(graph, tree, bub3);
        NetGraph ng4(graph, tree, bub4);
        
        assert(ng1.node_size() == 5);
        assert(ng2.node_size() == 4);
        assert(ng3.node_size() == 3);
        assert(ng4.node_size() == 3);
        
        // only doing the hardest one, cuz this is way too tedious
        {
            uint64_t n1, n2, n3, n4, n5;
            n1 = n2 = n3 = n4 = n5 = -1;
            uint64_t num_edges = 0;
            for (uint64_t i = 0; i < ng1.node_size(); ++i) {
                num_edges += ng1.next_size(i);
                if (ng1.label(i) == pair<uint64_t, bool>(0, false)) {
                    n1 = i;
                }
                else if (ng1.label(i) == pair<uint64_t, bool>(chain2, true)) {
                    n2 = i;
                }
                else if (ng1.label(i) == pair<uint64_t, bool>(7, false)) {
                    n3 = i;
                }
                else if (ng1.label(i) == pair<uint64_t, bool>(chain3, true)) {
                    n4 = i;
                }
                else if (ng1.label(i) == pair<uint64_t, bool>(11, false)) {
                    n5 = i;
                }
            }
            assert(n1 != -1);
            assert(n2 != -1);
            assert(n3 != -1);
            assert(n4 != -1);
            assert(n5 != -1);
            assert(num_edges == 6);
        }
    }
    
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
