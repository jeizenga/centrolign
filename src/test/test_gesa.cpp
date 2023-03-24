#include <cstdio>
#include <cstdlib>
#include <cstdlib>

#include <iostream>
#include <vector>
#include <algorithm>
#include <random>

#include "centrolign/gesa.hpp"
#include "centrolign/range_min_query.hpp"

using namespace std;
using namespace centrolign;

class TestGESA : public GESA {
public:
    using GESA::lcp_array;
    using GESA::child_array;
    using GESA::suffix_links;
    using GESA::edges;
    using GESA::children;
    using GESA::root;
    using GESA::link;
    using GESA::construct_child_array;
    using GESA::construct_suffix_links;
    using GESA::child_array_is_up;
    using GESA::child_array_is_l_index;
};

size_t UP = 0, DOWN = 1, NEXT = 2;

pair<vector<size_t>, vector<size_t>> manual_child_array(std::vector<size_t>& lcp_array) {
    
    std::vector<size_t> up(lcp_array.size(), -1);
    std::vector<size_t> down(lcp_array.size(), -1);
    std::vector<size_t> next(lcp_array.size(), -1);
    
    // conditions translated directly from definition in paper
    for (size_t i = 0; i < lcp_array.size(); ++i) {
        // up
        for (size_t q = 0; q < i; ++q) {
            if (lcp_array[q] > lcp_array[i]) {
                size_t k = q + 1;
                while (k < i && lcp_array[k] >= lcp_array[q]) {
                    ++k;
                }
                if (k == i) {
                    up[i] = q;
                    break;
                }
            }
        }
        // down
        for (size_t q = lcp_array.size() - 1; q > i; --q) {
            if (lcp_array[q] > lcp_array[i]) {
                size_t k = q - 1;
                while (k > i && lcp_array[k] > lcp_array[q]) {
                    --k;
                }
                if (k == i) {
                    down[i] = q;
                    break;
                }
            }
        }
        // next
        for (size_t q = i + 1; q < lcp_array.size(); ++q) {
            if (lcp_array[q] == lcp_array[i]) {
                next[i] = q;
                break;
            }
            if (lcp_array[q] < lcp_array[i]) {
                break;
            }
        }
    }
    
    vector<size_t> child_array(lcp_array.size(), -1);
    vector<size_t> tab_source(lcp_array.size(), -1);
    for (size_t i = 0; i < child_array.size(); ++i) {
        if (next[i] != -1) {
            child_array[i] = next[i];
            tab_source[i] = NEXT;
        }
        else if (i + 1 < up.size() && up[i + 1] != -1) {
            child_array[i] = up[i + 1];
            tab_source[i] = UP;
        }
        else if (down[i] != -1) {
            child_array[i] = down[i];
            tab_source[i] = DOWN;
        }
    }
    
    assert(child_array.back() == -1);
    child_array.pop_back();
    tab_source.pop_back();
    return make_pair(child_array, tab_source);
}

vector<pair<GESANode, GESANode>> manual_suffix_link(TestGESA& gesa) {

    // FIXME: this strategy doesn't really work for leaf nodes, so i'm just
    // removing them for now
    
    vector<pair<GESANode, GESANode>> links;
    
    // get all the LCP intervals
    vector<GESANode> stack;
    stack.push_back(gesa.root());
    while (!stack.empty()) {
        auto node = stack.back();
        stack.pop_back();
        for (auto child : gesa.children(node)) {
            if (child.is_leaf()) {
                // FIXME: what's a more robust way to identify leaf links manually?
                continue;
            }
            links.emplace_back(child, GESANode());
            stack.push_back(child);
        }
    }
    
    for (auto& record : links) {
        auto& node = record.first;
        auto& link = record.second;
        
        if (gesa.edges[node.begin].empty()) {
            link = gesa.root();
            continue;
        }
        
        size_t l = 10000000;
        for (size_t i = node.begin + 1; i <= node.end; ++i) {
            l = min(l, gesa.lcp_array[i]);
        }
        size_t link_l = l - 1;
        size_t next = gesa.edges[node.begin].front();
        size_t begin = next;
        while (begin != 0 && gesa.lcp_array[begin] >= link_l) {
            --begin;
        }
        size_t end = next;
        while (end + 1 < gesa.lcp_array.size() && gesa.lcp_array[end + 1] >= link_l) {
            ++end;
        }
        link.begin = begin;
        link.end = end;
        for (auto n : gesa.edges[node.begin]) {
            assert(n >= begin && n <= end);
        }
    }
    
    return links;
}

void test_lcp_interval_tree(std::vector<size_t>& lcp_array,
                            std::vector<std::vector<uint64_t>>& edges) {
    
    TestGESA gesa;
    gesa.lcp_array = lcp_array;
    gesa.edges = edges;
    
    gesa.construct_child_array();
    gesa.construct_suffix_links();
    
    vector<size_t> expected_child_array, tab_source;
    tie(expected_child_array, tab_source) = manual_child_array(lcp_array);
    
    bool failed = gesa.child_array != expected_child_array;
    
    for (size_t i = 0; i < gesa.child_array.size(); ++i) {
        if (tab_source[i] == UP) {
            assert(gesa.child_array_is_up(i));
        }
        else if (tab_source[i] == NEXT) {
            assert(gesa.child_array_is_l_index(i));
        }
        else {
            assert(tab_source[i] == DOWN);
            assert(!gesa.child_array_is_up(i) && !gesa.child_array_is_l_index(i));
        }
    }
    
    
    
    auto suffix_links = manual_suffix_link(gesa);
    for (auto record : suffix_links) {
        auto node = record.first;
        auto expected_link = record.second;
        auto link = gesa.link(node);
        if (link != expected_link) {
            failed = true;
            cerr << "incorrect link, expected " << node.begin << "," << node.end << " -> " << expected_link.begin << "," << expected_link.end << ", got " << link.begin << "," << link.end << '\n';
        }
    }
    
    if (failed) {
        cerr << "failed LCP interval tree tests on array\n";
        for (auto i : lcp_array) {
            cerr << ' ' << i;
        }
        cerr << '\n';
        cerr << "expected child array\n";
        for (auto i : expected_child_array) {
            cerr << ' ' << i;
        }
        cerr << '\n';
        cerr << "actual child array\n";
        for (auto i : gesa.child_array) {
            cerr << ' ' << i;
        }
        cerr << '\n';
        exit(1);
    }
}

int main(int argc, char* argv[]) {
    
    {
        // from the ESA paper
        vector<size_t> lcp_array{0, 2, 1, 3, 1, 2, 0, 2, 0, 1, 0};
        vector<vector<uint64_t>> edges{{1}, {3}, {6}, {7}, {8}, {9}, {0}, {4}, {5}, {10}, {}};
        test_lcp_interval_tree(lcp_array, edges);
    }
    
    {
        vector<size_t> lcp_array{0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 2, 2};
        vector<vector<uint64_t>> edges{{9, 11}, {7}, {}, {}, {2}, {3}, {4}, {8, 10}, {5}, {4}, {5}, {6}};
        test_lcp_interval_tree(lcp_array, edges);
    }
    
    {
        vector<size_t> lcp_array{0, 2, 1, 1, 0, 0};
        vector<vector<uint64_t>> edges{{2}, {3}, {4}, {5}, {}, {}};
        test_lcp_interval_tree(lcp_array, edges);
    }
    
    {
        BaseGraph graph1, graph2;
        graph1.add_node('T');
        graph1.add_node('A');
        graph1.add_node('A');
        graph1.add_edge(0, 1);
        graph1.add_edge(0, 2);
        graph1.add_edge(1, 2);
        
        graph2.add_node('C');
        graph2.add_node('T');
        graph2.add_node('G');
        graph2.add_node('A');
        graph2.add_edge(0, 1);
        graph2.add_edge(0, 2);
        graph2.add_edge(1, 3);
        graph2.add_edge(2, 3);
        
        vector<const BaseGraph*> graphs{&graph1, &graph2};
        
        GESA gesa(graphs);
    }
    
    cerr << "passed all tests!" << endl;
}
