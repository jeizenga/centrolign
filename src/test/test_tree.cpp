#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <unordered_set>
#include <vector>
#include <list>
#include <string>
#include <random>
#include <algorithm>

#include "centrolign/tree.hpp"

using namespace std;
using namespace centrolign;


class TestTree : public Tree {
public:
    TestTree(const string& str) : Tree(str) {}
    using typename Tree::Node;
    using Tree::nodes;
    using Tree::root;
};

void validate_tree(TestTree& tree) {
    if (tree.nodes.size() != 0) {
        assert(tree.root != -1);
    }
    for (size_t i = 0; i < tree.nodes.size(); ++i) {
        auto& node = tree.nodes[i];
        if (node.parent == -1) {
            assert(i == tree.root);
        }
        else {
            auto& parent = tree.nodes[node.parent];
            assert(find(parent.children.begin(), parent.children.end(), i) != parent.children.end());
        }
    }
}

void test_postorder(TestTree& tree) {
    auto order = tree.postorder();
    assert(order.size() == tree.nodes.size());
    for (size_t i = 0; i < order.size(); ++i) {
        for (auto j : tree.nodes[order[i]].children) {
            auto idx = find(order.begin(), order.end(), j) - order.begin();
            assert(idx < i);
        }
    }
}

void test_binarize(TestTree& tree) {
    tree.binarize();
    for (auto& node : tree.nodes) {
        assert(node.children.size() <= 2);
    }
}

void do_tests(TestTree& tree) {
    validate_tree(tree);
    test_postorder(tree);
    test_binarize(tree);
    validate_tree(tree);
    test_postorder(tree);
}

int main(int argc, char* argv[]) {
     
    random_device rd;
    default_random_engine gen(rd());
    
    // TODO: only testing that these don't crash and produce valid results, but the results
    // have been verified by eye using the debug output
    
    {
        TestTree tree("((A,B),(C,D));");
        do_tests(tree);
    }

    {
        TestTree tree("((A),(B,C));");
        do_tests(tree);
    }

    vector<string> wiki_examples{
        "(,,(,));",
        "(A,B,(C,D));",
        "(A,B,(C,D)E)F;",
        "(:0.1,:0.2,(:0.3,:0.4):0.5);",
        "(:0.1,:0.2,(:0.3,:0.4):0.5):0.0;",
        "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);",
        "(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;",
        "((B:0.2,(C:0.3,D:0.4)E:0.5)F:0.1)A;"
    };
    
    for (auto& newick_str : wiki_examples) {
        TestTree tree(newick_str);
        do_tests(tree);
    }
    
    cerr << "passed all tests!" << endl;
}
