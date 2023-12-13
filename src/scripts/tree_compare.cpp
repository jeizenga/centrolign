#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <tuple>
#include <utility>
#include <sstream>
#include <limits>

#include "centrolign/tree.hpp"
#include "centrolign/utility.hpp"

using namespace std;
using namespace centrolign;

vector<pair<double, vector<string>>> partition_table(const Tree& tree) {
    
    // compute min height by dynamic programming
    vector<double> height(tree.node_size(), numeric_limits<double>::max());
    // upward pass
    for (uint64_t node_id : tree.postorder()) {
        if (tree.is_leaf(node_id)) {
            height[node_id] = 0.0;
        }
        if (node_id != tree.get_root()) {
            auto parent_id = tree.get_parent(node_id);
            height[parent_id] = min(height[parent_id], height[node_id] + tree.distance(node_id));
        }
    }
    // downward pass
    // FIXME: this pass is only appropriate in an unrooted tree, but an ultrametric tree is robust either way
    for (uint64_t node_id : tree.preorder()) {
        if (node_id != tree.get_root()) {
            height[node_id] = min(height[node_id], height[tree.get_parent(node_id)] + tree.distance(node_id));
        }
    }
    
    // get the labels of all of the labels
    vector<string> labels;
    for (uint64_t node_id = 0; node_id < tree.node_size(); ++node_id) {
        if (tree.is_leaf(node_id)) {
            labels.push_back(tree.label(node_id));
        }
    }
    sort(labels.begin(), labels.end());
    
    // assumes sorted label sets
    auto flip_label_set = [&](const vector<string>& label_set) -> vector<string> {
        
        vector<string> flipped;
        size_t i = 0, j = 0;
        while (i < labels.size()) {
            if (j == label_set.size() || labels[i] < label_set[j]) {
                flipped.push_back(labels[i]);
                ++i;
            }
            else {
                ++i;
                ++j;
            }
            // note: since labels is a superset, we don't need to handle case of
            // label_set[j] < labels[i]
        }
        
        return flipped;
    };
    
    // destructive to the inputs, assumes sorted and disjoint
    auto merge = [](vector<string>& label_set1, vector<string>& label_set2) {
        
        vector<string> merged;
        merged.reserve(label_set1.size() + label_set2.size());
        
        size_t i = 0, j = 0;
        while (i < label_set1.size() || j < label_set2.size()) {
            if (i == label_set1.size() || (j != label_set2.size() && label_set2[j] < label_set1[i])) {
                merged.emplace_back(move(label_set2[j]));
                ++j;
            }
            else {
                merged.emplace_back(move(label_set1[i]));
                ++i;
            }
        }
        
        return merged;
    };
    
    vector<pair<double, vector<string>>> return_val;
    vector<vector<string>> label_sets(tree.node_size());
    
    for (uint64_t node_id : tree.postorder()) {
        if (tree.is_leaf(node_id)) {
            // base condition
            label_sets[node_id].push_back(tree.label(node_id));
        }
        else if (tree.get_children(node_id).size() + (node_id == tree.get_root() ? 0 : 1) > 2) {
            // this node could correspond to a non-trivial, non-redundant bipartition
            auto children = tree.get_children(node_id);
            assert(children.size() >= 2);
            label_sets[node_id] = merge(label_sets[children[0]], label_sets[children[1]]);
            for (size_t i = 2; i < children.size(); ++i) {
                auto tmp = merge(label_sets[node_id], label_sets[children[i]]);
                label_sets[node_id] = move(tmp);
            }
            
            // normalize to always include one arbitrary
            if (label_sets[node_id].front() == labels.front()) {
                return_val.emplace_back(height[node_id], label_sets[node_id]);
            }
            else {
                return_val.emplace_back(height[node_id], flip_label_set(label_sets[node_id]));
            }
        }
    }
    
    // we have to deduplicate all nodes along each nonbranching path (this shows up near the root)
    // FIXME: if i identified these regions beforehand, i wouldn't need the logn runtime penalty here
    sort(return_val.begin(), return_val.end(),
         [](const pair<double, vector<string>>& a, const pair<double, vector<string>>& b) {
        return a.second < b.second || (a.second == b.second && a.first < b.first);
    });
    auto unique_end = unique(return_val.begin(), return_val.end(),
                             [](const pair<double, vector<string>>& a, const pair<double, vector<string>>& b) {
        return a.second == b.second;
    });
    return_val.resize(unique_end - return_val.begin());
    
    return return_val;
}

namespace std {
template<typename T>
struct hash<std::vector<T>> {
    size_t operator()(const vector<T>& vec) const {
        size_t seed = 0;
        for (const T& val : vec) {
            hash_combine(seed, std::hash<T>()(val));
        }
        return seed;
    }
};
}

int main(int argc, char* argv[]) {
    
    static const bool debug = false;
    
    if (argc != 3) {
        cerr << "usage:\n";
        cerr << "tree_compare truth_tree.nwk compare_tree.nwk > subtree_correctness.tsv\n";
        return 1;
    }
    
    ifstream truth_in(argv[1]);
    if (!truth_in) {
        cerr << "error: could not open tree file " << argv[1] << '\n';
        return 1;
    }
    
    ifstream compare_in(argv[2]);
    if (!compare_in) {
        cerr << "error: could not open tree file " << argv[2] << '\n';
        return 1;
    }
    
    // read it from the file
    stringstream truth_sstrm;
    truth_sstrm << truth_in.rdbuf();
    string truth_newick = truth_sstrm.str();
    
    stringstream compare_sstrm;
    compare_sstrm << compare_in.rdbuf();
    string compare_newick = compare_sstrm.str();
    
    Tree truth(truth_newick);
    Tree compare(compare_newick);
    
    // make sure the trees have the same leaves
    size_t num_leaves = 0;
    for (uint64_t node_id = 0; node_id < truth.node_size(); ++node_id) {
        if (truth.is_leaf(node_id)) {
            assert(!truth.label(node_id).empty());
            assert(compare.has_label(truth.label(node_id)));
            assert(compare.is_leaf(compare.get_id(truth.label(node_id))));
            ++num_leaves;
        }
    }
    for (uint64_t node_id = 0; node_id < compare.node_size(); ++node_id) {
        if (compare.is_leaf(node_id)) {
            assert(!compare.label(node_id).empty());
            assert(truth.has_label(compare.label(node_id)));
            assert(truth.is_leaf(truth.get_id(compare.label(node_id))));
        }
    }
    
    
    auto truth_table = partition_table(truth);
    auto compare_table = partition_table(compare);
    
    if (debug) {
        for (auto p : {&truth_table, &compare_table}) {
            auto& tab = *p;
            cerr << (p == &truth_table ? "truth" : "comparison") << " table:\n";
            for (auto& r : tab) {
                cerr << r.first << '\t';
                for (size_t i = 0; i < r.second.size(); ++i) {
                    if (i) {
                        cerr << ',';
                    }
                    cerr << r.second[i];
                }
                cerr << '\n';
            }
        }
        cerr << '\n';
    }
    
    unordered_set<vector<string>> compare_rows;
    for (auto& partition : compare_table) {
        compare_rows.insert(move(partition.second));
    }
    
    for (auto& row : truth_table) {
        cout << row.first << '\t' << min<size_t>(row.second.size(), num_leaves - row.second.size()) << '\t' << compare_rows.count(row.second) << '\n';
    }
}
