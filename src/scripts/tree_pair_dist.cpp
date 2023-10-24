#include <unordered_map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <tuple>
#include <utility>
#include <sstream>

#include "centrolign/tree.hpp"
#include "centrolign/range_min_query.hpp"


using namespace std;
using namespace centrolign;

void euler_traversal(const Tree& tree, vector<uint64_t>& nodes_out,
                     vector<size_t>& depths_out) {
    
    
    // records of (node, next edge index)
    vector<pair<uint64_t, size_t>> stack;
    stack.emplace_back(tree.get_root(), 0);
    
    while (!stack.empty()) {
        auto& top = stack.back();
        nodes_out.push_back(top.first);
        depths_out.push_back(stack.size());
        
        if (top.second == tree.get_children(top.first).size()) {
            stack.pop_back();
        }
        else {
            ++top.second;
            stack.emplace_back(tree.get_children(top.first).at(top.second - 1), 0);
        }
    }
    
}

vector<double> get_depths(const Tree& tree) {
    
    bool has_distances = false;
    for (uint64_t node_id = 0; node_id < tree.node_size(); ++node_id) {
        if (tree.distance(node_id) != numeric_limits<double>::max()) {
            has_distances = true;
            break;
        }
    }
    
    vector<double> depths(tree.node_size());
    
    vector<pair<uint64_t, double>> stack;
    stack.emplace_back(tree.get_root(), 0.0);
    
    while (!stack.empty()) {
        auto top = stack.back();
        stack.pop_back();
        depths[top.first] = top.second;
        for (auto child_id : tree.get_children(top.first)) {
            stack.emplace_back(child_id, top.second + (has_distances ? tree.distance(child_id) : 1.0));
        }
    }
    
    return depths;
}

int main(int argc, char* argv[]) {
    
    if (argc != 2) {
        cerr << "usage:\n";
        cerr << "tree_pair_dist tree.nwk > pair_dists.tsv\n";
        return 1;
    }
    
    ifstream tree_in(argv[1]);
    if (!tree_in) {
        cerr << "error: could not open tree file " << argv[1] << '\n';
        return 1;
    }
    
    // read it from the file
    stringstream sstrm;
    sstrm << tree_in.rdbuf();
    string newick_string = sstrm.str();
    
    Tree tree(newick_string);
    
    // do the reduction of LCA to RMQ
    vector<uint64_t> euler_nodes;
    vector<size_t> euler_depths;
    euler_traversal(tree, euler_nodes, euler_depths);
    RMQ<size_t> euler_rmq(euler_depths);
    
    vector<uint64_t> leaves;
    for (uint64_t node_id = 0; node_id < tree.node_size(); ++node_id) {
        if (tree.is_leaf(node_id)) {
            leaves.push_back(node_id);
        }
    }
    
    // map the leaves to their position in the euler array
    unordered_map<uint64_t, size_t> leaf_euler_pos;
    for (auto leaf_id : leaves) {
        leaf_euler_pos[leaf_id] = -1;
    }
    for (size_t i = 0; i < euler_nodes.size(); ++i) {
        auto it = leaf_euler_pos.find(euler_nodes[i]);
        if (it != leaf_euler_pos.end()) {
            it->second = i;
        }
    }
    
    // compute depth from root for each node
    vector<double> depths = get_depths(tree);
    
    cout << "sample1" << '\t' << "sample2" << '\t' << "distance" << '\n';
    for (size_t i = 0; i < leaves.size(); ++i) {
        for (size_t j = i + 1; j < leaves.size(); ++j) {
            size_t lo = leaf_euler_pos.at(leaves[i]);
            size_t hi = leaf_euler_pos.at(leaves[j]);
            if (hi < lo) {
                swap(lo, hi);
            }
            uint64_t lca_id = euler_nodes[euler_rmq.range_arg_min(lo, hi + 1)];
            
            double distance = depths[leaves[i]] + depths[leaves[j]] - 2 * depths[lca_id];
            
            cout << tree.label(leaves[i]) << '\t' << tree.label(leaves[j]) << '\t' << distance << '\n';
        }
    }
}
