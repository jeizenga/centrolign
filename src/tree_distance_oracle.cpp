#include "centrolign/tree_distance_oracle.hpp"


namespace centrolign {

using namespace std;

TreeDistanceOracle::TreeDistanceOracle(const Tree& tree) {
    
    static const bool debug = false;
    
    if (debug) {
        cerr << "tree:\n";
        for (uint64_t n = 0; n < tree.node_size(); ++n) {
            cerr << n << ":";
            for (auto c : tree.get_children(n)) {
                cerr << ' ' << c;
            }
            cerr << '\n';
        }
    }
    
    // do an Euler traversal
    
    // records of (node, next edge index)
    vector<pair<uint64_t, size_t>> stack;
    stack.emplace_back(tree.get_root(), 0);
    
    euler_nodes.reserve(2 * tree.node_size());
    euler_depths.reserve(2 * tree.node_size());
    
    while (!stack.empty()) {
        auto& top = stack.back();
        euler_nodes.push_back(top.first);
        euler_depths.push_back(stack.size());
        
        if (top.second == tree.get_children(top.first).size()) {
            stack.pop_back();
        }
        else {
            ++top.second;
            stack.emplace_back(tree.get_children(top.first).at(top.second - 1), 0);
        }
    }
    
    if (debug) {
        cerr << "euler traversal:\n";
        for (size_t i = 0; i < euler_nodes.size(); ++i) {
            cerr << i << '\t' << euler_nodes[i] << '\t' << euler_depths[i] << '\n';
        }
    }
    
    // process the traversal for LCA retrieval
    euler_rmq = move(RMQ<size_t>(euler_depths));
    
    // find an occurrence of each node on the Euler traversal
    position.resize(tree.node_size());
    for (size_t i = 0; i < euler_nodes.size(); ++i) {
        position[euler_nodes[i]] = i;
    }
    
    // compute the depth of each node
    {
        bool has_distances = false;
        for (uint64_t node_id = 0; node_id < tree.node_size(); ++node_id) {
            if (tree.distance(node_id) != numeric_limits<double>::max()) {
                has_distances = true;
                break;
            }
        }
        
        depths.resize(tree.node_size());
        
        // records of (node, depth)
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
    }
}

}
