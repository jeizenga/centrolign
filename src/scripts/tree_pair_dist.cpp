#include <unordered_map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <tuple>
#include <utility>
#include <sstream>

#include "centrolign/tree_distance_oracle.hpp"


using namespace std;
using namespace centrolign;

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
    
    TreeDistanceOracle oracle(tree);
    
    vector<uint64_t> leaves;
    for (uint64_t node_id = 0; node_id < tree.node_size(); ++node_id) {
        if (tree.is_leaf(node_id)) {
            leaves.push_back(node_id);
        }
    }
    
    cout << "sample1" << '\t' << "sample2" << '\t' << "distance" << '\n';
    for (size_t i = 0; i < leaves.size(); ++i) {
        for (size_t j = i + 1; j < leaves.size(); ++j) {
            
            cout << tree.label(leaves[i]) << '\t' << tree.label(leaves[j]) << '\t' << oracle.distance(leaves[i], leaves[j]) << '\n';
        }
    }
}
