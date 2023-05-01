#include "centrolign/core.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <cassert>

#include "centrolign/utility.hpp"

namespace centrolign {

using namespace std;


Core::Core(const std::string& fasta_file, const std::string& tree_file) {
    
    ifstream fasta_fstream, tree_fstream;
    
    istream* fasta_stream = get_input(fasta_file, fasta_fstream);
    istream* tree_stream = get_input(tree_file, tree_fstream);
    
    auto sequences = parse_fasta(*fasta_stream);
    
    unordered_map<string, size_t> name_to_idx;
    for (size_t i = 0; i < sequences.size(); ++i) {
        auto& name = sequences[i].first;
        if (name_to_idx.count(name)) {
            throw runtime_error("FASTA contains duplicate name " + name);
        }
        name_to_idx[name] = i;
    }
    
    // read the newick stream into a tree
    stringstream sstrm;
    sstrm << tree_stream->rdbuf();
    tree = Tree(sstrm.str());
    
    // check the match between the fasta and the tree
    {
        std::vector<uint64_t> sequence_leaf_ids;
        for (size_t i = 0; i < sequences.size(); ++i) {
            
            const string& name = sequences[i].first;
            if (!tree.has_label(name)) {
                throw runtime_error("Guide tree does not include sequence " + name);
            }
            auto node_id = tree.get_id(name);
            if (!tree.is_leaf(node_id)) {
                throw runtime_error("Sequence " + name + " is not a leaf in the guide tree");
            }
            sequence_leaf_ids.push_back(node_id);
        }
        
        // get rid of samples we don't have the tree for
        tree.prune(sequence_leaf_ids);
    }
    
    // convert into a binary tree
    tree.binarize();
    
    // get rid of non-branching paths
    tree.compact();
    
    subproblems.resize(tree.node_size());
    
    for (uint64_t node_id = 0; node_id < tree.node_size(); ++node_id) {
        if (tree.is_leaf(node_id)) {
            const auto& name = tree.label(node_id);
            const auto& sequence = sequences[name_to_idx[name]].second;
            
            auto& subproblem = subproblems[node_id];
            
            subproblem.graph = make_base_graph(name, sequence);
            subproblem.tableau = add_sentinels(subproblem.graph, 5, 6);
            subproblem.name = name;
        }
    }
}

void Core::execute() {
    
    for (auto node_id : tree.postorder()) {
        if (tree.is_leaf(node_id)) {
            continue;
        }
        
        const auto& children = tree.get_children(node_id);
        assert(children.size() == 2);
        
        auto& subproblem1 = subproblems[children.front()];
        auto& subproblem2 = subproblems[children.back()];
        
        // compute chain merge data structures to use in the algorithm
        ChainMerge chain_merge1(subproblem1.graph, subproblem1.tableau);
        ChainMerge chain_merge2(subproblem2.graph, subproblem2.tableau);
        
        // give them unique sentinels for anchoring
        reassign_sentinels(subproblem1.graph, subproblem1.tableau, 5, 6);
        reassign_sentinels(subproblem2.graph, subproblem2.tableau, 7, 8);
        
        // compute the alignment
        auto anchors = anchorer.anchor_chain(subproblem1.graph, subproblem2.graph,
                                             chain_merge1, chain_merge2);
        Alignment alignment = stitcher.stitch(anchors,
                                              subproblem1.graph, subproblem2.graph,
                                              subproblem1.tableau, subproblem2.tableau,
                                              chain_merge1, chain_merge2);
        
        // record the results of this subproblem
        auto& next_problem = subproblems[node_id];
        if (preserve_subproblems) {
            next_problem.graph = subproblem1.graph;
        }
        else {
            next_problem.graph = move(subproblem1.graph);
        }
        next_problem.tableau = subproblem1.tableau;
        
        fuse(next_problem.graph, subproblem2.graph,
             next_problem.tableau, subproblem2.tableau,
             alignment);
        
        // determinize the graph
        BaseGraph determinized = determinize(next_problem.graph);
        next_problem.tableau = translate_tableau(determinized, next_problem.tableau);
        rewalk_paths(determinized, next_problem.tableau, next_problem.graph);
        next_problem.graph = move(determinized);
    }
    
}

}
