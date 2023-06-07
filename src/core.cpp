#include "centrolign/core.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <cassert>

#include "centrolign/utility.hpp"
#include "centrolign/determinize.hpp"
#include "centrolign/logging.hpp"

namespace centrolign {

using namespace std;


Core::Core(const std::string& fasta_file, const std::string& tree_file) {
    
    ifstream fasta_fstream, tree_fstream;
    
    istream* fasta_stream = get_input(fasta_file, fasta_fstream);
    istream* tree_stream = get_input(tree_file, tree_fstream);
    
    auto sequences = parse_fasta(*fasta_stream);
    
    // read the newick stream into a tree
    stringstream sstrm;
    sstrm << tree_stream->rdbuf();
    Tree parsed_tree(sstrm.str());
    
    init(move(sequences), move(parsed_tree));
}

Core::Core(std::vector<std::pair<std::string, std::string>>&& names_and_sequences,
           Tree&& tree) {
    
    init(move(names_and_sequences), move(tree));
}

void Core::init(std::vector<std::pair<std::string, std::string>>&& names_and_sequences,
                Tree&& tree_in) {
    
    auto sequences = move(names_and_sequences);
    tree = move(tree_in);
    
    unordered_map<string, size_t> name_to_idx;
    for (size_t i = 0; i < sequences.size(); ++i) {
        auto& name = sequences[i].first;
        if (name_to_idx.count(name)) {
            throw runtime_error("FASTA contains duplicate name " + name);
        }
        name_to_idx[name] = i;
    }
    
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
    
    logging::log(logging::Basic, "Initializing leaf subproblems");
    
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
    
    logging::log(logging::Minimal, "Beginning MSA");
    
    for (auto node_id : tree.postorder()) {
        if (tree.is_leaf(node_id)) {
            continue;
        }
        
        const auto& children = tree.get_children(node_id);
        assert(children.size() == 2);
        
        auto& subproblem1 = subproblems[children.front()];
        auto& subproblem2 = subproblems[children.back()];
        
        logging::log(logging::Basic, "Executing next subproblem");
        if (logging::level >= logging::Verbose) {
            std::vector<uint64_t> leaf_descendents;
            std::vector<uint64_t> stack = tree.get_children(node_id);
            while (!stack.empty()) {
                auto here = stack.back();
                stack.pop_back();
                if (tree.is_leaf(here)) {
                    leaf_descendents.push_back(here);
                }
                else {
                    for (auto next : tree.get_children(here)) {
                        stack.push_back(next);
                    }
                }
            }
            stringstream strm;
            strm << "Contains sequences:\n";
            for (auto leaf_id : leaf_descendents) {
                strm << '\t' << tree.label(leaf_id) << '\n';
            }
            logging::log(logging::Verbose, strm.str());
        }
        
        // give them unique sentinels for anchoring
        reassign_sentinels(subproblem1.graph, subproblem1.tableau, 5, 6);
        reassign_sentinels(subproblem2.graph, subproblem2.tableau, 7, 8);
        
        logging::log(logging::Verbose, "Simplifying complex regions");
        
        // simplify complex regions of the graph
        auto expanded1 = simplifier.simplify(subproblem1.graph, subproblem1.tableau);
        auto expanded2 = simplifier.simplify(subproblem2.graph, subproblem2.tableau);
        
        logging::log(logging::Verbose, "Finding matches");
        
        // find minimal rare matches
        auto matches = std::move(find_matches(expanded1, expanded2));
        
        logging::log(logging::Verbose, "Computing reachability");
        
        // compute chain merge data structures for later steps
        ChainMerge chain_merge1(subproblem1.graph, subproblem1.tableau);
        ChainMerge chain_merge2(subproblem2.graph, subproblem2.tableau);
        
        logging::log(logging::Verbose, "Anchoring");
        
        // anchor the alignment
        auto anchors = anchorer.anchor_chain(matches,
                                             subproblem1.graph, subproblem2.graph,
                                             chain_merge1, chain_merge2);
        
        logging::log(logging::Verbose, "Stitching anchors into alignment");
        
        // form a base-level alignment
        auto& next_problem = subproblems[node_id];
        next_problem.alignment = stitcher.stitch(anchors,
                                                 subproblem1.graph, subproblem2.graph,
                                                 subproblem1.tableau, subproblem2.tableau,
                                                 chain_merge1, chain_merge2);
        
        logging::log(logging::Verbose, "Fusing MSAs along the alignment");
        // record the results of this subproblem
        if (preserve_subproblems) {
            next_problem.graph = subproblem1.graph;
        }
        else {
            next_problem.graph = move(subproblem1.graph);
        }
        next_problem.tableau = subproblem1.tableau;
        
        fuse(next_problem.graph, subproblem2.graph,
             next_problem.tableau, subproblem2.tableau,
             next_problem.alignment);
        
        logging::log(logging::Verbose, "Determinizing the fused MSA");
        // determinize the graph
        BaseGraph determinized = determinize(next_problem.graph);
        next_problem.tableau = translate_tableau(determinized, next_problem.tableau);
        rewalk_paths(determinized, next_problem.tableau, next_problem.graph);
        next_problem.graph = move(determinized);
    }
}


std::vector<match_set_t>&& Core::find_matches(ExpandedGraph& expanded1,
                                              ExpandedGraph& expanded2) const {
    try {
        return std::move(match_finder.find_matches(expanded1.graph, expanded2.graph,
                                                   expanded1.back_translation,
                                                   expanded2.back_translation));
    } catch (GESASizeException& ex) {
        
        auto targets = simplifier.identify_target_nodes(ex.from_counts());
        
        size_t simplify_dist = (1 << ex.doubling_step());
        
        size_t pre_simplify_size1 = expanded1.graph.node_size();
        size_t pre_simplify_size2 = expanded2.graph.node_size();
        
        expanded1 = simplifier.targeted_simplify(expanded1.graph, expanded1.tableau,
                                                 targets[0], simplify_dist);
        expanded2 = simplifier.targeted_simplify(expanded2.graph, expanded2.tableau,
                                                 targets[1], simplify_dist);
        
        if (pre_simplify_size1 == expanded1.graph.node_size() &&
            pre_simplify_size2 == expanded2.graph.node_size()) {
            
            throw std::runtime_error("Simplification algorithm failed to simplify graph");
        }
        
        // recursively try again with a more simplified graph
        return std::move(find_matches(expanded1, expanded2));
    }
}

const Core::Subproblem& Core::root_subproblem() const {
    return subproblems[tree.get_root()];
}

const Core::Subproblem& Core::leaf_subproblem(const std::string& name) const {
    return subproblems[tree.get_id(name)];
}

const Core::Subproblem& Core::subproblem_covering(const std::vector<string>& names) const {
    // TODO: i'm sure there's an O(n) algorithm for this, but hopefully this O(n^2) is
    // good enough
    
    if (names.empty()){
        throw runtime_error("Cannor query an MSA subproblem with an empty set of sequences");
    }
    
    // walk up to the root on an arbitrary node
    vector<uint64_t> path_to_root;
    {
        uint64_t node_id = tree.get_id(names.front());
        path_to_root.push_back(node_id);
        while (node_id != tree.get_root()) {
            node_id = tree.get_parent(node_id);
            path_to_root.push_back(node_id);
        }
    }
    // we'll be removing from the bottom of this path
    reverse(path_to_root.begin(), path_to_root.end());
    
    // walk up to the root on each other node
    for (size_t i = 1; i < names.size(); ++i) {
        // TODO: repetitive
        
        unordered_set<uint64_t> next_path_to_root;
        
        uint64_t node_id = tree.get_id(names[i]);
        next_path_to_root.insert(node_id);
        while (node_id != tree.get_root()) {
            node_id = tree.get_parent(node_id);
            next_path_to_root.insert(node_id);
        }
        
        // remove from the first path until we find a shared node
        while (!next_path_to_root.count(path_to_root.back())) {
            path_to_root.pop_back();
        }
    }
    
    return subproblems[path_to_root.back()];
}


}
