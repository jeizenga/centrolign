
#include "centrolign/core.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <cassert>
#include <cmath>

#include "centrolign/chain_merge.hpp"
#include "centrolign/path_merge.hpp"
#include "centrolign/utility.hpp"
#include "centrolign/determinize.hpp"
#include "centrolign/logging.hpp"
#include "centrolign/gfa.hpp"
#include "centrolign/fuse.hpp"

namespace centrolign {

using namespace std;


Core::Core(const std::string& fasta_file, const std::string& tree_file) :
    anchorer(score_function), partitioner(score_function)
{
    
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
           Tree&& tree) :
    anchorer(score_function), partitioner(score_function)
{
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
        
        // get rid of samples we don't have the sequence for
        tree.prune(sequence_leaf_ids);
    }
    
    // get rid of non-branching paths
    tree.compact();
    
    if (logging::level >= logging::Debug) {
        Tree polytomized = tree;
        polytomized.polytomize();
        logging::log(logging::Debug, "Simplified and subsetted tree:\n" + polytomized.to_newick());
    }
    
    // convert into a binary tree
    tree.binarize();
    
    logging::log(logging::Debug, "Fully processed tree:\n" + tree.to_newick());
    
    logging::log(logging::Basic, "Initializing leaf subproblems.");
    
    subproblems.resize(tree.node_size());
    for (uint64_t node_id = 0; node_id < tree.node_size(); ++node_id) {
        if (tree.is_leaf(node_id)) {
            const auto& name = tree.label(node_id);
            const auto& sequence = sequences[name_to_idx[name]].second;
            
            auto& subproblem = subproblems[node_id];
            
            subproblem.graph = make_base_graph(name, sequence);
            subproblem.tableau = add_sentinels(subproblem.graph, 5, 6);
            subproblem.name = name;
            subproblem.complete = true;
        }
    }
}

std::vector<match_set_t> Core::get_matches(Subproblem& subproblem1, Subproblem& subproblem2,
                                           bool suppress_verbose_logging) const {
    
    // give them unique sentinels for anchoring
    reassign_sentinels(subproblem1.graph, subproblem1.tableau, 5, 6);
    reassign_sentinels(subproblem2.graph, subproblem2.tableau, 7, 8);
    
    ExpandedGraph expanded1;
    ExpandedGraph expanded2;
    if (match_finder.path_matches) {
        // no need to simplify if only querying paths, just borrow the graphs
        expanded1.graph = move(subproblem1.graph);
        expanded1.tableau = subproblem1.tableau;
        expanded2.graph = move(subproblem2.graph);
        expanded2.tableau = subproblem2.tableau;
    }
    else {
        logging::log(suppress_verbose_logging ? logging::Debug : logging::Verbose, "Simplifying complex regions.");
        
        // simplify complex regions of the graph
        expanded1 = move(simplifier.simplify(subproblem1.graph, subproblem1.tableau));
        expanded2 = move(simplifier.simplify(subproblem2.graph, subproblem2.tableau));
    }
    
    logging::log(suppress_verbose_logging ? logging::Debug : logging::Verbose , "Finding matches.");
    
    // find minimal rare matches
    auto matches = query_matches(expanded1, expanded2);
    
    if (match_finder.path_matches) {
        // return the graphs that we borrowed
        subproblem1.graph = move(expanded1.graph);
        subproblem2.graph = move(expanded2.graph);
    }
    
    return matches;
}

void Core::execute() {
    
    if (!skip_calibration) {
        calibrate_anchor_scores();
    }
    
    logging::log(logging::Minimal, "Beginning MSA.");
    
    for (auto node_id : tree.postorder()) {
        
        if (tree.is_leaf(node_id)) {
            continue;
        }
        
        logging::log(logging::Basic, "Executing next subproblem.");
        if (logging::level >= logging::Verbose) {
            
            stringstream strm;
            strm << "Contains sequences:\n";
            for (auto leaf_name : leaf_descendents(node_id)) {
                strm << '\t' << leaf_name << '\n';
            }
            logging::log(logging::Verbose, strm.str());
        }
        
        auto& next_problem = subproblems[node_id];
        
        if (next_problem.complete) {
            logging::log(logging::Verbose, "Problem already finished from restarted run.");
            continue;
        }
        
        const auto& children = tree.get_children(node_id);
        assert(children.size() == 2);
        
        auto& subproblem1 = subproblems[children.front()];
        auto& subproblem2 = subproblems[children.back()];
        
        auto matches = get_matches(subproblem1, subproblem2, false);
        
        {
            logging::log(logging::Verbose, "Computing reachability.");
            
            if (anchorer.chaining_algorithm == Anchorer::SparseAffine) {
                // use all paths for reachability to get better distance estimates
                
                PathMerge path_merge1(subproblem1.graph, subproblem1.tableau);
                PathMerge path_merge2(subproblem2.graph, subproblem2.tableau);
                
                next_problem.alignment = std::move(align(matches, subproblem1, subproblem2,
                                                         path_merge1, path_merge2));
                
            }
            else {
                // use non-overlapping chains for reachability for more efficiency
                ChainMerge chain_merge1(subproblem1.graph, subproblem1.tableau);
                ChainMerge chain_merge2(subproblem2.graph, subproblem2.tableau);
                
                next_problem.alignment = std::move(align(matches, subproblem1, subproblem2,
                                                         chain_merge1, chain_merge2));
            }
        }
        
        logging::log(logging::Verbose, "Fusing MSAs along the alignment.");
        
        // fuse either in place or in a copy
        BaseGraph fused_graph;
        if (preserve_subproblems) {
            fused_graph = subproblem1.graph;
        }
        else {
            fused_graph = move(subproblem1.graph);
        }
                
        fuse(fused_graph, subproblem2.graph,
             subproblem1.tableau, subproblem2.tableau,
             next_problem.alignment);
        
        if (!preserve_subproblems) {
            // we no longer need this, clobber it to save memory
            subproblem2.graph = BaseGraph();
        }
        
        if (match_finder.path_matches) {
            // just save the results
            next_problem.graph = move(fused_graph);
            next_problem.tableau = subproblem1.tableau;
        }
        else {
            logging::log(logging::Verbose, "Determinizing the fused MSA.");
            
            // determinize the graph and save the results
            next_problem.graph = determinize(fused_graph);
            next_problem.tableau = translate_tableau(next_problem.graph, subproblem1.tableau);
            rewalk_paths(next_problem.graph, next_problem.tableau, fused_graph);
            purge_uncovered_nodes(next_problem.graph, next_problem.tableau);
        }
        
        next_problem.complete = true;
        
        if (!subproblems_prefix.empty()) {
            emit_subproblem(node_id);
        }
    }
}

std::vector<std::string> Core::leaf_descendents(uint64_t tree_id) const {
    std::vector<std::string> descendents;
    std::vector<uint64_t> stack = tree.get_children(tree_id);
    while (!stack.empty()) {
        auto here = stack.back();
        stack.pop_back();
        if (tree.is_leaf(here)) {
            descendents.push_back(tree.label(here));
        }
        else {
            for (auto next : tree.get_children(here)) {
                stack.push_back(next);
            }
        }
    }
    return descendents;
}

void Core::calibrate_anchor_scores() {
    
    logging::log(logging::Basic, "Calibrating scale of anchoring parameters.");
    
    std::vector<double> intrinsic_scales;
    
    size_t num_leaves = (tree.node_size() + 1) / 2;
    size_t leaf_num = 1;
    for (uint64_t tree_id = 0; tree_id < tree.node_size(); ++tree_id) {
        if (!tree.is_leaf(tree_id)) {
            continue;
        }
        
        logging::log(logging::Verbose, "Estimating scale for sequence " + to_string(leaf_num++) + " of " + to_string(num_leaves) + ".");
        
        auto& subproblem = subproblems[tree_id];
        
        std::vector<match_set_t> matches;
        {
            // TODO: have to copy it for the irritating necessity of having distinct sentinels
            // maybe it would be okay to have the sentinels slip into the anchors for this case?
            auto subproblem_copy = subproblem;
            
            auto full_matches = get_matches(subproblem, subproblem_copy, true);
            
            // subset down to only matches on the main diagonal (retaining count for scoring)
            matches.reserve(full_matches.size());
            for (auto& match_set : full_matches) {
                for (auto& walk : match_set.walks1) {
                    matches.emplace_back();
                    auto& match = matches.back();
                    match.walks1.emplace_back(walk);
                    match.walks2.emplace_back(std::move(walk));
                    match.count1 = match_set.count1;
                    match.count2 = match_set.count2;
                }
            }
        }
        
        
        ChainMerge chain_merge(subproblem.graph, subproblem.tableau);
        
        double scale = anchorer.estimate_score_scale(matches, subproblem.graph, subproblem.graph,
                                                     subproblem.tableau, subproblem.tableau,
                                                     chain_merge, chain_merge);
        
        logging::log(logging::Debug, "Compute intrinsic scale of " + std::to_string(scale) + " for sequence " + tree.label(tree_id));
        intrinsic_scales.push_back(scale);
    }
    
    double mean = 0.0;
    for (auto scale : intrinsic_scales) {
        mean += scale;
    }
    mean /= intrinsic_scales.size();
    double var = 0.0;
    for (auto scale : intrinsic_scales) {
        double dev = scale - mean;
        var += dev * dev;
    }
    var /= intrinsic_scales.size() - 1;
    double std_dev = sqrt(var);
    
    logging::log(logging::Verbose, "Intrinsic sequence scales are centered at " + std::to_string(mean) + " +/- " + std::to_string(std_dev) + ".");
    
    partitioner.score_scale = mean;
}

std::string Core::subproblem_info_file_name() const {
    return subproblems_prefix + "_info.txt";
}

std::string Core::subproblem_file_name(uint64_t tree_id) const {
    auto seq_names = leaf_descendents(tree_id);
    sort(seq_names.begin(), seq_names.end());
    
    // create a hash digest of the sample names
    size_t hsh = 660422875706093811ull;
    for (auto& seq_name : seq_names) {
        hash_combine(hsh, size_t(2110260111091729000ull)); // spacer
        for (char c : seq_name) {
            hash_combine(hsh, c);
        }
    }
    
    string file_name = subproblems_prefix + "_" + to_hex(hsh) + ".gfa";
    return file_name;
}

void Core::emit_subproblem(uint64_t tree_id) const {
    
    auto gfa_file_name = subproblem_file_name(tree_id);
    auto info_file_name = subproblem_info_file_name();
    
    // check if the file already exists
    bool write_header = !(ifstream(info_file_name));
    
    ofstream info_out(info_file_name, ios_base::app);
    if (!info_out) {
        throw std::runtime_error("Failed to write to subproblem info file " + info_file_name);
    }
    ofstream gfa_out(gfa_file_name);
    if (!gfa_out) {
        throw std::runtime_error("Failed to write to subproblem file " + gfa_file_name);
    }
    if (!info_out) {
        throw std::runtime_error("Failed to write to subproblem info file " + info_file_name);
    }
    
    if (write_header) {
        info_out << "filename\tsequences\n";
    }
    auto sequences = leaf_descendents(tree_id);
    sort(sequences.begin(), sequences.end());
    info_out << gfa_file_name << '\t' << join(sequences, ",") << '\n';
    
    write_gfa(subproblems[tree_id].graph, subproblems[tree_id].tableau, gfa_out);
}

std::vector<match_set_t> Core::query_matches(ExpandedGraph& expanded1,
                                             ExpandedGraph& expanded2) const {
    
    std::vector<match_set_t> matches;
    try {
        if (expanded1.back_translation.empty() && expanded2.back_translation.empty()) {
            matches = move(match_finder.find_matches(expanded1.graph, expanded2.graph,
                                                     expanded1.tableau, expanded2.tableau));
        }
        else {
            matches = move(match_finder.find_matches(expanded1.graph, expanded2.graph,
                                                     expanded1.tableau, expanded2.tableau,
                                                     expanded1.back_translation,
                                                     expanded2.back_translation));
        }
    } catch (GESASizeException& ex) {

        logging::log(logging::Verbose, "Graph not simple enough to index, resimplifying.");

        auto targets = simplifier.identify_target_nodes(ex.from_counts());

        size_t simplify_dist = (1 << ex.doubling_step());

        size_t pre_simplify_size1 = expanded1.graph.node_size();
        size_t pre_simplify_size2 = expanded2.graph.node_size();

        auto expanded_more1 = simplifier.targeted_simplify(expanded1.graph, expanded1.tableau,
                                                           targets[0], simplify_dist);
        auto expanded_more2 = simplifier.targeted_simplify(expanded2.graph, expanded2.tableau,
                                                           targets[1], simplify_dist);
        
        for (auto& tr : expanded_more1.back_translation) {
            tr = expanded1.back_translation[tr];
        }
        for (auto& tr : expanded_more2.back_translation) {
            tr = expanded2.back_translation[tr];
        }
        
        expanded1 = std::move(expanded_more1);
        expanded2 = std::move(expanded_more2);

        if (pre_simplify_size1 == expanded1.graph.node_size() &&
            pre_simplify_size2 == expanded2.graph.node_size()) {

            throw std::runtime_error("Simplification algorithm failed to simplify graph");
        }

        // recursively try again with a more simplified graph
        matches = move(query_matches(expanded1, expanded2));
    }
    
    return matches;
}

void Core::restart() {
    
    int num_restarted = 0;
    for (auto node_id : tree.preorder()) {
        
        if (subproblems[node_id].complete) {
            continue;
        }
        
        auto file_name = subproblem_file_name(node_id);
        
        ifstream gfa_in(file_name);
        
        if (gfa_in) {
            // we have the results of this alignment problem saved
            
            ++num_restarted;
            logging::log(logging::Debug, "Loading previously completed subproblem " + file_name);
            
            auto& subproblem = subproblems[node_id];
            
            subproblem.graph = read_gfa(gfa_in);
            subproblem.tableau = add_sentinels(subproblem.graph, 5, 6);
            subproblem.complete = true;
            // FIXME: we dont' save the subproblems alignments, but for now that's not a problem
            //subproblem.alignment = ?
            
            // mark descendents as complete
            auto stack = tree.get_children(node_id);
            while (!stack.empty()) {
                auto top = stack.back();
                stack.pop_back();
                subproblems[top].complete = true;
                if (!preserve_subproblems) {
                    // clear out descendents
                    subproblems[top].graph = BaseGraph();
                }
                for (auto child_id : tree.get_children(top)) {
                    stack.push_back(child_id);
                }
            }
            
        }
    }
    
    logging::log(logging::Basic, "Loaded results for " + std::to_string(num_restarted) + " subproblem(s) from previously completed run.");
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
