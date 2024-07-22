
#include "centrolign/core.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <cassert>
#include <cmath>
#include <limits>

#include "centrolign/chain_merge.hpp"
#include "centrolign/path_merge.hpp"
#include "centrolign/utility.hpp"
#include "centrolign/determinize.hpp"
#include "centrolign/logging.hpp"
#include "centrolign/gfa.hpp"
#include "centrolign/fuse.hpp"
#include "centrolign/step_index.hpp"

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
    match_finder(score_function), anchorer(score_function), partitioner(score_function)
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
    
    log_memory_usage(logging::Debug);
    
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
    
    log_memory_usage(logging::Debug);
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
    
    std::vector<Alignment> bond_alignments;
    if (!skip_calibration || cyclize_tandem_duplications) {
        bond_alignments = std::move(calibrate_anchor_scores_and_identify_bonds());
    }
    
    logging::log(logging::Minimal, "Beginning MSA.");
    log_memory_usage(logging::Debug);
    
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
        
        if (logging::level >= logging::Debug) {
            size_t graph_mem_size = 0;
            size_t aln_mem_size = 0;
            for (const auto& subproblem : subproblems) {
                graph_mem_size += subproblem.graph.memory_size();
                aln_mem_size += subproblem.alignment.capacity() * sizeof(decltype(subproblem.alignment)::value_type);
            }
            logging::log(logging::Debug, "In-memory graphs and alignments are occupying " + format_memory_usage(graph_mem_size) + " and " + format_memory_usage(aln_mem_size) + ", respectively.");
            logging::log(logging::Debug, "Current memory use is " + format_memory_usage(current_memory_usage()));
        }
        
        const auto& children = tree.get_children(node_id);
        assert(children.size() == 2);
        
        auto& subproblem1 = subproblems[children.front()];
        auto& subproblem2 = subproblems[children.back()];
        
        auto matches = get_matches(subproblem1, subproblem2, false);

        log_memory_usage(logging::Debug);

        {
            logging::log(logging::Verbose, "Computing reachability.");
            
            if (anchorer.chaining_algorithm == Anchorer::SparseAffine) {
                // use all paths for reachability to get better distance estimates
                
                #define _gen_path_merge(UIntSize, UIntChain) \
                    PathMerge<UIntSize, UIntChain> path_merge1(subproblem1.graph, subproblem1.tableau); \
                    PathMerge<UIntSize, UIntChain> path_merge2(subproblem2.graph, subproblem2.tableau); \
                    next_problem.alignment = std::move(align(matches, subproblem1, subproblem2, \
                                                             path_merge1, path_merge2))
                
                size_t max_nodes = std::max(subproblem1.graph.node_size(), subproblem2.graph.node_size());
                size_t max_paths = std::max(subproblem1.graph.path_size(), subproblem2.graph.path_size());
                if (max_nodes < std::numeric_limits<uint32_t>::max() && max_paths < std::numeric_limits<uint8_t>::max()) {
                    _gen_path_merge(uint32_t, uint8_t);
                }
                else if (max_nodes < std::numeric_limits<uint32_t>::max()) {
                    _gen_path_merge(uint32_t, uint16_t);
                }
                else {
                    _gen_path_merge(uint64_t, uint16_t);
                }
                
                #undef _gen_path_merge
                
            }
            else {
                // use non-overlapping chains for reachability for more efficiency
                ChainMerge chain_merge1(subproblem1.graph, subproblem1.tableau);
                ChainMerge chain_merge2(subproblem2.graph, subproblem2.tableau);
                                
                next_problem.alignment = std::move(align(matches, subproblem1, subproblem2,
                                                         chain_merge1, chain_merge2));
            }
        }
        
        log_memory_usage(logging::Debug);
        
        // we do this now in case we're not preserving the graphs in the subproblems
        if (!subalignments_filepath.empty()) {
            emit_subalignment(node_id);
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
            // we no longer need these, clobber them to save memory
            BaseGraph dummy_graph = std::move(subproblem2.graph);
            Alignment dummy_aln1 = std::move(subproblem1.alignment);
            Alignment dummy_aln2 = std::move(subproblem2.alignment);
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
            
            logging::log(logging::Verbose, "Fusing MSAs along the alignment.");
        }
        
        next_problem.complete = true;
        
        if (!subproblems_prefix.empty()) {
            emit_subproblem(node_id);
        }
        
        log_memory_usage(logging::Verbose);
    }
    
    if (cyclize_tandem_duplications) {
        apply_bonds(bond_alignments);
    }
}

std::vector<std::string> Core::leaf_descendents(uint64_t tree_id) const {
    std::vector<std::string> descendents;
    if (tree.is_leaf(tree_id)) {
        descendents.push_back(tree.label(tree_id));
    }
    else {
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
    }
    return descendents;
}

std::vector<Alignment> Core::calibrate_anchor_scores_and_identify_bonds() {
    
    std::string msg;
    if (cyclize_tandem_duplications && !skip_calibration) {
        msg = "Calibrating scale of anchoring parameters and identifying tandem duplications.";
    }
    else if (cyclize_tandem_duplications) {
        msg = "Identifying tandem duplications.";
    }
    else {
        msg = "Calibrating scale of anchoring parameters.";
    }
    logging::log(logging::Basic, msg);
    log_memory_usage(logging::Debug);
    
    std::vector<double> intrinsic_scales;
    
    std::vector<Alignment> bond_alns;
    
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
        
        PathMerge<> path_merge(subproblem.graph, subproblem.tableau);
        
        std::vector<anchor_t> chain;
        double scale = anchorer.estimate_score_scale(matches, subproblem.graph, subproblem.graph,
                                                     subproblem.tableau, subproblem.tableau,
                                                     path_merge, path_merge, &chain);
        
        logging::log(logging::Debug, "Compute intrinsic scale of " + std::to_string(scale) + " for sequence " + tree.label(tree_id));
        log_memory_usage(logging::Debug);
        
        if (!skip_calibration) {
            intrinsic_scales.push_back(scale);
        }
        
        if (cyclize_tandem_duplications) {
            
            // we need to reanchor and
            
            auto mask = generate_diagonal_mask(matches);
            
            for (size_t iter = 0; iter < max_tandem_duplication_search_rounds; ++iter) {
                
                logging::log(logging::Verbose, "Beginning round " + std::to_string(iter) + " of tandem duplication detection for sequence " + tree.label(tree_id) + ".");
                
                // get the next-best unmasked chain
                auto secondary_chain = anchorer.anchor_chain(matches, subproblem.graph, subproblem.graph,
                                                             subproblem.tableau, subproblem.tableau,
                                                             path_merge, path_merge, &mask);
                
                // identify high-enough scoring segments
                auto bonds = bonder.identify_bonds(subproblem.graph, subproblem.graph,
                                                   subproblem.tableau, subproblem.tableau,
                                                   path_merge, path_merge,
                                                   chain, secondary_chain);
                
                logging::log(logging::Verbose, "Found " + std::to_string(bonds.size()) + " tandem duplications.");
                
                if (bonds.empty()) {
                    // if we didn't find any bonds this round, we're unlikely to in the future
                    break;
                }
                
                // stitch the bonds into alignments
                for (auto& bond : bonds) {
                    
                    auto bond_chain = bonds_to_chain(subproblem.graph, bond);
                    
                    bond_alns.emplace_back(stitcher.internal_stitch(bond_chain, subproblem.graph, path_merge));
                }
                
                if (iter != max_tandem_duplication_search_rounds) {
                    // mask out anchors that overlap this chain's matches
                    update_mask(matches, secondary_chain, mask, true);
                }
            }
        }
    }
    
    if (!skip_calibration) {
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
        
        score_function.score_scale = mean;
    }
    
    return bond_alns;
}

std::unordered_set<std::tuple<size_t, size_t, size_t>> Core::generate_diagonal_mask(const std::vector<match_set_t>& matches) const {
    
    std::unordered_set<std::tuple<size_t, size_t, size_t>> mask;
    for (size_t i = 0; i < matches.size(); ++i) {
        
        const auto& match_set = matches[i];
        std::unordered_map<uint64_t, size_t> start_to_idx;
        for (size_t j = 0; j < match_set.walks1.size(); ++j) {
            start_to_idx[match_set.walks1[j].front()] = j;
        }
        for (size_t k = 0; k < match_set.walks2.size(); ++k) {
            auto it = start_to_idx.find(match_set.walks2[k].front());
            if (it != start_to_idx.end()) {
                mask.emplace(i, it->second, k);
            }
        }
    }
    
    return mask;
}

void Core::update_mask(const std::vector<match_set_t>& matches, const std::vector<anchor_t>& chain,
                       std::unordered_set<std::tuple<size_t, size_t, size_t>>& masked_matches, bool mask_reciprocal) const {
    
    logging::log(logging::Debug, "Updating mask");
    
    // record the node pairings that we're going to mask
    std::unordered_map<uint64_t, uint64_t> paired_node_ids;
    for (const auto& anchor : chain) {
        for (size_t i = 0; i < anchor.walk1.size(); ++i) {
            paired_node_ids[anchor.walk1[i]] = anchor.walk2[i];
            if (mask_reciprocal) {
                paired_node_ids[anchor.walk2[i]] = anchor.walk1[i];
            }
        }
    }
    
    for (size_t i = 0; i < matches.size(); ++i) {
        
        const auto& match_set = matches[i];
        
        // for each position in walk, for each node ID, the walk indexes that it matches
        std::vector<std::unordered_map<uint64_t, std::vector<size_t>>> walk2_node(match_set.walks1.front().size());
        for (size_t k = 0; k < match_set.walks2.size(); ++k) {
            const auto& walk2 = match_set.walks2[k];
            for (size_t l = 0; l < walk2.size(); ++l) {
                walk2_node[l][walk2[l]].push_back(k);
            }
        }
        
        // check for walk2s that have a paired node ID in the matching position
        for (size_t j = 0; j < match_set.walks1.size(); ++j) {
            const auto& walk1 = match_set.walks1[j];
            for (size_t l = 0; l < walk1.size(); ++l) {
                auto it = paired_node_ids.find(walk1[l]);
                if (it != paired_node_ids.end()) {
                    auto it2 = walk2_node[l].find(it->second);
                    if (it2 != walk2_node[l].end()) {
                        for (auto k : it2->second) {
                            masked_matches.emplace(i, j, k);
                        }
                    }
                }
            }
        }
    }
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

void Core::emit_subalignment(uint64_t tree_id) const {
    
    // TODO: this is fragile in that we need to separately use the same order for
    // the children here and in the alignment routine
    const auto& subproblem = subproblems[tree_id];
    const auto& children = tree.get_children(tree_id);
    const auto& graph1 = subproblems[children.front()].graph;
    const auto& graph2 = subproblems[children.back()].graph;
    
    ofstream out(subalignments_filepath, ios_base::app);
    if (!out) {
        throw std::runtime_error("Failed to write to subalignment file " + subalignments_filepath);
    }
    
    out << "# sequence set 1\n";
    for (const auto& seq_name : leaf_descendents(children.front())) {
        out << seq_name << '\n';
    }
    out << "# sequence set 2\n";
    for (const auto& seq_name : leaf_descendents(children.back())) {
        out << seq_name << '\n';
    }
    
    StepIndex step_index1(graph1);
    StepIndex step_index2(graph2);
    out << "# alignment\n";
    for (const auto& aln_pair : subproblem.alignment) {
        if (aln_pair.node_id1 == AlignedPair::gap) {
            out << "-\t-\t-";
        }
        else {
            uint64_t path_id;
            size_t step;
            tie(path_id, step) = step_index1.path_steps(aln_pair.node_id1).front();
            out << graph1.path_name(path_id) << '\t' << step << '\t' << decode_base(graph1.label(graph1.path(path_id)[step]));
        }
        out << '\t';
        if (aln_pair.node_id2 == AlignedPair::gap) {
            out << "-\t-\t-";
        }
        else {
            uint64_t path_id;
            size_t step;
            tie(path_id, step) = step_index2.path_steps(aln_pair.node_id2).front();
            out << graph2.path_name(path_id) << '\t' << step << '\t' << decode_base(graph2.label(graph2.path(path_id)[step]));
        }
        out << '\n';
        
    }
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
    }
    catch (GESASizeException& ex) {

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

void Core::apply_bonds(const std::vector<Alignment>& bond_alignments) {
    
    auto& root_subproblem = subproblems[tree.get_root()];
    
    SentinelTableau cyclized_tableau;
    BaseGraph cyclized = internal_fuse(root_subproblem.graph, bond_alignments, &root_subproblem.tableau, &cyclized_tableau);
    
    // TODO: clean up the alignment
    
    root_subproblem.graph = std::move(cyclized);
    root_subproblem.tableau = cyclized_tableau;
}

void Core::log_memory_usage(logging::LoggingLevel level) const {
    if (logging::level >= level) {
        int64_t peak_mem = max_memory_usage();
        if (peak_mem < 0) {
            logging::log(level, "Failed to measure peak memory usage.");
        }
        else {
            logging::log(level, "Peak memory usage so far: " + format_memory_usage(peak_mem) + ".");
        }
        if (level == logging::Debug) {
            int64_t curr_mem = current_memory_usage();
            if (curr_mem < 0) {
                logging::log(level, "Failed to measure current memory usage.");
            }
            else {
                logging::log(level, "Current memory usage: " + format_memory_usage(curr_mem) + ".");
            }
        }
    }
}

void Core::restart() {
    
    int num_restarted = 0;
    int num_pruned = 0;
    for (auto node_id : tree.preorder()) {
        
        if (subproblems[node_id].complete) {
            if (!tree.is_leaf(node_id)) {
                ++num_pruned;
            }
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
            log_memory_usage(logging::Debug);
            
            // mark descendents as complete
            auto stack = tree.get_children(node_id);
            while (!stack.empty()) {
                auto top = stack.back();
                stack.pop_back();
                subproblems[top].complete = true;
                if (!tree.is_leaf(top)) {
                    ++num_pruned;
                }
                if (!preserve_subproblems) {
                    // clear out descendents
                    BaseGraph dummy = std::move(subproblems[top].graph);
                }
                for (auto child_id : tree.get_children(top)) {
                    stack.push_back(child_id);
                }
            }
            
        }
    }
    
    logging::log(logging::Basic, "Loaded results for " + std::to_string(num_restarted) + " subproblem(s) from previously completed run and pruned " + std::to_string(num_pruned) + " of their children as unnecessary.");
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
