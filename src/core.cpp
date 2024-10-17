
#include "centrolign/core.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <cassert>
#include <cmath>
#include <limits>
#include <functional>

#include "centrolign/chain_merge.hpp"
#include "centrolign/path_merge.hpp"
#include "centrolign/utility.hpp"
#include "centrolign/determinize.hpp"
#include "centrolign/logging.hpp"
#include "centrolign/gfa.hpp"
#include "centrolign/fuse.hpp"
#include "centrolign/step_index.hpp"
#include "centrolign/induced_match_finder.hpp"

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
    path_match_finder(score_function), anchorer(score_function), partitioner(score_function)
{
    init(move(names_and_sequences), move(tree));
}

void Core::init(std::vector<std::pair<std::string, std::string>>&& names_and_sequences,
                Tree&& tree_in) {
    
    main_execution = std::move(Execution(std::move(names_and_sequences), std::move(tree_in)));
}


void Core::execute() {
    
    std::vector<std::pair<std::string, Alignment>> bond_alignments;
    if (!skip_calibration || (cyclize_tandem_duplications && !restarted_bond_alignments.get())) {
        bond_alignments = std::move(calibrate_anchor_scores_and_identify_bonds());
    }
    
    if (restarted_bond_alignments.get()) {
        bond_alignments = std::move(*restarted_bond_alignments.get());
    }
    else if (cyclize_tandem_duplications && !subproblems_prefix.empty()) {
        emit_restart_bonds(bond_alignments);
    }
    
    logging::log(logging::Minimal, "Beginning MSA.");
    log_memory_usage(logging::Debug);
    
    do_execution(main_execution, this->path_match_finder, true);
    
    if (!induced_pairwise_prefix.empty()) {
        output_pairwise_alignments(false);
    }
    
    if (cyclize_tandem_duplications) {
        apply_bonds(bond_alignments);
        
        if (!induced_pairwise_prefix.empty()) {
            output_pairwise_alignments(true);
        }
    }
}

std::vector<std::pair<std::string, Alignment>> Core::calibrate_anchor_scores_and_identify_bonds() {
    
    std::string msg;
    if ((cyclize_tandem_duplications && !restarted_bond_alignments.get()) && !skip_calibration) {
        msg = "Calibrating scale of anchoring parameters and identifying tandem duplications.";
    }
    else if (cyclize_tandem_duplications && !restarted_bond_alignments.get()) {
        msg = "Identifying tandem duplications.";
    }
    else {
        msg = "Calibrating scale of anchoring parameters.";
    }
    logging::log(logging::Basic, msg);
    log_memory_usage(logging::Debug);
    
    std::vector<double> intrinsic_scales;
    
    std::vector<std::pair<std::string, Alignment>> bond_alns;
    
    auto leaves = main_execution.leaf_subproblems();
    std::vector<std::pair<std::vector<match_set_t>, std::vector<anchor_t>>> match_query_memo;
    if (cyclize_tandem_duplications) {
        // we'll save the results
        match_query_memo.resize(leaves.size());
    }
    
    for (size_t i = 0; i < leaves.size(); ++i) {
        
        logging::log(logging::Verbose, "Estimating scale for sequence " + to_string(i + 1) + " of " + to_string(leaves.size()) + ".");
        
        auto& subproblem = *leaves[i];
        
        reassign_sentinels(subproblem.graph, subproblem.tableau, 5, 6);
        SentinelTableau dummy_tableau = subproblem.tableau;
        dummy_tableau.src_sentinel = 7;
        dummy_tableau.snk_sentinel = 8;
        std::vector<match_set_t> matches = path_match_finder.find_matches(subproblem.graph, subproblem.graph,
                                                                          subproblem.tableau, dummy_tableau);
        
        // subset down to only matches on the main diagonal (retaining count for scoring)
        std::vector<match_set_t> diagonal_matches;
        diagonal_matches.reserve(matches.size());
        for (const auto& match_set : matches) {
            for (const auto& walk : match_set.walks1) {
                diagonal_matches.emplace_back();
                auto& match = diagonal_matches.back();
                match.walks1.emplace_back(walk);
                match.walks2.emplace_back(walk);
                match.count1 = match_set.count1;
                match.count2 = match_set.count2;
                match.full_length = match_set.full_length;
            }
        }
        
        ChainMerge chain_merge(subproblem.graph, subproblem.tableau);
        
        // compute the chain and scale
        std::vector<anchor_t> chain;
        double scale = anchorer.estimate_score_scale(diagonal_matches, subproblem.graph, subproblem.graph,
                                                     subproblem.tableau, subproblem.tableau,
                                                     chain_merge, chain_merge, &chain);
        
        {
            // clear the diagonal restricted matches out, we don't need them anymore
            auto dummy = std::move(diagonal_matches);
        }
        std::cerr << "got anchors\n";
        
        intrinsic_scales.push_back(scale);
        
        logging::log(logging::Debug, "Compute intrinsic scale of " + std::to_string(scale) + " for sequence " + subproblem.name);
        
        if (cyclize_tandem_duplications && !restarted_bond_alignments.get()) {
            // save the results to use in cyclizing
            match_query_memo[i].first = std::move(matches);
            match_query_memo[i].second = std::move(chain);
        }
        
        log_memory_usage(logging::Debug);
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
    
    if (cyclize_tandem_duplications && !restarted_bond_alignments.get()) {
        // we need to reanchor and find tandem duplications
        
        size_t scale_idx = 0;
        for (size_t i = 0; i < leaves.size(); ++i) {
            
            auto& subproblem = *leaves[i];
            
            PathMerge<> path_merge(subproblem.graph, subproblem.tableau);
            
            // pull the matches that we queried
            auto matches = std::move(match_query_memo[i].first);
            auto chain = std::move(match_query_memo[i].second);
            
            auto mask = generate_diagonal_mask(matches);
            logging::log(logging::Debug, "Initial mask consists of " + std::to_string(mask.size()) + " matches");
            
            StepIndex step_index;
            
            size_t bonds_identified = 0;
            
            for (size_t iter = 0; iter < max_tandem_duplication_search_rounds; ++iter) {
                
                logging::log(logging::Verbose, "Beginning round " + std::to_string(iter + 1) + " of tandem duplication detection of (at maximum) " + std::to_string(max_tandem_duplication_search_rounds) + " for sequence " + subproblem.name + ".");
                
                // get the next-best unmasked chain
                auto secondary_chain = anchorer.anchor_chain(matches, subproblem.graph, subproblem.graph,
                                                             subproblem.tableau, subproblem.tableau,
                                                             path_merge, path_merge, &mask, &intrinsic_scales[scale_idx]);
                
                // identify high-enough scoring segments
                auto bonds = bonder.identify_bonds(subproblem.graph, subproblem.graph,
                                                   subproblem.tableau, subproblem.tableau,
                                                   path_merge, path_merge,
                                                   chain, secondary_chain);
                
                bonder.deduplicate_self_bonds(bonds);
                
                log_memory_usage(logging::Debug);
                
                logging::log(logging::Verbose, "Found " + std::to_string(bonds.size()) + " tandem duplications in this round.");
                
                if (bonds.empty()) {
                    // if we didn't find any bonds this round, we're unlikely to in the future
                    break;
                }
                
                if (iter == 0) {
                    // initialize a step index
                    step_index = std::move(StepIndex(subproblem.graph));
                }
                
                // stitch the bonds into alignments
                for (auto& bond : bonds) {
                    
                    auto bond_chain = bonds_to_chain(subproblem.graph, bond);
                    
                    static const bool instrument_bond_partition = false;
                    if (instrument_bond_partition) {
                        auto copy_chain = bond_chain; // because it will steal the anchors
                        auto partitioned = partitioner.partition_anchors(copy_chain, subproblem.graph, subproblem.graph,
                                                                         subproblem.tableau, subproblem.tableau,
                                                                         path_merge, path_merge, true);
                        std::cerr << "partitioned bond chain:\n";
                        for (size_t i = 0; i < partitioned.size(); ++i) {
                            std::cerr << '}' << '\t' << partitioned[i].front().walk1.front() << '\t' << partitioned[i].back().walk1.back() << '\t' << partitioned[i].front().walk2.front() << '\t' << partitioned[i].back().walk2.back() << '\n';
                        }
                    }
                    
                    // get the alignment with node IDs
                    bond_alns.emplace_back(subproblem.graph.path_name(0),
                                           stitcher.internal_stitch(bond_chain, subproblem.graph, path_merge));
                    
                    if (!bonds_prefix.empty()) {
                        output_bond_alignment(bond_alns.back().second, subproblem.graph, 0, bonds_identified);
                    }
                    
                    // convert it to an alignment with path positions
                    for (auto& aln_pair : bond_alns.back().second) {
                        if (aln_pair.node_id1 != AlignedPair::gap) {
                            aln_pair.node_id1 = step_index.path_steps(aln_pair.node_id1).front().second;
                        }
                        if (aln_pair.node_id2 != AlignedPair::gap) {
                            aln_pair.node_id2 = step_index.path_steps(aln_pair.node_id2).front().second;
                        }
                    }
                    
                    ++bonds_identified;
                }
                
                if (iter != max_tandem_duplication_search_rounds) {
                    // mask out anchors that overlap this chain's matches
                    update_mask(matches, secondary_chain, mask, true);
                    
                    logging::log(logging::Debug, "Updated mask consists of " + std::to_string(mask.size()) + " matches");
                }
            }
            ++scale_idx;
        }
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

std::string Core::subproblem_bond_file_name() const {
    return subproblems_prefix + "_bonds.txt";
}

std::string Core::subproblem_file_name(const Subproblem& subproblem) const {
    return subproblems_prefix + "_" + to_hex(main_execution.subproblem_hash(subproblem)) + ".gfa";
}


std::string Core::get_subpath_name(const std::string& path_name, size_t begin, size_t end) const {
    return path_name + ":" + std::to_string(begin) + "-" + std::to_string(end);
}

std::tuple<std::string, size_t, size_t> Core::parse_subpath_name(const std::string& subpath_name)  const {
    size_t sep = subpath_name.find_last_of(':');
    auto it = std::find(subpath_name.begin() + sep + 1, subpath_name.end(), '-');
    return std::tuple<std::string, size_t, size_t>(subpath_name.substr(0, sep),
                                                   parse_int(std::string(subpath_name.begin() + sep + 1, it)),
                                                   parse_int(std::string(it + 1, subpath_name.end())));
    
    ;
}

void Core::emit_subproblem(const Subproblem& subproblem) const {
    
    auto gfa_file_name = subproblem_file_name(subproblem);
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
    
    if (write_header) {
        info_out << "filename\tsequences\n";
    }
    auto sequences = main_execution.leaf_descendents(subproblem);
    sort(sequences.begin(), sequences.end());
    info_out << gfa_file_name << '\t' << join(sequences, ",") << '\n';
    
    write_gfa(subproblem.graph, subproblem.tableau, gfa_out);
}

void Core::emit_subalignment() const {
    
    // TODO: this is fragile in that we need to separately use the same order for
    // the children here and in the alignment routine
    auto problem_ptrs = main_execution.current();
    const auto& subproblem = *get<0>(problem_ptrs);
    const auto& child1 = *get<1>(problem_ptrs);
    const auto& child2 = *get<2>(problem_ptrs);
    const auto& graph1 = child1.graph;
    const auto& graph2 = child2.graph;
    
    ofstream out(subalignments_filepath, ios_base::app);
    if (!out) {
        throw std::runtime_error("Failed to write to subalignment file " + subalignments_filepath);
    }
    
    out << "# sequence set 1\n";
    for (const auto& seq_name : main_execution.leaf_descendents(child1)) {
        out << seq_name << '\n';
    }
    out << "# sequence set 2\n";
    for (const auto& seq_name : main_execution.leaf_descendents(child2)) {
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

void Core::emit_restart_bonds(const std::vector<std::pair<std::string, Alignment>>& bond_alignments) const {
    
    ofstream out(subproblem_bond_file_name());
    if (!out) {
        throw std::runtime_error("Couldn't write subproblem bonds to file '" + subproblem_bond_file_name() + "'.");
    }
    
    for (size_t i = 0; i < bond_alignments.size(); ++i) {
        out << '#' << bond_alignments[i].first << '\n';
        for (const auto& aln_pair : bond_alignments[i].second) {
            out << (int64_t) aln_pair.node_id1 << '\t' << (int64_t) aln_pair.node_id2 << '\n';
        }
    }
}

void Core::restart_bonds() {
    
    if (!cyclize_tandem_duplications) {
        return;
    }
    
    restarted_bond_alignments.reset(new std::vector<std::pair<std::string, Alignment>>());
    auto& restart_alignments = *restarted_bond_alignments;
    
    ifstream in(subproblem_bond_file_name());
    if (!in) {
        throw std::runtime_error("Couldn't open tandem duplication bonds to restart from file " + subproblem_bond_file_name() + ".");
    }
    
    string line;
    while (in) {
        getline(in, line);
        if (line.empty()) {
            continue;
        }
        if (line.front() == '#') {
            restart_alignments.emplace_back();
            restart_alignments.back().first = std::move(line.substr(1, line.size()));
        }
        else {
            auto tokens = tokenize(line);
            assert(tokens.size() == 2);
            restart_alignments.back().second.emplace_back(parse_int(tokens[0]), parse_int(tokens[1]));
        }
    }
}

void Core::output_pairwise_alignments(bool cyclic) const {
    
    assert(!induced_pairwise_prefix.empty());
    
    const auto& graph = main_execution.final_subproblem().graph;
    
    for (uint64_t path_id1 = 0; path_id1 < graph.path_size(); ++path_id1) {
        for (uint64_t path_id2 = path_id1 + 1; path_id2 < graph.path_size(); ++path_id2) {
            
            auto path_name1 = graph.path_name(path_id1);
            auto path_name2 = graph.path_name(path_id2);
            // get rid of slashes in path names that look like subdirectories
            std::replace(path_name1.begin(), path_name1.end(), '/', '_');
            std::replace(path_name2.begin(), path_name2.end(), '/', '_');
            
            auto out_filename = induced_pairwise_prefix + "_" + path_name1 + "_" + path_name2 + (cyclic ? ".maf" : ".txt");
            ofstream out_file(out_filename);
            if (!out_file) {
                throw runtime_error("could not write to induced pairwise alignment file " + out_filename + "\n");
                
            }
            if (cyclic) {
                output_maf(out_file, induced_cyclic_pairwise_alignment(graph, path_id1, path_id2),
                           graph, path_id1, path_id2);
            }
            else {
                out_file << explicit_cigar(induced_pairwise_alignment(graph, path_id1, path_id2),
                                           path_to_string(graph, graph.path(path_id1)),
                                           path_to_string(graph, graph.path(path_id2))) << '\n';
            }
        }
    }
}

void Core::apply_bonds(std::vector<std::pair<std::string, Alignment>>& bond_alignments) {
    
    if (bond_alignments.empty()) {
        return;
    }
    
    logging::log(logging::Basic, "Cyclizing the final graph.");
    
    logging::log(logging::Verbose, "Merging " + std::to_string(bond_alignments.size()) + " tandem duplications.");
    if (logging::level >= logging::Debug) {
        size_t aln_len = 0;
        for (const auto& aln : bond_alignments) {
            for (const auto& aln_pair : aln.second) {
                aln_len += (aln_pair.node_id1 != AlignedPair::gap && aln_pair.node_id2 != AlignedPair::gap);
            }
        }
        logging::log(logging::Debug, "Tandem duplication alignments total " + std::to_string(aln_len) + " aligned bases.");
    }
    
    auto& root_subproblem = main_execution.final_subproblem();
        
    // convert from path positions to node IDs
    std::vector<Alignment> alignments_to_fuse;
    alignments_to_fuse.reserve(bond_alignments.size());
    for (auto& bond_aln : bond_alignments) {
        uint64_t path_id = root_subproblem.graph.path_id(bond_aln.first);
        for (auto& aln_pair : bond_aln.second) {
            if (aln_pair.node_id1 != AlignedPair::gap) {
                aln_pair.node_id1 = root_subproblem.graph.path(path_id)[aln_pair.node_id1];
            }
            if (aln_pair.node_id2 != AlignedPair::gap) {
                aln_pair.node_id2 = root_subproblem.graph.path(path_id)[aln_pair.node_id2];
            }
        }
        alignments_to_fuse.emplace_back(std::move(bond_aln.second));
    }
    
    // merge along all of the bond alignments
    SentinelTableau cyclized_tableau;
    Alignment cyclized_alignment;
    BaseGraph cyclized = internal_fuse(root_subproblem.graph, alignments_to_fuse,
                                       &root_subproblem.tableau, &cyclized_tableau,
                                       &root_subproblem.alignment, &cyclized_alignment);
    simplify_bubbles(cyclized, cyclized_tableau);
    
    logging::log(logging::Debug, "Cyclized graph reduces from " + std::to_string(root_subproblem.graph.node_size()) + " to " + std::to_string(cyclized.node_size()) + " nodes after merging.");
        
    root_subproblem.graph = std::move(cyclized);
    root_subproblem.tableau = cyclized_tableau;
    // FIXME: this currently breaks under bubble merging
    //root_subproblem.alignment = std::move(cyclized_alignment);
    root_subproblem.alignment.clear();
        
    // FIXME: not a real fix -- do something more principled
    //partitioner.minimum_segment_score = 0.0;
    
    polish_cyclized_graph(root_subproblem);
}

void Core::polish_cyclized_graph(Subproblem& subproblem) const {
        
    auto inconsistencies = inconsistency_identifier.identify_inconsistencies(subproblem.graph, subproblem.tableau);
    
    StepIndex step_index(subproblem.graph);
    
    static const bool instrument_inconsistencies = false;
    if (instrument_inconsistencies) {
        static const bool all_locations = false;
        std::cerr << "found " << inconsistencies.size() << " potential inconsistencies\n";
        for (const auto& bounds : inconsistencies) {
            std::cerr << '}' << '\t' << subproblem.graph.path_name(step_index.path_steps(bounds.first).front().first) << '\t' << step_index.path_steps(bounds.first).front().second << '\t' << subproblem.graph.path_name(step_index.path_steps(bounds.second).front().first) << '\t' << step_index.path_steps(bounds.second).front().second << '\n';
            if (all_locations) {
                for (auto step : step_index.path_steps(bounds.first)) {
                    std::cerr << '\t' << 'l' << '\t' << subproblem.graph.path_name(step.first) << '\t' << step.second << '\n';
                }
                for (auto step : step_index.path_steps(bounds.second)) {
                    std::cerr << '\t' << 'r' << '\t' << subproblem.graph.path_name(step.first) << '\t' << step.second << '\n';
                }
            }
        }
    }
    
    reassign_sentinels(subproblem.graph, subproblem.tableau, 5, 6);
    
    // the path ESA doesn't require the graph to actually have the sentinel values as node labels,
    // so we just make a fictitious tableau here instead of copying it
    // TODO: poor segmentation
    SentinelTableau dummy_tableau = subproblem.tableau;
    dummy_tableau.src_sentinel = 7;
    dummy_tableau.snk_sentinel = 8;
    
    // query self-matches in the full graph to get realistic global frequencies for anchors
    std::vector<match_set_t> full_match_set = path_match_finder.find_matches(subproblem.graph, subproblem.graph,
                                                                             subproblem.tableau, dummy_tableau);
    
    // this will synthesize matches for the inconsistency subproblems based on the global frequencies
    InducedMatchFinder induced_match_finder(subproblem.graph, full_match_set, inconsistencies, step_index);
    
    std::vector<Subproblem> realigned;
    
    for (size_t i = 0; i < inconsistencies.size(); ++i) {
        
        auto inconsistency = inconsistencies[i];
                
        // organize the steps by their path
        std::unordered_map<uint64_t, std::pair<std::vector<size_t>, std::vector<size_t>>> path_locations;
        for (auto step : step_index.path_steps(inconsistency.first)) {
            path_locations[step.first].first.push_back(step.second);
        }
        for (auto step : step_index.path_steps(inconsistency.second)) {
            path_locations[step.first].second.push_back(step.second);
        }
        
        std::vector<uint64_t> path_ids;
        // make sure that the occurrences of this interval are properly paired
        for (auto it = path_locations.begin(); it != path_locations.end(); ++it) {
            path_ids.push_back(it->first);
            if (!std::is_sorted(it->second.first.begin(), it->second.first.end())) {
                std::sort(it->second.first.begin(), it->second.first.end());
            }
            if (!std::is_sorted(it->second.second.begin(), it->second.second.end())) {
                std::sort(it->second.second.begin(), it->second.first.end());
            }
        }
        // to get an implementation-independent ordering over path IDs
        std::sort(path_ids.begin(), path_ids.end());
        
        // records of (path id, begin, end)
        std::vector<std::tuple<uint64_t, size_t, size_t>> subpath_intervals;
        // names and sequences
        std::vector<std::pair<std::string, std::string>> subpaths;
        
        // convert the intervals into the form we want
        for (auto path_id : path_ids) {
            const auto& locations = path_locations[path_id];
            if (locations.first.size() != locations.second.size()) {
                throw std::runtime_error("Path starts or ends in the middle of a cycle realignment interval");
            }
            for (size_t i = 0; i < locations.first.size(); ++i) {
                subpath_intervals.emplace_back(path_id, locations.first[i], locations.second[i]);
                subpaths.emplace_back();
                subpaths.back().first = get_subpath_name(subproblem.graph.path_name(path_id), locations.first[i], locations.second[i]);
                for (size_t j = locations.first[i], n = locations.second[i]; j <= n; ++j) {
                    // note: we have to decode because the BaseGraph convenience constructor assumes that the sequences are raw
                    subpaths.back().second.push_back(decode_base(subproblem.graph.label(subproblem.graph.path(path_id)[j])));
                }
            }
        }
        
        auto expanded_tree = make_copy_expanded_tree(subpath_intervals, subpaths);
        
        Execution realignment(std::move(subpaths), std::move(expanded_tree));
        
        do_execution(realignment, induced_match_finder.component_view(i), false);
        
        realigned.emplace_back(std::move(realignment.final_subproblem()));
    }
    
    integrate_polished_subgraphs(subproblem, realigned);
}

Tree Core::make_copy_expanded_tree(const std::vector<std::tuple<uint64_t, size_t, size_t>>& subpath_intervals,
                                   const std::vector<std::pair<std::string, std::string>>& subpaths) const {
    
    static const bool debug = false;
    if (debug) {
        std::cerr << "making expanded tree for subpaths:\n";
        assert(subpath_intervals.size() == subpaths.size());
        for (size_t i = 0; i < subpath_intervals.size(); ++i) {
            std::cerr << subpaths[i].first << '\t' << std::get<0>(subpath_intervals[i]) << '\t' << std::get<1>(subpath_intervals[i]) << '\t' << std::get<2>(subpath_intervals[i]) << '\n';
        }
    }
    
    const auto& tree = main_execution.get_tree();
    
    std::unordered_map<std::string, std::vector<std::string>> copies;
    {
        // make sure that we add these in order by sequence interval
        auto indexes = range_vector(subpath_intervals.size());
        std::sort(indexes.begin(), indexes.end(), [&](size_t i, size_t j) {
            return subpath_intervals[i] < subpath_intervals[j];
        });
        for (auto i : indexes) {
            const auto& subpath = subpaths[i];
            copies[std::get<0>(parse_subpath_name(subpath.first))].push_back(subpath.first);
        }
    }
    
    // identify subtrees that have the same copy count across all paths
    std::vector<uint64_t> subtree_copy_count(tree.node_size(), 0);
    for (const auto& copy_record : copies) {
        subtree_copy_count[tree.get_id(copy_record.first)] = copy_record.second.size();
    }
    for (auto node_id : tree.postorder()) {
        if (tree.is_leaf(node_id)) {
            continue;
        }
        
        // check if the copy count of all the observed children is consistent
        uint64_t last_count = -2; // sentinel for unobserved
        for (auto child_id : tree.get_children(node_id)) {
            if (subtree_copy_count[child_id] == -1 ||
                (last_count != -2 && subtree_copy_count[child_id] != last_count)) {
                last_count = -1; // sentinel for inconsistent
                break;
            }
            if (subtree_copy_count[child_id] != 0) {
                last_count = subtree_copy_count[child_id];
            }
        }
        if (last_count != -2) {
            subtree_copy_count[node_id] = last_count;
        }
    }
    if (debug) {
        std::cerr << "subtree copy counts:\n";
        for (uint64_t node_id = 0; node_id < tree.node_size(); ++node_id) {
            std::cerr << '\t' << node_id << '\t' << subtree_copy_count[node_id] << '\t' << (int64_t) tree.get_parent(node_id) << '\t' << tree.label(node_id) << '\n';
        }
    }
    
    // now we construct a Newick string for this tree
    // TODO: repetitive with Tree::to_newick
    // TODO: it feels hokey to construct the tree using a newick string rather than directly
    std::stringstream strm;
    // records of (node ID, which copy, children, next child)
    std::vector<std::tuple<uint64_t, size_t, std::vector<std::pair<uint64_t, size_t>>, size_t>> stack;
    
    if (subtree_copy_count[tree.get_root()] == 0) {
        throw std::runtime_error("Root is not included in induced subpath tree");
    }
    
    stack.emplace_back();
    if (subtree_copy_count[tree.get_root()] == -1) {
        // no consistent copy count;
        std::get<0>(stack.back()) = tree.get_root();
        std::get<1>(stack.back()) = -1;
        for (auto child_id : tree.get_children(tree.get_root())) {
            std::get<2>(stack.back()).emplace_back(child_id, -1);
        }
    }
    else {
        // has a consistent copy count at the root (must be the first encountered consistent copy count)
        // make a virtual node
        std::get<0>(stack.back()) = -1;
        std::get<1>(stack.back()) = -1;
        for (size_t i = 0; i < subtree_copy_count[tree.get_root()]; ++i) {
            std::get<2>(stack.back()).emplace_back(tree.get_root(), i);
        }
    }
    std::get<3>(stack.back()) = 0;
    
    while (!stack.empty()) {
        
        auto& top = stack.back();
        
        if (std::get<3>(top) == std::get<2>(top).size()) {
            // we've traversed the last of this node's edges
            if (!std::get<2>(top).empty()) {
                // it has children, close out their list
                strm << ')';
            }
            if (std::get<0>(top) != -1 && tree.is_leaf(std::get<0>(top))) {
                if (std::get<1>(top) == -1) {
                    throw std::runtime_error("Leaf of induced subpath tree was not marked as having consistent count");
                }
                // output the label corresponding to this copy of the leaf
                strm << '"' << copies.at(tree.label(std::get<0>(top)))[std::get<1>(top)] << '"';
            }
            
            // put virtual nodes at distance 0
            double dist = std::get<0>(top) == -1 ? 0.0 : tree.distance(std::get<0>(top));
            if (dist != std::numeric_limits<double>::max()) {
                strm << ':' << dist;
            }
            
            stack.pop_back();
            continue;
        }
        
        if (std::get<3>(top) == 0) {
            // this node has downward edges, but we haven't traversed any yet
            strm << '(';
        }
        else {
            // separate the edges
            strm << ',';
        }
        
        uint64_t next_id;
        size_t which_copy;
        std::tie(next_id, which_copy) = std::get<2>(top)[std::get<3>(top)++];
        stack.emplace_back();
        auto& next = stack.back();
        if (which_copy == -1 && subtree_copy_count[next_id] != -1) {
            // this is the first instance of a copy consistent subtree we've encountered, make a virtual node
            // to house the copies
            std::get<0>(next) = -1;
            std::get<1>(next) = -1;
            for (size_t i = 0; i < subtree_copy_count[next_id]; ++i) {
                if (subtree_copy_count[next_id] != 0) {
                    std::get<2>(next).emplace_back(next_id, i);
                }
            }
        }
        else {
            // we either haven't encountered a copy consistent subtree, or we are already in one--pass that status down
            std::get<0>(next) = next_id;
            std::get<1>(next) = which_copy;
            for (auto child_id : tree.get_children(next_id)) {
                if (subtree_copy_count[next_id] != 0) {
                    std::get<2>(next).emplace_back(child_id, which_copy);
                }
            }
        }
        std::get<3>(next) = 0;
    }
    
    strm << ';';
    auto newick = strm.str();
    
    if (debug) {
        std::cerr << "input newick string:\n" << newick << '\n';
    }
    
    Tree expanded(newick);
    
    expanded.compact();
    
    expanded.binarize();
    
    if (debug) {
        std::cerr << "processed newick string:\n" << expanded.to_newick() << '\n';
    }
    
    return expanded;
}

void Core::integrate_polished_subgraphs(Subproblem& root, const std::vector<Subproblem>& realigned_graphs) const {
        
    static const bool debug = false;
    
    for (const auto& realigned : realigned_graphs) {
        
        if (debug) {
            std::cerr << "adding new realigned graph with " << realigned.graph.node_size() << " nodes to parent graph with " << root.graph.node_size() << " nodes\n";
        }
        
        // add nodes
        std::vector<uint64_t> inject_trans(realigned.graph.node_size(), -1);
        for (uint64_t node_id = 0; node_id < realigned.graph.node_size(); ++node_id) {
            if (node_id != realigned.tableau.src_id && node_id != realigned.tableau.snk_id) {
                inject_trans[node_id] = root.graph.add_node(realigned.graph.label(node_id));
            }
        }
        // add edges
        for (uint64_t node_id = 0; node_id < realigned.graph.node_size(); ++node_id) {
            if (node_id == realigned.tableau.src_id || node_id == realigned.tableau.snk_id) {
                continue;
            }
            for (auto next_id : realigned.graph.next(node_id)) {
                if (next_id == realigned.tableau.src_id || next_id == realigned.tableau.snk_id) {
                    continue;
                }
                root.graph.add_edge(inject_trans[node_id], inject_trans[next_id]);
            }
        }
        // rewire the paths
        std::unordered_set<std::pair<uint64_t, uint64_t>> path_adjacencies;
        for (uint64_t path_id = 0; path_id < realigned.graph.path_size(); ++path_id) {
            
            if (debug) {
                std::cerr << "path " << realigned.graph.path_name(path_id) << " projects to nodes between " << inject_trans[realigned.graph.path(path_id).front()] << " and " << inject_trans[realigned.graph.path(path_id).back()] << '\n';
            }
            
            std::string path_name;
            size_t begin, end;
            std::tie(path_name, begin, end) = parse_subpath_name(realigned.graph.path_name(path_id));
            
            if (begin == end) {
                continue;
            }
            
            uint64_t root_path_id = root.graph.path_id(path_name);
            
            if (debug) {
                std::cerr << "current neighbors are " << root.graph.path(root_path_id)[begin - 1] << " and " << root.graph.path(root_path_id)[end + 1] << '\n';
            }
            
            // add edges to the next node on the path
            uint64_t prev_id, next_id;
            if (begin == 0) {
                prev_id = root.tableau.src_id;
            }
            else {
                prev_id = root.graph.path(root_path_id)[begin - 1];
            }
            if (end + 1 == root.graph.path(root_path_id).size()) {
                next_id = root.tableau.snk_id;
            }
            else {
                next_id = root.graph.path(root_path_id)[end + 1];
            }
            
            if (path_adjacencies.emplace(prev_id, inject_trans[realigned.graph.path(path_id).front()]).second) {
                root.graph.add_edge(prev_id, inject_trans[realigned.graph.path(path_id).front()]);
                if (debug) {
                    std::cerr << "add backward connection edge " << prev_id << " -> " << inject_trans[realigned.graph.path(path_id).front()] << '\n';
                }
            }
            if (path_adjacencies.emplace(inject_trans[realigned.graph.path(path_id).back()], next_id).second) {
                root.graph.add_edge(inject_trans[realigned.graph.path(path_id).back()], next_id);
                if (debug) {
                    std::cerr << "add forward connection edge " << inject_trans[realigned.graph.path(path_id).back()] << " -> " << next_id << '\n';
                }
            }
            
            // swap over to the new node IDs
            std::vector<uint64_t> new_ids;
            new_ids.reserve(end - begin);
            for (auto node_id : realigned.graph.path(path_id)) {
                new_ids.emplace_back(inject_trans[node_id]);
            }
            root.graph.reassign_subpath(root_path_id, begin, new_ids);
        }
    }
    
    // remove the nodes corresponding to the old paths
    purge_uncovered_nodes(root.graph, root.tableau);
}

void Core::restart() {
    
    std::function<std::string(const Subproblem&)> get_file_name = [&](const Subproblem& subproblem) -> std::string {
        return subproblem_file_name(subproblem);
    };
    main_execution.restart(get_file_name, preserve_subproblems || !skip_calibration, preserve_subproblems);
    
    if (cyclize_tandem_duplications) {
        restart_bonds();
    }
}

const Subproblem& Core::root_subproblem() const {
    return main_execution.final_subproblem();
}

const Subproblem& Core::leaf_subproblem(const std::string& name) const {
    return main_execution.leaf_subproblem(name);
}

//const Core::Subproblem& Core::subproblem_covering(const std::vector<string>& names) const {
//    // TODO: i'm sure there's an O(n) algorithm for this, but hopefully this O(n^2) is
//    // good enough
//
//    if (names.empty()){
//        throw runtime_error("Cannor query an MSA subproblem with an empty set of sequences");
//    }
//
//    // walk up to the root on an arbitrary node
//    vector<uint64_t> path_to_root;
//    {
//        uint64_t node_id = tree.get_id(names.front());
//        path_to_root.push_back(node_id);
//        while (node_id != tree.get_root()) {
//            node_id = tree.get_parent(node_id);
//            path_to_root.push_back(node_id);
//        }
//    }
//    // we'll be removing from the bottom of this path
//    reverse(path_to_root.begin(), path_to_root.end());
//
//    // walk up to the root on each other node
//    for (size_t i = 1; i < names.size(); ++i) {
//        // TODO: repetitive
//
//        unordered_set<uint64_t> next_path_to_root;
//
//        uint64_t node_id = tree.get_id(names[i]);
//        next_path_to_root.insert(node_id);
//        while (node_id != tree.get_root()) {
//            node_id = tree.get_parent(node_id);
//            next_path_to_root.insert(node_id);
//        }
//
//        // remove from the first path until we find a shared node
//        while (!next_path_to_root.count(path_to_root.back())) {
//            path_to_root.pop_back();
//        }
//    }
//
//    return subproblems[path_to_root.back()];
//}


}
