#ifndef centrolign_core_hpp
#define centrolign_core_hpp

#include <vector>
#include <string>
#include <memory>
#include <cmath>

#include "centrolign/anchorer.hpp"
#include "centrolign/stitcher.hpp"
#include "centrolign/tree.hpp"
#include "centrolign/simplifier.hpp"
#include "centrolign/match_finder.hpp"
#include "centrolign/partitioner.hpp"
#include "centrolign/score_function.hpp"
#include "centrolign/utility.hpp"
#include "centrolign/bonder.hpp"
#include "centrolign/inconsistency_identifier.hpp"
#include "centrolign/execution.hpp"
#include "centrolign/path_merge.hpp"
#include "centrolign/packed_path_merge.hpp"
#include "centrolign/chain_merge.hpp"
#include "centrolign/fuse.hpp"

namespace centrolign {

/*
 * The core runner for centrolign
 */
class Core {
public:
    
    // parse files to construct core (either may be - for stdin)
    Core(const std::string& fasta_file, const std::string& tree_file);
    
    // construct core from already-parsed inputs (consumes the inputs)
    Core(std::vector<std::pair<std::string, std::string>>&& names_and_sequences,
         Tree&& tree);
    
    Core() = default;
    ~Core() = default;
    
    // trigger the MSA
    void execute();
    
    /*
     * Configurable submodules
     */
    
    // a function for scoring anchors
    ScoreFunction score_function;
    // queries matches between the input graphs
    PathMatchFinder path_match_finder;
    // makes a chain of alignment anchors using matches
    Anchorer anchorer;
    // partitions anchor chains into well-anchored and poorly-anchored portions
    Partitioner partitioner;
    // aligns in between anchors to them stitch into a base-level alignment
    Stitcher stitcher;
    // identifies sequences to bond together as cycles (if cyclizing)
    Bonder bonder;
    // flags graph regions with potential cyclization-induced artifacts for normalization
    InconsistencyIdentifier inconsistency_identifier;
    
    
    // don't calibrate the scale of the scoring function before executing
    bool skip_calibration = false;
    
    // preserve subproblems whose parent problems have been completed
    bool preserve_subproblems = true;
    
    // merge tandem duplications into cycles in the final graph
    bool cyclize_tandem_duplications = false;
    
    // number of threads for parallel sections
    uint64_t threads = 1;
    
    // switch to slower, more memory-efficient data structures when the (graph size * num sequences) hits this amount
    size_t memory_restraint_size = 1 << 30;
    
    // the maximum number of overlapping duplications that will be found
    size_t max_tandem_duplication_search_rounds = 3;
    
    // if non-empty, prefix to give GFA output for all suproblems
    std::string subproblems_prefix;
    
    // if non-empty, file to write suproblem alignments to
    std::string subalignments_filepath;
    
    // if non-empty, write a file for each induced pairwise alignment after completion
    std::string induced_pairwise_prefix;
    
    // if non-empty, write a file for each bond alignment after identifying them
    std::string bonds_prefix;
    
    // load alignments from the prefix and start where they left off
    void restart();
    
    // the root subproblem
    const Subproblem& root_subproblem() const;
    
    // the leaf subproblem that corresponds to a sequences
    const Subproblem& leaf_subproblem(const std::string& name) const;
    
protected:
    
    void init(std::vector<std::pair<std::string, std::string>>&& names_and_sequences,
              Tree&& tree_in);
    
    template<class MFinder>
    void do_execution(Execution& execution, const MFinder& match_finder, bool is_main_execution) const;
    
    std::string subproblem_file_name(const Subproblem& subproblem) const;
    
    std::string subproblem_info_file_name() const;
    
    std::string subproblem_bond_file_name() const;
    
    std::string get_subpath_name(const std::string& path_name, size_t begin, size_t end) const;
    
    std::tuple<std::string, size_t, size_t> parse_subpath_name(const std::string& subpath_name) const;
    
    void emit_subproblem(const Subproblem& subproblem) const;
    
    void emit_subalignment() const;
    
    void emit_restart_bonds(const std::vector<std::pair<std::string, Alignment>>& bond_alignments) const;
    
    void restart_bonds();
    
    std::vector<match_set_t> get_matches(Subproblem& subproblem1, Subproblem& subproblem2,
                                         bool suppress_verbose_logging) const;
    
    std::vector<match_set_t> query_matches(ExpandedGraph& expanded1,
                                           ExpandedGraph& expanded2) const;
    
    template<class XMerge>
    Alignment align(std::vector<match_set_t>& matches,
                    const Subproblem& subproblem1, const Subproblem& subproblem2,
                    XMerge& xmerge1, XMerge& xmerge2, bool is_main_execution) const;
    
    std::unordered_set<std::tuple<size_t, size_t, size_t>> generate_diagonal_mask(const std::vector<match_set_t>& matches) const;
    
    void update_mask(const std::vector<match_set_t>& matches, const std::vector<anchor_t>& chain,
                     std::unordered_set<std::tuple<size_t, size_t, size_t>>& masked_matches, bool mask_reciprocal = false) const;
    
    template<class BGraph>
    std::vector<anchor_t> bonds_to_chain(const BGraph& graph, const bond_interval_t& bond_interval) const;
    
    void apply_bonds(std::vector<std::pair<std::string, Alignment>>& bond_alignments);
    
    void polish_cyclized_graph(Subproblem& subproblem) const;
    
    Tree make_copy_expanded_tree(const std::vector<std::tuple<uint64_t, size_t, size_t>>& subpath_intervals,
                                 const std::vector<std::pair<std::string, std::string>>& subpaths) const;
    
    void integrate_polished_subgraphs(Subproblem& subproblem, const std::vector<Subproblem>& realigned_graphs) const;
    
    void output_pairwise_alignments(bool cyclic) const;
    
    template<class BGraph>
    void output_bond_alignment(const Alignment& bond_alignment, const BGraph& graph, uint64_t path_id, size_t bond_number) const;
    
    // learn the intrinsic scale of the anchor scoring function on these sequences
    std::vector<std::pair<std::string, Alignment>> calibrate_anchor_scores_and_identify_bonds();
        
    // the primary MSA problem
    Execution main_execution;
    
    // TODO: ugly
    std::unique_ptr<std::vector<std::pair<std::string, Alignment>>> restarted_bond_alignments;
};





/*
 * Template and inline implementations
 */

template<class XMerge>
Alignment Core::align(std::vector<match_set_t>& matches,
                      const Subproblem& subproblem1, const Subproblem& subproblem2,
                      XMerge& xmerge1, XMerge& xmerge2, bool is_main_execution) const {
    
    if (logging::level >= logging::Debug) {
        size_t merge_size = xmerge1.memory_size() + xmerge2.memory_size();
        logging::log(logging::Debug, "Merge data structures of type " + std::string(typeid(XMerge).name()) + " are consuming " + format_memory_usage(merge_size) + " of memory.");
        log_memory_usage(logging::Debug);
    }
    
    // get the best anchor chain
    bool restrain_memory = (subproblem1.graph.path_size() * subproblem2.graph.path_size() * anchorer.max_num_match_pairs * log2(anchorer.max_num_match_pairs) > memory_restraint_size);
    auto anchors = anchorer.anchor_chain(matches, subproblem1.graph, subproblem2.graph,
                                         subproblem1.tableau, subproblem2.tableau,
                                         xmerge1, xmerge2, restrain_memory);
    
    log_memory_usage(logging::Debug);
    
    static const bool output_anchors = false;
    if (output_anchors) {
        std::cerr << "outputting anchors\n";
        for (const auto& a : anchors) {
            std::cout << a.walk1.front() << '\t' << a.walk2.front() << '\t' << a.walk1.size() << '\t' << 0 << '\n';
        }
        
        std::cerr << "making initial mask\n";
        
        std::unordered_set<std::tuple<size_t, size_t, size_t>> mask;
        update_mask(matches, anchors, mask);
        
        size_t num_reanchors = 1;
        bool mask_reciprocal = true;
        for (size_t i = 1; i <= num_reanchors; ++i) {
            
            auto anchors_secondary = anchorer.anchor_chain(matches, subproblem1.graph, subproblem2.graph,
                                                           subproblem1.tableau, subproblem2.tableau,
                                                           xmerge1, xmerge2, false, &mask);
            update_mask(matches, anchors_secondary, mask, mask_reciprocal);
            
            for (const auto& a : anchors_secondary) {
                std::cout << a.walk1.front() << '\t' << a.walk2.front() << '\t' << a.walk1.size() << '\t' << i << '\n';
            }
            
            auto bonds = bonder.identify_bonds(subproblem1.graph, subproblem2.graph,
                                               subproblem1.tableau, subproblem2.tableau,
                                               xmerge1, xmerge2, anchors, anchors_secondary);
        }
        
        exit(0);
    }
        
    // partition the anchor chain into good and bad segments
    auto anchor_segments = partitioner.partition_anchors(anchors, subproblem1.graph, subproblem2.graph,
                                                         subproblem1.tableau, subproblem2.tableau,
                                                         xmerge1, xmerge2,
                                                         !is_main_execution); // assume significant boundaries for fill-in problems
    
    log_memory_usage(logging::Debug);
    
    logging::log(logging::Verbose, "Stitching anchors into alignment.");
    
    for (auto& anchor_segment : anchor_segments) {
        stitcher.despecify_indel_breakpoints(anchor_segment);
    }
    
    // form a base-level alignment
    Alignment alignment = stitcher.stitch(anchor_segments, subproblem1.graph, subproblem2.graph,
                                          subproblem1.tableau, subproblem2.tableau,
                                          xmerge1, xmerge2);
    
    return alignment;
}

template<class MFinder>
void Core::do_execution(Execution& execution, const MFinder& match_finder, bool is_main_execution) const {
    
    // reduce the logging for non-main executions
    // TODO: very ugly
    logging::LoggingLevel current_log_level = logging::level;
    if (!is_main_execution) {
        if (current_log_level != logging::Silent && logging::level != logging::Debug) {
            logging::level = logging::Minimal;
        }
    }
    
    while (!execution.finished()) {
        
        auto problem_ptrs = execution.next();
        
        auto& next_problem = *std::get<0>(problem_ptrs);
        
        if (next_problem.complete) {
            logging::log(logging::Verbose, "Problem already finished from restarted run.");
            continue;
        }
        
        if (logging::level >= logging::Debug) {
            logging::log(logging::Debug, "In-memory graphs and alignments are occupying " + format_memory_usage(execution.memory_size()) + ".");
            logging::log(logging::Debug, "Current memory use is " + format_memory_usage(current_memory_usage()));
        }
        
        auto& subproblem1 = *std::get<1>(problem_ptrs);
        auto& subproblem2 = *std::get<2>(problem_ptrs);
        
        reassign_sentinels(subproblem1.graph, subproblem1.tableau, 5, 6);
        reassign_sentinels(subproblem2.graph, subproblem2.tableau, 7, 8);
        auto matches = match_finder.find_matches(subproblem1.graph, subproblem2.graph,
                                                 subproblem1.tableau, subproblem2.tableau);
        
        log_memory_usage(logging::Debug);
        
        logging::log(logging::Verbose, "Computing reachability.");
        
        if (anchorer.chaining_algorithm == Anchorer::SparseAffine) {
            // use all paths for reachability to get better distance estimates
            
//#define __FAST_BUILD
#ifdef __FAST_BUILD
            PathMerge<size_t, size_t> path_merge1(subproblem1.graph, subproblem1.tableau);
            PathMerge<size_t, size_t> path_merge2(subproblem2.graph, subproblem2.tableau);
            next_problem.alignment = std::move(align(matches, subproblem1, subproblem2,
                                                     path_merge1, path_merge2, is_main_execution));
#else
            size_t max_nodes = std::max(subproblem1.graph.node_size(), subproblem2.graph.node_size());
            size_t max_paths = std::max(subproblem1.graph.path_size(), subproblem2.graph.path_size());
            size_t total_size = (subproblem1.graph.node_size() * subproblem1.graph.path_size()
                                 + subproblem2.graph.node_size() + subproblem2.graph.path_size());
            if (total_size > memory_restraint_size) {
                #define _gen_packed_path_merge(UIntSize, UIntChain) \
                    PackedPathMerge<UIntSize, UIntChain, 2048, 127> path_merge1(subproblem1.graph, subproblem1.tableau); \
                    PackedPathMerge<UIntSize, UIntChain, 2048, 127> path_merge2(subproblem2.graph, subproblem2.tableau); \
                    next_problem.alignment = std::move(align(matches, subproblem1, subproblem2, \
                                                             path_merge1, path_merge2, is_main_execution))
                
                if (max_nodes < std::numeric_limits<uint32_t>::max() && max_paths < std::numeric_limits<uint8_t>::max()) {
                    _gen_packed_path_merge(uint32_t, uint8_t);
                }
                else if (max_nodes < std::numeric_limits<uint32_t>::max() && max_paths < std::numeric_limits<uint16_t>::max()) {
                    _gen_packed_path_merge(uint32_t, uint16_t);
                }
                else {
                    _gen_packed_path_merge(uint64_t, uint16_t);
                }
                #undef _gen_packed_path_merge
            }
            else {
                    
                #define _gen_path_merge(UIntSize, UIntChain) \
                    PathMerge<UIntSize, UIntChain> path_merge1(subproblem1.graph, subproblem1.tableau); \
                    PathMerge<UIntSize, UIntChain> path_merge2(subproblem2.graph, subproblem2.tableau); \
                    next_problem.alignment = std::move(align(matches, subproblem1, subproblem2, \
                                                             path_merge1, path_merge2, is_main_execution))
                
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
#endif
        }
        else {
            // use non-overlapping chains for reachability for more efficiency
            ChainMerge chain_merge1(subproblem1.graph, subproblem1.tableau);
            ChainMerge chain_merge2(subproblem2.graph, subproblem2.tableau);
            
            next_problem.alignment = std::move(align(matches, subproblem1, subproblem2,
                                                     chain_merge1, chain_merge2, is_main_execution));
        }
        
        log_memory_usage(logging::Debug);
        
        // we do this now in case we're not preserving the graphs in the subproblems
        if (!subalignments_filepath.empty() && is_main_execution) {
            emit_subalignment();
        }
        
        logging::log(logging::Verbose, "Fusing MSAs along the alignment.");
        
        // fuse either in place or in a copy
        BaseGraph fused_graph;
        if (preserve_subproblems) {
            fused_graph = subproblem1.graph;
        }
        else {
            fused_graph = std::move(subproblem1.graph);
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
        
        next_problem.graph = std::move(fused_graph);
        next_problem.tableau = subproblem1.tableau;
        
        next_problem.complete = true;
        
        if (!subproblems_prefix.empty() && is_main_execution) {
            emit_subproblem(next_problem);
        }
        
        log_memory_usage(logging::Verbose);
    }
    
    if (!is_main_execution) {
        logging::level = current_log_level;
    }
}

template<class BGraph>
std::vector<anchor_t> Core::bonds_to_chain(const BGraph& graph, const bond_interval_t& bond_interval) const {
    
    std::vector<anchor_t> chain(bond_interval.size());
    for (size_t i = 0; i < bond_interval.size(); ++i) {
        const auto& bond = bond_interval[i];
        auto& anchor = chain[i];
        
        auto path_id1 = graph.path_id(bond.path1);
        auto path_id2 = graph.path_id(bond.path2);
        
        for (size_t j = 0; j < bond.length; ++j) {
            anchor.walk1.push_back(graph.path(path_id1)[bond.offset1 + j]);
            anchor.walk2.push_back(graph.path(path_id2)[bond.offset2 + j]);
        }
        anchor.score = bond.score;
    }
    
    return chain;
}


template<class BGraph>
void Core::output_bond_alignment(const Alignment& bond_alignment, const BGraph& graph, uint64_t path_id, size_t bond_number) const {
    
    std::string bond_aln_filename = bonds_prefix + "_" + graph.path_name(path_id) + "_cigar_" + std::to_string(bond_number) + ".txt";
    std::ofstream out(bond_aln_filename);
    if (!out) {
        throw std::runtime_error("Could not write bond alignment to " + bond_aln_filename);
    }
    
    // find the boundaries
    uint64_t first1 = -1, first2 = -1, last1 = -1, last2 = -1;
    for (size_t i = 0; i < bond_alignment.size() && (first1 == -1 || first2 == -1); ++i) {
        if (first1 == -1 && bond_alignment[i].node_id1 != AlignedPair::gap) {
            first1 = bond_alignment[i].node_id1;
        }
        if (first2 == -1 && bond_alignment[i].node_id2 != AlignedPair::gap) {
            first2 = bond_alignment[i].node_id2;
        }
    }
    for (size_t i = bond_alignment.size() - 1; i < bond_alignment.size() && (last1 == -1 || last2 == -1); --i) {
        if (last1 == -1 && bond_alignment[i].node_id1 != AlignedPair::gap) {
            last1 = bond_alignment[i].node_id1;
        }
        if (last2 == -1 && bond_alignment[i].node_id2 != AlignedPair::gap) {
            last2 = bond_alignment[i].node_id2;
        }
    }
    
    if (first1 == -1) {
        // the alignment is empty
        out << '\n';
        return;
    }
    
    // add the leading indels to a full-sequence alignment
    Alignment indel_padded;
    for (size_t i = 0; i < graph.path(path_id).size(); ++i) {
        if (graph.path(path_id)[i] == first1) {
            break;
        }
        indel_padded.emplace_back(graph.path(path_id)[i], AlignedPair::gap);
    }
    for (size_t i = 0; i < graph.path(path_id).size(); ++i) {
        if (graph.path(path_id)[i] == first2) {
            break;
        }
        indel_padded.emplace_back(AlignedPair::gap, graph.path(path_id)[i]);
    }
    
    // add the actual aligned portion
    for (const auto& aln_pair : bond_alignment) {
        indel_padded.emplace_back(aln_pair);
    }
    
    // add the lagging indels (in reverse)
    size_t to_reverse = 0;
    for (size_t i = graph.path(path_id).size() - 1; i < graph.path(path_id).size(); --i) {
        if (graph.path(path_id)[i] == last2) {
            break;
        }
        indel_padded.emplace_back(AlignedPair::gap, graph.path(path_id)[i]);
        ++to_reverse;
    }
    for (size_t i = graph.path(path_id).size() - 1; i < graph.path(path_id).size(); --i) {
        if (graph.path(path_id)[i] == last1) {
            break;
        }
        indel_padded.emplace_back(graph.path(path_id)[i], AlignedPair::gap);
        ++to_reverse;
    }
    
    // put the lagging indels in forward order
    std::reverse(indel_padded.end() - to_reverse, indel_padded.end());
    
    out << explicit_cigar(indel_padded, graph, graph) << '\n';
}

}

#endif /* centrolign_core_hpp */
