#ifndef centrolign_core_hpp
#define centrolign_core_hpp

#include <vector>
#include <string>

#include "centrolign/modify_graph.hpp"
#include "centrolign/anchorer.hpp"
#include "centrolign/stitcher.hpp"
#include "centrolign/tree.hpp"
#include "centrolign/simplifier.hpp"
#include "centrolign/match_finder.hpp"
#include "centrolign/partitioner.hpp"
#include "centrolign/score_function.hpp"
#include "centrolign/utility.hpp"
#include "centrolign/bonder.hpp"

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
    // simplifies graph topology in advance of querying matches
    Simplifier simplifier;
    // queries matches between the input graphs
    MatchFinder match_finder;
    // makes a chain of alignment anchors using matches
    Anchorer anchorer;
    // partitions anchor chains into well-anchored and poorly-anchored portions
    Partitioner partitioner;
    // aligns in between anchors to them stitch into a base-level alignment
    Stitcher stitcher;
    // identifies sequences to bond together as cycles (if cyclizing)
    Bonder bonder;
    
    
    /*
     * An alignment subproblem in the process of making the full MSA
     */
    struct Subproblem {
        Subproblem() noexcept = default;
        Subproblem(const Subproblem& other) noexcept = default;
        Subproblem(Subproblem&& other) noexcept = default;
        ~Subproblem() = default;
        Subproblem& operator=(const Subproblem& other) noexcept = default;
        Subproblem& operator=(Subproblem&& other) noexcept = default;
        
        BaseGraph graph;
        SentinelTableau tableau;
        Alignment alignment;
        std::string name;
        bool complete = false;
    };
    
    // don't calibrate the scale of the scoring function before executing
    bool skip_calibration = false;
    
    // preserve subproblems whose parent problems have been completed
    bool preserve_subproblems = true;
    
    // if non-empty, prefix to give GFA output for all suproblems
    std::string subproblems_prefix;
    
    // if non-empty, file to write suproblem alignments to
    std::string subalignments_filepath;
    
    // load alignments from the prefix and start where they left off
    void restart();
    
    // the root subproblem
    const Subproblem& root_subproblem() const;
    
    // the leaf subproblem that corresponds to a sequences
    const Subproblem& leaf_subproblem(const std::string& name) const;
    
    // the narrowest subproblem that includes all these sequences
    const Subproblem& subproblem_covering(const std::vector<std::string>& names) const;
    
    // learn the intrinsic scale of the anchor scoring function on these sequences
    void calibrate_anchor_scores();
    
private:
    
    void init(std::vector<std::pair<std::string, std::string>>&& names_and_sequences,
              Tree&& tree_in);
    
    std::string subproblem_file_name(uint64_t tree_id) const;
    
    std::string subproblem_info_file_name() const;
    
    void emit_subproblem(uint64_t tree_id) const;
    
    void emit_subalignment(uint64_t tree_id) const;
    
    std::vector<std::string> leaf_descendents(uint64_t tree_id) const;
    
    std::vector<match_set_t> get_matches(Subproblem& subproblem1, Subproblem& subproblem2,
                                         bool suppress_verbose_logging) const;
    
    std::vector<match_set_t> query_matches(ExpandedGraph& expanded1,
                                           ExpandedGraph& expanded2) const;
    
    template<class XMerge>
    Alignment align(std::vector<match_set_t>& matches,
                    const Subproblem& subproblem1, const Subproblem& subproblem2,
                    XMerge& xmerge1, XMerge& xmerge2) const;
    
    void log_memory_usage(logging::LoggingLevel level) const;
    
    // the guide tree
    Tree tree;
    
    // the individual alignment subproblems (including single-sequence leaves)
    std::vector<Subproblem> subproblems;
    
};





/*
 * Template and inline implementations
 */

template<class XMerge>
Alignment Core::align(std::vector<match_set_t>& matches,
                      const Subproblem& subproblem1, const Subproblem& subproblem2,
                      XMerge& xmerge1, XMerge& xmerge2) const {
    
    if (logging::level >= logging::Debug) {
        size_t merge_size = xmerge1.memory_size() + xmerge2.memory_size();
        logging::log(logging::Debug, "Merge data structures are consuming " + format_memory_usage(merge_size) + " of memory.");
        log_memory_usage(logging::Debug);
    }
    
    // get the best anchor chain
    auto anchors = anchorer.anchor_chain(matches, subproblem1.graph, subproblem2.graph,
                                         subproblem1.tableau, subproblem2.tableau,
                                         xmerge1, xmerge2);
    
    log_memory_usage(logging::Debug);
    
    static const bool output_anchors = false;
    if (output_anchors) {
        std::cerr << "outputting anchors\n";
        for (const auto& a : anchors) {
            std::cout << a.walk1.front() << '\t' << a.walk2.front() << '\t' << a.walk1.size() << '\t' << 0 << '\n';
        }
        
        std::cerr << "making initial mask\n";
        
        std::unordered_set<std::tuple<size_t, size_t, size_t>> mask;
        anchorer.update_mask(matches, anchors, mask);
        
        size_t num_reanchors = 1;
        bool mask_reciprocal = true;
        for (size_t i = 1; i <= num_reanchors; ++i) {
            
            auto anchors_secondary = anchorer.anchor_chain(matches, subproblem1.graph, subproblem2.graph,
                                                           subproblem1.tableau, subproblem2.tableau,
                                                           xmerge1, xmerge2, &mask);
            anchorer.update_mask(matches, anchors_secondary, mask, mask_reciprocal);
            
            std::unordered_map<uint64_t, std::pair<uint64_t, size_t>> paired_node_ids;
            for (size_t j = 0; j < anchors.size(); ++j) {
                const auto& anchor = anchors[j];
                for (size_t i = 0; i < anchor.walk1.size(); ++i) {
                    paired_node_ids[anchor.walk1[i]] = std::make_pair(anchor.walk2[i], j);
                }
            }
            for (size_t j = 0; j < anchors_secondary.size(); ++j) {
                const auto& anchor = anchors_secondary[j];
                for (size_t i = 0; i < anchor.walk1.size(); ++i) {
                    if (paired_node_ids.count(anchor.walk1[i])) {
                        if (paired_node_ids[anchor.walk1[i]].first == anchor.walk2[i]) {
                            std::cerr << "!! WARNING: prohibited match on masked anchor at nodes " << anchor.walk1[i] << ", " << anchor.walk2[i] << " between opt anchor " << paired_node_ids[anchor.walk1[i]].second << " and secondary anchor " << j << '\n';
                            exit(1);
                        }
                    }
                }
            }
            
            for (const auto& a : anchors_secondary) {
                std::cout << a.walk1.front() << '\t' << a.walk2.front() << '\t' << a.walk1.size() << '\t' << i << '\n';
            }
            
            auto bonds = bonder.identify_bonds(subproblem1.graph, subproblem2.graph,
                                               subproblem1.tableau, subproblem2.tableau,
                                               xmerge1, xmerge2, anchors, anchors_secondary);
            
//            std::cerr << "instrumenting partition for chain " << i << "\n";
//            partitioner.partition_anchors(anchors_secondary, subproblem1.graph, subproblem2.graph,
//                                          subproblem1.tableau, subproblem2.tableau,
//                                          xmerge1, xmerge2);
        }
        
        exit(0);
    }
        
    // partition the anchor chain into good and bad segments
    auto anchor_segments = partitioner.partition_anchors(anchors, subproblem1.graph, subproblem2.graph,
                                                         subproblem1.tableau, subproblem2.tableau,
                                                         xmerge1, xmerge2);
    
    log_memory_usage(logging::Debug);
    
    logging::log(logging::Verbose, "Stitching anchors into alignment.");
    
    // form a base-level alignment
    Alignment alignment = stitcher.stitch(anchor_segments, subproblem1.graph, subproblem2.graph,
                                          subproblem1.tableau, subproblem2.tableau,
                                          xmerge1, xmerge2);
    
    return alignment;
}


}

#endif /* centrolign_core_hpp */
