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
#include "centrolign/minmax_distance.hpp"
#include "centrolign/partitioner.hpp"
#include "centrolign/score_function.hpp"

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
    
    // simplifies graph topology in advance of querying matches
    Simplifier simplifier;
    // queries matches between the input graphs
    MatchFinder match_finder;
    // a function for scoring anchors
    ScoreFunction score_function;
    // makes a chain of alignment anchors using matches
    Anchorer anchorer;
    // partitions anchor chains into well-anchored and poorly-anchored portions
    Partitioner partitioner;
    // aligns in between anchors to them stitch into a base-level alignment
    Stitcher stitcher;
    
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
                    const XMerge& xmerge1, const XMerge& xmerge2) const;
    
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
                      const XMerge& xmerge1, const XMerge& xmerge2) const {
    
    // get the best anchor chain
    auto anchors = anchorer.anchor_chain(matches, subproblem1.graph, subproblem2.graph,
                                         subproblem1.tableau, subproblem2.tableau,
                                         xmerge1, xmerge2);
    
    static const bool output_anchors = false;
    if (output_anchors) {
        std::unordered_set<std::tuple<size_t, size_t, size_t>> mask;
        mask.reserve(anchors.size());
        for (const auto& a : anchors) {
            std::cout << a.walk1.front() << '\t' << a.walk2.front() << '\t' << a.walk1.size() << '\t' << 0 << '\n';
            mask.emplace(a.match_set, a.idx1, a.idx2);
        }
    }
    
    // partition the anchor chain into good and bad segments
    auto anchor_segments = partitioner.partition_anchors(anchors, subproblem1.graph, subproblem2.graph,
                                                         subproblem1.tableau, subproblem2.tableau,
                                                         xmerge1, xmerge2);

    if (output_anchors) {
        // wait to do this until after partitioning so we can instrument the partitioner
        exit(0);
    }
    
    logging::log(logging::Verbose, "Stitching anchors into alignment.");
    
    // form a base-level alignment
    return stitcher.stitch(anchor_segments, subproblem1.graph, subproblem2.graph,
                           subproblem1.tableau, subproblem2.tableau,
                           xmerge1, xmerge2);
}


}

#endif /* centrolign_core_hpp */
