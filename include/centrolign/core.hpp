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
    
    // configurable submodules
    Simplifier simplifier;
    MatchFinder match_finder;
    Anchorer anchorer;
    Stitcher stitcher;
    
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
    
    // preserve subproblems whose parent problems have been completed
    bool preserve_subproblems = true;
    
    // if non-empty prefix to give GFA output for all suproblems
    std::string subproblems_prefix;
    
    // load alignments from the prefix and start where they left off
    void restart();
    
    // the root subproblem
    const Subproblem& root_subproblem() const;
    
    const Subproblem& leaf_subproblem(const std::string& name) const;
    
    // the narrowest subproblem that includes all these sequences
    const Subproblem& subproblem_covering(const std::vector<std::string>& names) const;
    
    void calibrate_anchor_scores();
    
private:
    
    void init(std::vector<std::pair<std::string, std::string>>&& names_and_sequences,
              Tree&& tree_in);
    
    std::string subproblem_file_name(uint64_t tree_id) const;
    
    void emit_subproblem(uint64_t tree_id) const;
    
    std::vector<std::string> leaf_descendents(uint64_t tree_id) const;
    
    std::vector<match_set_t> get_matches(Subproblem& subproblem1, Subproblem& subproblem2,
                                         bool suppress_verbose_logging) const;
    
    std::vector<match_set_t> query_matches(ExpandedGraph& expanded1,
                                           ExpandedGraph& expanded2) const;
    
    std::vector<size_t> assign_reanchor_budget(const std::vector<std::pair<SubGraphInfo, SubGraphInfo>>& stitch_graphs) const;
    
    template<class XMerge>
    std::vector<anchor_t> anchor(std::vector<match_set_t>& matches, double anchor_scale,
                                 const Subproblem& subproblem1, const Subproblem& subproblem2,
                                 const XMerge& xmerge1, const XMerge& xmerge2, bool gap_calibrate) const;
    
    template<class XMerge>
    Alignment align(std::vector<match_set_t>& matches,
                    const Subproblem& subproblem1, const Subproblem& subproblem2,
                    const XMerge& xmerge1, const XMerge& xmerge2) const;
    
    template<class XMerge>
    double estimate_gap_cost_scale(std::vector<match_set_t>& matches,
                                   const Subproblem& subproblem1, const Subproblem& subproblem2,
                                   const XMerge& xmerge1, const XMerge& xmerge2) const;
    
    Tree tree;
    
    std::vector<Subproblem> subproblems;
    
};





/*
 * Template and inline implementations
 */

template<class XMerge>
std::vector<anchor_t> Core::anchor(std::vector<match_set_t>& matches, double anchor_scale,
                                   const Subproblem& subproblem1, const Subproblem& subproblem2,
                                   const XMerge& xmerge1, const XMerge& xmerge2, bool gap_calibrate) const {
    
    logging::log(gap_calibrate ? logging::Debug : logging::Verbose, "Anchoring.");
    
    // anchor the alignment
    std::vector<anchor_t> anchors;
    if (gap_calibrate) {
        // override with algorithm with unscored gaps
        anchors = std::move(anchorer.anchor_chain(matches, subproblem1.graph, subproblem2.graph,
                                                  xmerge1, xmerge2, Anchorer::Sparse));
    }
    else {
        anchors = std::move(anchorer.anchor_chain(matches,
                                                  subproblem1.graph, subproblem2.graph,
                                                  xmerge1, xmerge2, anchor_scale));
    }
    
    {
        logging::log(gap_calibrate ? logging::Debug : logging::Verbose, "Extracting subgraphs for fill-in anchoring");
        
        auto stitch_graphs = stitcher.extract_stitch_graphs(anchors, subproblem1.graph, subproblem2.graph,
                                                            subproblem1.tableau, subproblem2.tableau,
                                                            xmerge1, xmerge2);
        
        stitcher.project_paths(subproblem1.graph, subproblem2.graph, stitch_graphs);
        
        auto stitch_matches = stitcher.divvy_matches(matches, subproblem1.graph, subproblem2.graph,
                                                     stitch_graphs);
        
        auto budgets = assign_reanchor_budget(stitch_graphs);
        
        logging::log(gap_calibrate ? logging::Debug : logging::Verbose, "Performing anchoring subproblems for fill-in anchoring");
        
        std::vector<std::vector<anchor_t>> stitch_anchors(stitch_graphs.size());
        for (size_t i = 0; i < stitch_graphs.size(); ++i) {
            XMerge stitch_xmerge1(stitch_graphs[i].first.subgraph);
            XMerge stitch_xmerge2(stitch_graphs[i].second.subgraph);
            if (gap_calibrate) {
                // override with algorithm with unscored gaps
                stitch_anchors[i] = std::move(anchorer.global_anchor_chain(stitch_matches[i],
                                                                           stitch_graphs[i].first.subgraph,
                                                                           stitch_graphs[i].second.subgraph,
                                                                           stitch_xmerge1, stitch_xmerge2,
                                                                           stitch_graphs[i].first.sources,
                                                                           stitch_graphs[i].second.sources,
                                                                           stitch_graphs[i].first.sinks,
                                                                           stitch_graphs[i].second.sinks,
                                                                           Anchorer::Sparse,
                                                                           budgets[i]));
            }
            else {
                stitch_anchors[i] = std::move(anchorer.global_anchor_chain(stitch_matches[i],
                                                                           stitch_graphs[i].first.subgraph,
                                                                           stitch_graphs[i].second.subgraph,
                                                                           stitch_xmerge1, stitch_xmerge2,
                                                                           stitch_graphs[i].first.sources,
                                                                           stitch_graphs[i].second.sources,
                                                                           stitch_graphs[i].first.sinks,
                                                                           stitch_graphs[i].second.sinks,
                                                                           anchor_scale,
                                                                           budgets[i]));
            }
        }
        
        stitcher.merge_stitch_chains(anchors, stitch_anchors, stitch_graphs);
        
        logging::log(logging::Debug, "Filled-in anchor chain consists of " + std::to_string(anchors.size()) + " anchors");
    }
    
    return anchors;
}

template<class XMerge>
Alignment Core::align(std::vector<match_set_t>& matches,
                      const Subproblem& subproblem1, const Subproblem& subproblem2,
                      const XMerge& xmerge1, const XMerge& xmerge2) const {
    
    double scale = anchorer.global_scale;
    if (anchorer.chaining_algorithm == Anchorer::SparseAffine) {
        logging::log(logging::Verbose, "Pre-anchoring to calibrate parameter scale");
        
        scale = estimate_gap_cost_scale(matches, subproblem1, subproblem2, xmerge1, xmerge2);
        
        logging::log(logging::Debug, "Estimate anchoring scale of " + std::to_string(scale));
    }
    
    auto anchors = anchor(matches, scale, subproblem1, subproblem2, xmerge1, xmerge2, false);
    
    static const bool instrument = true;
    if (instrument) {
        anchorer.instrument_anchor_chain(anchors, scale, subproblem1.graph, subproblem2.graph, xmerge1, xmerge2);
    }
    
    logging::log(logging::Verbose, "Stitching anchors into alignment");
    
    // form a base-level alignment
    return stitcher.stitch(anchors, subproblem1.graph, subproblem2.graph,
                           subproblem1.tableau, subproblem2.tableau,
                           xmerge1, xmerge2);
}

template<class XMerge>
double Core::estimate_gap_cost_scale(std::vector<match_set_t>& matches,
                                     const Subproblem& subproblem1, const Subproblem& subproblem2,
                                     const XMerge& xmerge1, const XMerge& xmerge2) const {
    
    // get an anchoring with unscored gaps
    auto anchors = anchor(matches, anchorer.global_scale, subproblem1, subproblem2, xmerge1, xmerge2, true);
    
    // measure its weight
    double total_weight = 0.0;
    size_t total_length = 0;
    for (const auto& anchor : anchors) {
        total_weight += anchorer.anchor_weight(anchor);
        total_length += anchor.walk1.size();
    }
    
    // estimate the length of the in-between bits
    auto stitch_graphs = stitcher.extract_stitch_graphs(anchors, subproblem1.graph, subproblem2.graph,
                                                        subproblem1.tableau, subproblem2.tableau,
                                                        xmerge1, xmerge2);
    for (const auto& stitch_pair : stitch_graphs) {
        
        size_t stitch_length = std::numeric_limits<size_t>::max();
        for (auto subgraph_ptr : {&stitch_pair.first, &stitch_pair.second}) {
            // compute the minimum distance across either of the stitch graphs
            const auto& subgraph = *subgraph_ptr;
            if (subgraph.subgraph.node_size() == 0) {
                stitch_length = 0;
            }
            else {
                auto minmax_dists = minmax_distance(subgraph.subgraph, &subgraph.sources);
                for (auto sink : subgraph.sinks) {
                    stitch_length = std::min<size_t>(stitch_length, minmax_dists[sink].first);
                }
            }
        }
        
        total_length += stitch_length;
    }
    
    return total_weight / total_length;
}


}

#endif /* centrolign_core_hpp */
