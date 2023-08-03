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
        Subproblem() = default;
        ~Subproblem() = default;
        
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
    
private:
    
    void init(std::vector<std::pair<std::string, std::string>>&& names_and_sequences,
              Tree&& tree_in);
    
    std::string subproblem_file_name(uint64_t tree_id) const;
    
    void emit_subproblem(uint64_t tree_id) const;
    
    std::vector<std::string> leaf_descendents(uint64_t tree_id) const;
    
    std::vector<match_set_t> find_matches(ExpandedGraph& expanded1,
                                          ExpandedGraph& expanded2) const;
    
    std::vector<size_t> assign_reanchor_budget(const std::vector<std::pair<SubGraphInfo, SubGraphInfo>>& stitch_graphs) const;
    
    template<class XMerge>
    Alignment align(std::vector<match_set_t>& matches,
                    const Subproblem& subproblem1, const Subproblem& subproblem2,
                    const XMerge& xmerge1, const XMerge& xmerge2) const;
    
    
    
    Tree tree;
    
    std::vector<Subproblem> subproblems;
    
};





/*
 * Template and inline implementations
 */

template<class XMerge>
Alignment Core::align(std::vector<match_set_t>& matches,
                      const Subproblem& subproblem1, const Subproblem& subproblem2,
                      const XMerge& xmerge1, const XMerge& xmerge2) const {
    
    
    logging::log(logging::Verbose, "Anchoring");
    
    // anchor the alignment
    auto anchors = anchorer.anchor_chain(matches,
                                         subproblem1.graph, subproblem2.graph,
                                         xmerge1, xmerge2);
    
    {
        logging::log(logging::Verbose, "Extracting subgraphs for fill-in anchoring");
        
        auto stitch_graphs = stitcher.extract_stitch_graphs(anchors, subproblem1.graph, subproblem2.graph,
                                                            subproblem1.tableau, subproblem2.tableau,
                                                            xmerge1, xmerge2);
        
        stitcher.project_paths(subproblem1.graph, subproblem2.graph, stitch_graphs);
        
        auto stitch_matches = stitcher.divvy_matches(matches, subproblem1.graph, subproblem2.graph,
                                                     stitch_graphs);
        
        auto budgets = assign_reanchor_budget(stitch_graphs);
        
        logging::log(logging::Verbose, "Performing anchoring subproblems for fill-in anchoring");
        
        // TODO: make an optional override for max_num_match_pairs in Anchorer
        std::vector<std::vector<anchor_t>> stitch_anchors(stitch_graphs.size());
        for (size_t i = 0; i < stitch_graphs.size(); ++i) {
            XMerge stitch_xmerge1(stitch_graphs[i].first.subgraph);
            XMerge stitch_xmerge2(stitch_graphs[i].second.subgraph);
            stitch_anchors[i] = std::move(anchorer.anchor_chain(stitch_matches[i],
                                                                stitch_graphs[i].first.subgraph,
                                                                stitch_graphs[i].second.subgraph,
                                                                stitch_xmerge1, stitch_xmerge2,
                                                                budgets[i], true));
        }
        
        stitcher.merge_stitch_chains(anchors, stitch_anchors, stitch_graphs);
        
        logging::log(logging::Debug, "Filled-in anchor chain consists of " + std::to_string(anchors.size()) + " anchors");
    }
    
    static const bool instrument = true;
    if (instrument) {
        anchorer.instrument_anchor_chain(anchors, subproblem1.graph, subproblem2.graph, xmerge1, xmerge2);
    }
    
    logging::log(logging::Verbose, "Stitching anchors into alignment");
    
    // form a base-level alignment
    return stitcher.stitch(anchors, subproblem1.graph, subproblem2.graph,
                           subproblem1.tableau, subproblem2.tableau,
                           xmerge1, xmerge2);
}


}

#endif /* centrolign_core_hpp */
