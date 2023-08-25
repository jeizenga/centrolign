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
#include "centrolign/partition.hpp"

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
    
    auto anchors = anchorer.anchor_chain(matches, subproblem1.graph, subproblem2.graph,
                                         subproblem1.tableau, subproblem2.tableau,
                                         xmerge1, xmerge2);
    
//    // TODO: normalize the threshold for the scoring function
//    // - both postive area under curve (over lengths) and peak value produce similar scales (up to a few % different)
//    static const bool instrument_partition = true;
//    if (instrument_partition) {
//
//
//        auto stitch_graphs = stitcher.extract_stitch_graphs(anchors, subproblem1.graph, subproblem2.graph,
//                                                            subproblem1.tableau, subproblem2.tableau,
//                                                            xmerge1, xmerge2);
//
//
//        std::vector<std::pair<double, double>> partition_data(anchors.size() + stitch_graphs.size());
//
//        for (size_t i = 0; i < partition_data.size(); ++i) {
//            if (i % 2 == 0) {
//                const auto& stitch_pair = stitch_graphs[i / 2];
//                size_t stitch_length = std::numeric_limits<size_t>::max();
//                for (auto subgraph_ptr : {&stitch_pair.first, &stitch_pair.second}) {
//                    // compute the minimum distance across either of the stitch graphs
//                    const auto& subgraph = *subgraph_ptr;
//                    if (subgraph.subgraph.node_size() == 0) {
//                        stitch_length = 0;
//                    }
//                    else {
//                        auto minmax_dists = minmax_distance(subgraph.subgraph, &subgraph.sources);
//                        for (auto sink : subgraph.sinks) {
//                            stitch_length = std::min<size_t>(stitch_length, minmax_dists[sink].first);
//                        }
//                    }
//                }
//
//                partition_data[i].first = 0.0;
//                partition_data[i].second = std::max<double>(1.0, stitch_length);
//            }
//            else {
//                partition_data[i].first = anchorer.anchor_weight(anchors[i / 2]);
//                partition_data[i].second = anchors[i / 2].walk1.size();
//            }
//        }
//
//        auto partition = average_constrained_partition(partition_data, 0.15 / anchorer.global_scale, 4000.0 / anchorer.global_scale);
//
//        for (size_t i = 0; i < partition.size(); ++i) {
//            auto p = partition[i];
//            double weight = 0.0;
//            for (size_t j = p.first / 2; j < p.second / 2; ++j) {
//                weight += anchorer.anchor_weight(anchors[j]);
//            }
//            std::cerr << '|' << '\t' << (p.first / 2) << '\t' << (p.second / 2) << '\t' << weight << '\t' << anchors[p.first / 2].walk1.front() << '\t' << anchors[p.second / 2 - 1].walk1.back() << '\t' << anchors[p.first / 2].walk2.front() << '\t' << anchors[p.second / 2 - 1].walk2.back() << '\t' << (anchors[p.second / 2 - 1].walk1.back() - anchors[p.first / 2].walk1.front()) << '\t' << (anchors[p.second / 2 - 1].walk2.back() - anchors[p.first / 2].walk2.front());
//            if (i != 0) {
//                auto q = partition[i - 1];
//                std::cerr << '\t' << (anchors[p.first / 2].walk1.front() - anchors[q.second / 2 - 1].walk1.back()) << '\t' << (anchors[p.first / 2].walk2.front() - anchors[q.second / 2 - 1].walk2.back());
//            }
//            else {
//                std::cerr << '\t' << 0 << '\t' << 0;
//            }
//            std::cerr << '\n';
//        }
//    }
    
    
//    for (const auto& a : anchors) {
//        std::cout << a.walk1.front() << '\t' << a.walk2.front() << '\t' << a.walk1.size() << '\n';
//    }
//    exit(0);
    
    logging::log(logging::Verbose, "Stitching anchors into alignment");
    
    // form a base-level alignment
    return stitcher.stitch(anchors, subproblem1.graph, subproblem2.graph,
                           subproblem1.tableau, subproblem2.tableau,
                           xmerge1, xmerge2);
}


}

#endif /* centrolign_core_hpp */
