#ifndef centrolign_anchorer_hpp
#define centrolign_anchorer_hpp

#include <vector>
#include <cstdint>
#include <cmath>
#include <iostream>
#include <limits>
#include <array>
#include <memory>

#include "centrolign/chain_merge.hpp"
#include "centrolign/modify_graph.hpp"
#include "centrolign/gesa.hpp"
#include "centrolign/max_search_tree.hpp"
#include "centrolign/orthogonal_max_search_tree.hpp"
#include "centrolign/topological_order.hpp"
#include "centrolign/match_finder.hpp"
#include "centrolign/subgraph_extraction.hpp"
#include "centrolign/step_index.hpp"
#include "centrolign/score_function.hpp"
#include "centrolign/utility.hpp"
#include "centrolign/alignment.hpp"
#include "centrolign/superbubbles.hpp"
#include "centrolign/structure_distances.hpp"
#include "centrolign/match_bank.hpp"
#include "centrolign/packed_match_bank.hpp"

namespace centrolign {

// a pair of walks of the same sequence in two graphs
struct anchor_t {
    anchor_t() noexcept = default;
    anchor_t(const anchor_t& other) noexcept = default;
    anchor_t(anchor_t&& other) noexcept = default;
    ~anchor_t() = default;
    anchor_t& operator=(const anchor_t& other) noexcept = default;
    anchor_t& operator=(anchor_t&& other) noexcept = default;
    
    std::vector<uint64_t> walk1;
    std::vector<uint64_t> walk2;
    size_t count1 = 0;
    size_t count2 = 0;
    size_t full_length = 0;
    double score = 0.0;
    int64_t gap_before = 0;
    int64_t gap_after = 0;
    double gap_score_before = 0.0;
    double gap_score_after = 0.0;
    size_t match_set = -1;
    size_t idx1 = -1;
    size_t idx2 = -1;
};

/*
 * An interface used by objects that handle subgraphs between anchors
 */
class Extractor {
protected:
    Extractor() = default;
    ~Extractor() = default;
    
    // pull out the graphs between anchors
    template<class BGraph1, class BGraph2, class XMerge1, class XMerge2>
    static std::vector<std::pair<SubGraphInfo, SubGraphInfo>>
    extract_graphs_between(const std::vector<anchor_t>& anchor_chain,
                           const BGraph1& graph1, const BGraph2& graph2,
                           const SentinelTableau& tableau1, const SentinelTableau& tableau2,
                           const XMerge1& chain_merge1, const XMerge2& chain_merge2);
    
    // pull out the graphs between anchors, but not the ones before and after the chain
    template<class BGraph1, class BGraph2, class XMerge1, class XMerge2>
    static std::vector<std::pair<SubGraphInfo, SubGraphInfo>>
    extract_graphs_between(const std::vector<anchor_t>& anchor_chain,
                           const BGraph1& graph1, const BGraph2& graph2,
                           const XMerge1& chain_merge1, const XMerge2& chain_merge2);
    
    // pull out the graph between anchor segments, consists of a vector of between-graphs for
    // each segment and a vector of between-graphs that is between each segment pair
    template<class BGraph1, class BGraph2, class XMerge1, class XMerge2>
    static std::pair<std::vector<std::vector<std::pair<SubGraphInfo, SubGraphInfo>>>,
                     std::vector<std::pair<SubGraphInfo, SubGraphInfo>>>
    extract_graphs_between(const std::vector<std::vector<anchor_t>>& anchor_segments,
                           const BGraph1& graph1, const BGraph2& graph2,
                           const SentinelTableau& tableau1, const SentinelTableau& tableau2,
                           const XMerge1& chain_merge1, const XMerge2& chain_merge2);
    
    // add paths onto stitch graphs
    template<class BGraph1, class BGraph2>
    static void project_paths(const BGraph1& graph1, const BGraph2& graph2,
                              std::vector<std::pair<SubGraphInfo, SubGraphInfo>>& stitch_graphs);
    
    // (min, max) distance from any source to any sink
    static std::pair<int64_t, int64_t> source_sink_minmax(const SubGraphInfo& extraction);
    
    static std::vector<size_t> get_logging_indexes(size_t size);
    
private:
    
    
    template<class BGraph1, class BGraph2, class XMerge1, class XMerge2>
    static std::vector<std::pair<SubGraphInfo, SubGraphInfo>>
    extract_graphs_between_internal(const std::vector<anchor_t>& anchor_chain,
                                    const BGraph1& graph1, const BGraph2& graph2,
                                    const SentinelTableau* tableau1, const SentinelTableau* tableau2,
                                    const XMerge1& chain_merge1, const XMerge2& chain_merge2);
    
    template<class BGraph1, class BGraph2, class XMerge1, class XMerge2>
    static std::pair<SubGraphInfo, SubGraphInfo>
    do_extraction(uint64_t from1, uint64_t to1, uint64_t from2, uint64_t to2,
                  const BGraph1& graph1, const BGraph2& graph2,
                  const XMerge1& chain_merge1, const XMerge2& chain_merge2);
    
    template<class BGraph>
    static void do_project(const BGraph& graph, SubGraphInfo& subgraph,
                           const StepIndex& step_index);
};


/*
 * Data structure finding anchors between two graphs
 */
class Anchorer : public Extractor {
public:
    Anchorer(const ScoreFunction& score_function) : score_function(&score_function) {}
    Anchorer() = default;
    ~Anchorer() = default;
    
    enum ChainAlgorithm { Exhaustive = 0, Sparse = 1, SparseAffine = 2 };
    
    template<class BGraph, class XMerge>
    std::vector<anchor_t> anchor_chain(std::vector<match_set_t>& matches,
                                       const BGraph& graph1,
                                       const BGraph& graph2,
                                       const SentinelTableau& tableau1,
                                       const SentinelTableau& tableau2,
                                       const XMerge& xmerge1,
                                       const XMerge& xmerge2,
                                       bool restrain_memory,
                                       std::unordered_set<std::tuple<size_t, size_t, size_t>>* masked_matches = nullptr,
                                       double* override_scale = nullptr) const;
    
    /*
     * Configurable parameters
     */
    
    // select which algorithm to use
    ChainAlgorithm chaining_algorithm = SparseAffine;
    
    // follow anchoring by reanchoring between anchors using previously excluded matches
    bool do_fill_in_anchoring = true;
    // try to adjust gap penalties to the scale of the anchor scores
    bool autocalibrate_gap_penalties = true;
    // force anchors to start at the beginning and end sentinels
    bool global_anchoring = true;
    // affine gap parameters
    std::array<double, 3> gap_open{1.25, 50.0, 5000.0};
    std::array<double, 3> gap_extend{2.5, 0.1, 0.0015};
    // the max number of match pairs we will use for anchoring
    size_t max_num_match_pairs = 1000000;
    
    // split anchors at branch positions in the graph to avoid reachability artifacts
    bool split_matches_at_branchpoints = true;
    // only split anchors within this distance of the ends of the anchor
    size_t anchor_split_limit = 5;
    // only split anchors that are at least this long
    size_t min_split_length = 128;
    // only split anchors that border bubbles with path lengths that differ by at least this much
    size_t min_path_length_spread = 50;
    // only split anchors from match sets of at most this size
    size_t max_split_match_set_size = 16;
    
protected:
    
    const ScoreFunction* const score_function = nullptr;
    
    template<class BGraph, class XMerge>
    std::vector<anchor_t> anchor_chain(std::vector<match_set_t>& matches,
                                       const BGraph& graph1,
                                       const BGraph& graph2,
                                       const SentinelTableau& tableau1,
                                       const SentinelTableau& tableau2,
                                       const XMerge& xmerge1,
                                       const XMerge& xmerge2,
                                       bool restrain_memory,
                                       ChainAlgorithm local_chaining_algorithm,
                                       bool suppress_verbose_logging,
                                       double anchor_scale,
                                       std::unordered_set<std::tuple<size_t, size_t, size_t>>* masked_matches) const;
    
    
    // compute a heaviest weight anchoring of a set of matches
    // may reorder matches
    // if sources and sinks are provided, performs global anchoring as opposed to
    // local anchoring
    template<class BGraph, class XMerge>
    std::vector<anchor_t> anchor_chain(std::vector<match_set_t>& matches,
                                       const BGraph& graph1,
                                       const BGraph& graph2,
                                       const XMerge& chain_merge1,
                                       const XMerge& chain_merge2,
                                       bool restrain_memory,
                                       const std::vector<uint64_t>* sources1,
                                       const std::vector<uint64_t>* sources2,
                                       const std::vector<uint64_t>* sinks1,
                                       const std::vector<uint64_t>* sinks2,
                                       size_t local_max_num_match_pairs,
                                       bool suppress_verbose_logging,
                                       ChainAlgorithm local_chaining_algorithm,
                                       double anchor_scale,
                                       std::unordered_set<std::tuple<size_t, size_t, size_t>>* masked_matches) const;

    static const bool debug_anchorer;
    
    
    // a pair of two of the occurrence of a match
    struct AnchorNode {
        AnchorNode(size_t set, size_t idx1, size_t idx2, double weight, double initial_weight = 0.0, double final_weight = 0.0) :
            set(set), idx1(idx1), idx2(idx2), weight(weight), in_degree(0), initial_weight(initial_weight), final_weight(final_weight) { }
        AnchorNode() = default;
        ~AnchorNode() = default;
        size_t set = 0;
        size_t idx1 = 0;
        size_t idx2 = 0;
        double weight = 0.0;
        size_t in_degree = 0; // for the sake of satisfying the topological_order interface
        std::vector<size_t> edges;
        std::vector<double> edge_weights;
        double initial_weight = 0.0;
        double final_weight = 0.0;
    };
    
    // mutual reachability graph over anchor pairs
    class AnchorGraph {
    public:
        AnchorGraph() = default;
        ~AnchorGraph() = default;
        
        uint64_t add_node(size_t set, size_t idx1, size_t idx2, double weight, double initial_weight = 0.0, double final_weight = 0.0);
        void add_edge(uint64_t from_id, uint64_t to_id, double weight = 0.0);
        
        std::vector<uint64_t> heaviest_weight_path(double min_score = 0.0) const;
        // get (set, idx1, idx2)
        std::tuple<size_t, size_t, size_t> label(uint64_t node_id) const;
        
        size_t node_size() const;
        const std::vector<size_t>& next(uint64_t node_id) const;
        size_t next_size(uint64_t node_id) const;
        size_t previous_size(uint64_t node_id) const;
        
    private:
        
        std::vector<AnchorNode> nodes;
    };
    
    // a dynamic programming value and backpointer of (anchor set, walk1 index, walk2 index)
    template<typename UIntSet, typename UIntMatch>
    using dp_entry_t = std::tuple<double, UIntSet, UIntMatch, UIntMatch>;
    
    // assumes that the graphs have already been given unique sentinels
    // note: these will never show up in anchors because they can't match
    
    template<typename UIntSet, typename UIntMatch, class BGraph, class XMerge>
    std::vector<anchor_t> exhaustive_chain_dp(const std::vector<match_set_t>& match_sets,
                                              const BGraph& graph1,
                                              const BGraph& graph2,
                                              const XMerge& chain_merge1,
                                              const XMerge& chain_merge2,
                                              bool score_edges,
                                              double local_scale,
                                              size_t num_match_sets,
                                              const std::vector<uint64_t>* sources1 = nullptr,
                                              const std::vector<uint64_t>* sources2 = nullptr,
                                              const std::vector<uint64_t>* sinks1 = nullptr,
                                              const std::vector<uint64_t>* sinks2 = nullptr,
                                              const std::unordered_set<std::tuple<size_t, size_t, size_t>>* masked_matches = nullptr) const;
    
    template<typename UIntDist, typename UIntSet, typename UIntMatch, typename UIntAnchor, typename ScoreFloat,
             class DistMatchVector, class AnchorVector, class MBank, class BGraph, class XMerge>
    std::vector<anchor_t> sparse_chain_dp(const std::vector<match_set_t>& match_sets,
                                          const BGraph& graph1,
                                          const XMerge& chain_merge1,
                                          const XMerge& chain_merge2,
                                          size_t num_match_sets,
                                          bool suppress_verbose_logging,
                                          const std::vector<uint64_t>* sources1 = nullptr,
                                          const std::vector<uint64_t>* sources2 = nullptr,
                                          const std::vector<uint64_t>* sinks1 = nullptr,
                                          const std::vector<uint64_t>* sinks2 = nullptr,
                                          const std::unordered_set<std::tuple<size_t, size_t, size_t>>* masked_matches = nullptr) const;
    
    template<typename UIntSet, typename UIntMatch, typename UIntDist, typename IntShift, typename UIntAnchor, typename ScoreFloat,
             class ShiftMatchVector, class DistMatchVector, class DistVector, class AnchorVector, class MBank, class BGraph, class XMerge, size_t NumPW>
    std::vector<anchor_t> sparse_affine_chain_dp(const std::vector<match_set_t>& match_sets,
                                                 const BGraph& graph1,
                                                 const BGraph& graph2,
                                                 const XMerge& xmerge1,
                                                 const XMerge& xmerge2,
                                                 const std::array<double, NumPW>& gap_open,
                                                 const std::array<double, NumPW>& gap_extend,
                                                 double local_scale,
                                                 size_t num_match_sets,
                                                 bool suppress_verbose_logging,
                                                 const std::vector<uint64_t>* sources1 = nullptr,
                                                 const std::vector<uint64_t>* sources2 = nullptr,
                                                 const std::vector<uint64_t>* sinks1 = nullptr,
                                                 const std::vector<uint64_t>* sinks2 = nullptr,
                                                 const std::unordered_set<std::tuple<size_t, size_t, size_t>>* masked_matches = nullptr) const;
    
    // the shortest distance along any chain to each node after jumping from a given chain
    template<typename UIntDist, class BGraph, class XMerge>
    std::vector<std::vector<UIntDist>> post_switch_distances(const BGraph& graph, const XMerge& xmerge) const;
    
    template<class BGraph>
    std::pair<std::vector<bool>, std::vector<bool>> generate_forward_edge_masks(const BGraph& graph1, const std::vector<match_set_t>& match_sets, size_t num_match_sets) const;
    
    template<typename ScoreFloat, class FinalFunc, class MBank>
    std::vector<anchor_t> traceback_sparse_dp(const std::vector<match_set_t>& match_sets,
                                              const MBank& match_bank,
                                              const FinalFunc& final_function, ScoreFloat min_score,
                                              bool suppress_verbose_logging) const;
    
    
    template<class BGraph1, class BGraph2>
    void split_branching_matches(std::vector<match_set_t>& matches, const BGraph1& graph1, const BGraph2& graph2,
                                 const SentinelTableau& tableau1, const SentinelTableau& tableau2,
                                 std::unordered_set<std::tuple<size_t, size_t, size_t>>* masked_matches) const;
    
    // split up matches
    template<class BGraph1, class BGraph2>
    std::vector<std::vector<match_set_t>>
    divvy_matches(const std::vector<match_set_t>& matches,
                  const BGraph1& graph1, const BGraph2& graph2,
                  const std::vector<std::pair<SubGraphInfo, SubGraphInfo>>& stitch_graphs,
                  std::vector<std::vector<std::pair<size_t, std::pair<std::vector<size_t>, std::vector<size_t>>>>>& origins_out) const;
    
    std::vector<size_t> assign_reanchor_budget(const std::vector<std::pair<SubGraphInfo, SubGraphInfo>>& stitch_graphs) const;
    
    void merge_fill_in_chains(std::vector<anchor_t>& anchors,
                              const std::vector<std::vector<anchor_t>> stitch_chains,
                              const std::vector<std::pair<SubGraphInfo, SubGraphInfo>>& stitch_graphs,
                              const std::vector<std::vector<std::pair<size_t, std::pair<std::vector<size_t>, std::vector<size_t>>>>>& match_origin) const;
    
    template<class BGraph, class XMerge>
    void fill_in_anchor_chain(std::vector<anchor_t>& anchors,
                              std::vector<match_set_t>& matches,
                              const BGraph& graph1,
                              const BGraph& graph2,
                              const SentinelTableau& tableau1,
                              const SentinelTableau& tableau2,
                              const XMerge& xmerge1,
                              const XMerge& xmerge2,
                              bool restrain_memory,
                              ChainAlgorithm local_chaining_algorithm,
                              double anchor_scale,
                              const std::unordered_set<std::tuple<size_t, size_t, size_t>>* masked_matches) const;
    
    inline void annotate_scores(std::vector<anchor_t>& anchors) const;
    
    inline double anchor_weight(const anchor_t& anchor) const;
    
    // affine edge score, assumes that both pairs of nodes are reachable
    template<class XMerge>
    double edge_weight(uint64_t from_id1, uint64_t to_id1, uint64_t from_id2, uint64_t to_id2,
                       double local_scale,
                       const XMerge& xmerge1, const XMerge& xmerge2,
                       const std::vector<std::vector<size_t>>& switch_dists1,
                       const std::vector<std::vector<size_t>>& switch_dists2) const;
    
    template<class BGraph, class XMerge>
    void instrument_anchor_chain(const std::vector<anchor_t>& chain, double local_scale,
                                 const BGraph& graph1, const BGraph& graph2,
                                 const XMerge& xmerge1, const XMerge& xmerge2) const;
    
    
public:
    
    
    template<class BGraph, class XMerge>
    double estimate_score_scale(std::vector<match_set_t>& matches,
                                const BGraph& graph1, const BGraph& graph2,
                                const SentinelTableau& tableau1, const SentinelTableau& tableau2,
                                const XMerge& xmerge1, const XMerge& xmerge2,
                                bool restrain_memory,
                                std::vector<anchor_t>* chain_out = nullptr,
                                std::unordered_set<std::tuple<size_t, size_t, size_t>>* masked_matches = nullptr) const;
    
};






/*
 * Template implementations
 */


template<class BGraph1, class BGraph2, class XMerge1, class XMerge2>
std::pair<SubGraphInfo, SubGraphInfo>
Extractor::do_extraction(uint64_t from1, uint64_t to1, uint64_t from2, uint64_t to2,
                         const BGraph1& graph1, const BGraph2& graph2,
                         const XMerge1& chain_merge1, const XMerge2& chain_merge2) {
    
    return std::make_pair(extract_connecting_graph(graph1, from1, to1, chain_merge1),
                          extract_connecting_graph(graph2, from2, to2, chain_merge2));
}


template<class BGraph1, class BGraph2, class XMerge1, class XMerge2>
std::vector<std::pair<SubGraphInfo, SubGraphInfo>>
Extractor::extract_graphs_between(const std::vector<anchor_t>& anchor_chain,
                                  const BGraph1& graph1, const BGraph2& graph2,
                                  const SentinelTableau& tableau1, const SentinelTableau& tableau2,
                                  const XMerge1& xmerge1, const XMerge2& xmerge2) {
    
    return extract_graphs_between_internal(anchor_chain, graph1, graph2, &tableau1, &tableau2, xmerge1, xmerge2);
}


template<class BGraph1, class BGraph2, class XMerge1, class XMerge2>
std::vector<std::pair<SubGraphInfo, SubGraphInfo>>
Extractor::extract_graphs_between(const std::vector<anchor_t>& anchor_chain,
                                  const BGraph1& graph1, const BGraph2& graph2,
                                  const XMerge1& xmerge1, const XMerge2& xmerge2) {
    
    return extract_graphs_between_internal(anchor_chain, graph1, graph2, nullptr, nullptr, xmerge1, xmerge2);
}

template<class BGraph1, class BGraph2, class XMerge1, class XMerge2>
std::vector<std::pair<SubGraphInfo, SubGraphInfo>>
Extractor::extract_graphs_between_internal(const std::vector<anchor_t>& anchor_chain,
                                           const BGraph1& graph1, const BGraph2& graph2,
                                           const SentinelTableau* tableau1, const SentinelTableau* tableau2,
                                           const XMerge1& xmerge1, const XMerge2& xmerge2) {
    
    size_t next_log_idx = 0;
    std::vector<size_t> logging_indexes;
    if (logging::level >= logging::Debug) {
        logging_indexes = get_logging_indexes(anchor_chain.size());
        
        logging::log(logging::Debug, "Extracting graphs in between chain of " + std::to_string(anchor_chain.size()) + " anchors");
    }
    
    std::vector<std::pair<SubGraphInfo, SubGraphInfo>> stitch_pairs;
    stitch_pairs.reserve(anchor_chain.size() + 1);
    
    if (anchor_chain.empty() && tableau1 && tableau2) {
        stitch_pairs.emplace_back(do_extraction(tableau1->src_id, tableau1->snk_id, tableau2->src_id, tableau2->snk_id,
                                                graph1, graph2, xmerge1, xmerge2));
    }
    else {
        if (tableau1 && tableau2) {
            stitch_pairs.emplace_back(do_extraction(tableau1->src_id, anchor_chain.front().walk1.front(),
                                                    tableau2->src_id, anchor_chain.front().walk2.front(),
                                                    graph1, graph2, xmerge1, xmerge2));
        }
        
        for (size_t i = 1; i < anchor_chain.size(); ++i) {
            
            if (next_log_idx < logging_indexes.size() && i == logging_indexes[next_log_idx]) {
                logging::log(logging::Debug, "Graph extraction iteration " + std::to_string(i + 1) + " of " + std::to_string(anchor_chain.size()));
                ++next_log_idx;
            }
            
            const auto& prev_anchor = anchor_chain[i - 1];
            const auto& anchor = anchor_chain[i];
            
            stitch_pairs.emplace_back(do_extraction(prev_anchor.walk1.back(), anchor.walk1.front(),
                                                    prev_anchor.walk2.back(), anchor.walk2.front(),
                                                    graph1, graph2, xmerge1, xmerge2));
        }
        
        if (tableau1 && tableau2) {
            stitch_pairs.emplace_back(do_extraction(anchor_chain.back().walk1.back(), tableau1->snk_id,
                                                    anchor_chain.back().walk2.back(), tableau2->snk_id,
                                                    graph1, graph2, xmerge1, xmerge2));
        }
    }
    
    if (logging::level >= logging::Debug) {
        size_t mem = 0;
        for (const auto& stitch_pair : stitch_pairs) {
            mem += stitch_pair.first.memory_size() + stitch_pair.second.memory_size();
        }
        logging::log(logging::Debug, "Extracted subgraphs are occupying " + format_memory_usage(mem) + " of memory.");
    }
    
    return stitch_pairs;
}

// TODO: this is very repetitive with the other one...
template<class BGraph1, class BGraph2, class XMerge1, class XMerge2>
std::pair<std::vector<std::vector<std::pair<SubGraphInfo, SubGraphInfo>>>,
          std::vector<std::pair<SubGraphInfo, SubGraphInfo>>>
Extractor::extract_graphs_between(const std::vector<std::vector<anchor_t>>& anchor_segments,
                                  const BGraph1& graph1, const BGraph2& graph2,
                                  const SentinelTableau& tableau1, const SentinelTableau& tableau2,
                                  const XMerge1& xmerge1, const XMerge2& xmerge2) {
    size_t next_log_idx = 0;
    std::vector<size_t> logging_indexes;
    size_t size = 0;
    if (logging::level >= logging::Debug) {
        for (const auto& segment : anchor_segments) {
            size += segment.size() - 1;
        }
        logging_indexes = get_logging_indexes(size);
        
        logging::log(logging::Debug, "Extracting graphs in between segmented chain of " + std::to_string(size) + " anchors");
    }
    
    std::pair<std::vector<std::vector<std::pair<SubGraphInfo, SubGraphInfo>>>,
              std::vector<std::pair<SubGraphInfo, SubGraphInfo>>> return_val;
    
    auto& within_segment_graphs = return_val.first;
    auto& between_segment_graphs = return_val.second;
    
    size_t idx = 0;
    if (anchor_segments.empty()) {
        // empty, extract whole graph
        between_segment_graphs.emplace_back(do_extraction(tableau1.src_id, tableau1.snk_id, tableau2.src_id, tableau2.snk_id,
                                                          graph1, graph2, xmerge1, xmerge2));
    }
    else {
        
        // before first segment
        between_segment_graphs.emplace_back(do_extraction(tableau1.src_id, anchor_segments.front().front().walk1.front(),
                                                          tableau2.src_id, anchor_segments.front().front().walk2.front(),
                                                          graph1, graph2, xmerge1, xmerge2));
        
        for (size_t i = 0; i < anchor_segments.size(); ++i) {
            
            const auto& segment = anchor_segments[i];
            
            if (i != 0) {
                // between segments
                between_segment_graphs.emplace_back(do_extraction(anchor_segments[i - 1].back().walk1.back(), segment.front().walk1.front(),
                                                                  anchor_segments[i - 1].back().walk2.back(), segment.front().walk2.front(),
                                                                  graph1, graph2, xmerge1, xmerge2));
            }
            
            // within the segment
            within_segment_graphs.emplace_back();
            auto& segment_graphs = within_segment_graphs.back();
            for (size_t j = 1; j < segment.size(); ++j) {
                
                const auto& prev_anchor = segment[j - 1];
                const auto& anchor = segment[j];
                
                segment_graphs.emplace_back(do_extraction(prev_anchor.walk1.back(), anchor.walk1.front(),
                                                          prev_anchor.walk2.back(), anchor.walk2.front(),
                                                          graph1, graph2, xmerge1, xmerge2));
                
                if (next_log_idx < logging_indexes.size() && idx == logging_indexes[next_log_idx]) {
                    logging::log(logging::Debug, "Graph extraction iteration " + std::to_string(idx + 1) + " of " + std::to_string(size));
                    ++next_log_idx;
                }
                ++idx;
            }
        }
        
        // after last segment
        between_segment_graphs.emplace_back(do_extraction(anchor_segments.back().back().walk1.back(), tableau1.snk_id,
                                                          anchor_segments.back().back().walk2.back(), tableau2.snk_id,
                                                          graph1, graph2, xmerge1, xmerge2));
    }
    
    if (logging::level >= logging::Debug) {
        size_t within_mem = 0;
        for (const auto& segment : within_segment_graphs) {
            for (const auto& subgraph_pair : segment) {
                within_mem += subgraph_pair.first.memory_size() + subgraph_pair.second.memory_size();
            }
        }
        size_t between_mem = 0;
        for (const auto& subgraph_pair : between_segment_graphs) {
            between_mem += subgraph_pair.first.memory_size() + subgraph_pair.second.memory_size();
        }
        logging::log(logging::Debug, "Extracted subgraphs within segments are occupying " + format_memory_usage(within_mem) + " of memory.");
        logging::log(logging::Debug, "Extracted subgraphs between segments are occupying " + format_memory_usage(between_mem) + " of memory.");
    }
    
    return return_val;
}


template<class BGraph1, class BGraph2>
void Extractor::project_paths(const BGraph1& graph1, const BGraph2& graph2,
                              std::vector<std::pair<SubGraphInfo, SubGraphInfo>>& stitch_graphs) {
    
    StepIndex step_index1(graph1);
    StepIndex step_index2(graph2);
    for (auto& stitch_pair : stitch_graphs) {
        do_project(graph1, stitch_pair.first, step_index1);
        do_project(graph2, stitch_pair.second, step_index2);
    }
}


template<class BGraph>
void Extractor::do_project(const BGraph& graph, SubGraphInfo& subgraph,
                           const StepIndex& step_index) {
    
    std::unordered_map<uint64_t, uint64_t> path_ids;
    for (auto node_id : topological_order(subgraph.subgraph)) {
        for (const auto& step : step_index.path_steps(subgraph.back_translation[node_id])) {
            auto it = path_ids.find(step.first);
            if (it == path_ids.end()) {
                auto new_path_id = subgraph.subgraph.add_path(graph.path_name(step.first));
                it = path_ids.insert(std::make_pair(step.first, new_path_id)).first;
            }
            subgraph.subgraph.extend_path(it->second, node_id);
        }
    }
}


template<class BGraph, class XMerge>
void Anchorer::fill_in_anchor_chain(std::vector<anchor_t>& anchors,
                                    std::vector<match_set_t>& matches,
                                    const BGraph& graph1,
                                    const BGraph& graph2,
                                    const SentinelTableau& tableau1,
                                    const SentinelTableau& tableau2,
                                    const XMerge& xmerge1,
                                    const XMerge& xmerge2,
                                    bool restrain_memory,
                                    ChainAlgorithm local_chaining_algorithm,
                                    double anchor_scale,
                                    const std::unordered_set<std::tuple<size_t, size_t, size_t>>* masked_matches) const {
    
    if (anchors.empty()) {
        logging::log(logging::Debug, "Skipping fill-in anchoring on an empty chain");
        return;
    }
    
    logging::log(logging::Debug, "Extracting subgraphs for fill-in anchoring");
    
    auto fill_in_graphs = extract_graphs_between(anchors, graph1, graph2,
                                                 tableau1, tableau2,
                                                 xmerge1, xmerge2);
    
    project_paths(graph1, graph2, fill_in_graphs);
    
    logging::log(logging::Debug, "Assigning matches to fill-in subproblems");
    
    std::vector<std::vector<std::pair<size_t, std::pair<std::vector<size_t>, std::vector<size_t>>>>> match_origin;
    auto fill_in_matches = divvy_matches(matches, graph1, graph2, fill_in_graphs, match_origin);
    
    auto budgets = assign_reanchor_budget(fill_in_graphs);
    
    logging::log(logging::Debug, "Performing anchoring subproblems for fill-in anchoring");
    
    std::vector<std::vector<anchor_t>> fill_in_anchors(fill_in_graphs.size());
    
    for (size_t i = 0; i < fill_in_graphs.size(); ++i) {
        
        XMerge fill_in_xmerge1(fill_in_graphs[i].first.subgraph);
        XMerge fill_in_xmerge2(fill_in_graphs[i].second.subgraph);
        
        std::unordered_set<std::tuple<size_t, size_t, size_t>> fill_in_masked_matches;
        if (masked_matches) {
            // translate the masked indexes for the new set
            const auto& fill_in_origin = match_origin[i];
            for (size_t set = 0; set < fill_in_origin.size(); ++set) {
                size_t orig_set = fill_in_origin[set].first;
                const auto& walks = fill_in_origin[set].second;
                for (size_t idx1 = 0; idx1 < walks.first.size(); ++idx1) {
                    size_t orig_idx1 = walks.first[idx1];
                    for (size_t idx2 = 0; idx2 < walks.second.size(); ++idx2) {
                        size_t orig_idx2 = walks.second[idx2];
                        if (masked_matches->count(std::make_tuple(orig_set, orig_idx1, orig_idx2))) {
                            fill_in_masked_matches.emplace(set, idx1, idx2);
                        }
                    }
                }
            }
        }
        
        fill_in_anchors[i] = std::move(anchor_chain(fill_in_matches[i],
                                                    fill_in_graphs[i].first.subgraph,
                                                    fill_in_graphs[i].second.subgraph,
                                                    fill_in_xmerge1, fill_in_xmerge2,
                                                    restrain_memory,
                                                    &fill_in_graphs[i].first.sources,
                                                    &fill_in_graphs[i].second.sources,
                                                    &fill_in_graphs[i].first.sinks,
                                                    &fill_in_graphs[i].second.sinks,
                                                    budgets[i], true,
                                                    local_chaining_algorithm, anchor_scale,
                                                    &fill_in_masked_matches));
    }
    
    merge_fill_in_chains(anchors, fill_in_anchors, fill_in_graphs, match_origin);
    
    logging::log(logging::Debug, "Filled-in anchor chain consists of " + std::to_string(anchors.size()) + " anchors");
        
}

template<class BGraph1, class BGraph2>
std::vector<std::vector<match_set_t>>
Anchorer::divvy_matches(const std::vector<match_set_t>& matches,
                        const BGraph1& graph1, const BGraph2& graph2,
                        const std::vector<std::pair<SubGraphInfo, SubGraphInfo>>& fill_in_graphs,
                        std::vector<std::vector<std::pair<size_t, std::pair<std::vector<size_t>, std::vector<size_t>>>>>& origins_out) const {
    
    static const bool debug = false;
    
    // identify which nodes are in the stitch graphs
    std::vector<std::pair<size_t, uint64_t>> forward_trans1(graph1.node_size(), std::pair<size_t, uint64_t>(-1, -1));
    std::vector<std::pair<size_t, uint64_t>> forward_trans2(graph2.node_size(), std::pair<size_t, uint64_t>(-1, -1));
    for (size_t i = 0; i < fill_in_graphs.size(); ++i) {
        const auto& back_trans1 = fill_in_graphs[i].first.back_translation;
        for (uint64_t fwd_id1 = 0; fwd_id1 < back_trans1.size(); ++fwd_id1) {
            forward_trans1[back_trans1[fwd_id1]] = std::make_pair(i, fwd_id1);
        }
        const auto& back_trans2 = fill_in_graphs[i].second.back_translation;
        for (uint64_t fwd_id2 = 0; fwd_id2 < back_trans2.size(); ++fwd_id2) {
            forward_trans2[back_trans2[fwd_id2]] = std::make_pair(i, fwd_id2);
        }
    }
    
    if (debug) {
        std::cerr << "forward translations to stitch graphs:\n";
        std::cerr << "graph1:\n";
        for (uint64_t n = 0; n < graph1.node_size(); ++n) {
            std::cerr << n << ":\t" << forward_trans1[n].first << '\t' << forward_trans1[n].second << '\n';
        }
        std::cerr << "graph2:\n";
        for (uint64_t n = 0; n < graph2.node_size(); ++n) {
            std::cerr << n << ":\t" << forward_trans2[n].first << '\t' << forward_trans2[n].second << '\n';
        }
    }
    
    std::vector<std::vector<match_set_t>> divvied(fill_in_graphs.size());
    origins_out.clear();
    origins_out.resize(fill_in_graphs.size());
    
    
    for (size_t i = 0; i < matches.size(); ++i) {
        const auto& match_set = matches[i];
        std::unordered_set<size_t> stitch_sets_initialized;
        for (size_t j = 0; j < match_set.walks1.size(); ++j) {
            const auto& walk1 = match_set.walks1[j];
            size_t stitch_idx = forward_trans1[walk1.front()].first;
            if (stitch_idx != -1 && stitch_idx == forward_trans1[walk1.back()].first) {
                // both ends of the path are in the same stitch graph
                
                if (!stitch_sets_initialized.count(stitch_idx)) {
                    // we haven't initialized a corresponding match set for the stitch graph yet
                    divvied[stitch_idx].emplace_back();
                    origins_out[stitch_idx].emplace_back();
                    origins_out[stitch_idx].back().first = i;
                    // the counts and length get retained from the original match sets
                    divvied[stitch_idx].back().count1 = match_set.count1;
                    divvied[stitch_idx].back().count2 = match_set.count2;
                    divvied[stitch_idx].back().full_length = match_set.full_length;
                    stitch_sets_initialized.insert(stitch_idx);
                }
                
                // copy and translate the walk
                origins_out[stitch_idx].back().second.first.push_back(j);
                divvied[stitch_idx].back().walks1.emplace_back();
                auto& stitch_walk = divvied[stitch_idx].back().walks1.back();
                for (auto node_id : walk1) {
                    stitch_walk.push_back(forward_trans1[node_id].second);
                }
            }
        }
        
        for (size_t k = 0; k < match_set.walks2.size(); ++k) {
            const auto& walk2 = match_set.walks2[k];
            size_t stitch_idx = forward_trans2[walk2.front()].first;
            if (stitch_sets_initialized.count(stitch_idx) && stitch_idx == forward_trans2[walk2.back()].first) {
                // both ends of the path are in the same stitch graph
                // note: checking whether it's initialized also checks whether it's -1 because
                // of the previous loop
                
                // copy and translate the walk
                origins_out[stitch_idx].back().second.second.push_back(k);
                divvied[stitch_idx].back().walks2.emplace_back();
                auto& stitch_walk = divvied[stitch_idx].back().walks2.back();
                for (auto node_id : walk2) {
                    stitch_walk.push_back(forward_trans2[node_id].second);
                }
            }
        }
        
        // clear out any initialized sets that didn't end up getting a walk from graph2
        for (auto stitch_idx : stitch_sets_initialized) {
            if (divvied[stitch_idx].back().walks2.empty()) {
                origins_out[stitch_idx].pop_back();
                divvied[stitch_idx].pop_back();
            }
        }
    }
    
    return divvied;
}


template<class BGraph1, class BGraph2>
void Anchorer::split_branching_matches(std::vector<match_set_t>& matches, const BGraph1& graph1, const BGraph2& graph2,
                                       const SentinelTableau& tableau1, const SentinelTableau& tableau2,
                                       std::unordered_set<std::tuple<size_t, size_t, size_t>>* masked_matches) const {
    
    // handle this as a special case so we can assume it's positive later
    if (anchor_split_limit == 0) {
        return;
    }
    
    logging::log(logging::Verbose, "Splitting matches at branch points.");
    
    
    std::unordered_map<size_t, std::vector<std::pair<size_t, size_t>>> set_masked_matches;
    if (masked_matches) {
        for (auto mask : *masked_matches) {
            set_masked_matches[std::get<0>(mask)].emplace_back(std::get<1>(mask), std::get<2>(mask));
        }
    }
    
    SuperbubbleTree bubbles1(graph1, tableau1), bubbles2(graph2, tableau2);
    SuperbubbleDistances bub_dists1(bubbles1, graph1), bub_dists2(bubbles2, graph2);
    
    size_t num_original_sets = matches.size();
    size_t num_new_split_pairs = 0;
    size_t num_original_pairs = 0;
    for (size_t i = 0; i < num_original_sets; ++i) {
        
        // don't split short anchors or anchors from very large match sets
        if (matches[i].walks1.size() * matches[i].walks2.size() > max_split_match_set_size || matches[i].walks1.front().size() <  min_split_length) {
            continue;
        }
        
        // find the indexes where one of the matches is at a branch
        std::vector<size_t> division_idxs;
        {
            const auto& match_set = matches[i];
            // note: we stop early, since a branch after the final position is not problematic
            for (size_t j = 0; j < match_set.walks1.front().size(); ++j) {
                
                if (j == anchor_split_limit && j + anchor_split_limit < match_set.walks1.front().size()) {
                    // skip to the suffix of the match
                    j = match_set.walks1.front().size() - anchor_split_limit;
                }
                
                // don't make two splits at the same index
                if (j != 0 && (division_idxs.empty() || division_idxs.back() != j)) {
                    // look for splits facing backwards
                    bool found_branch = false;
                    for (size_t k = 0; k < match_set.walks1.size() && !found_branch; ++k) {
                        auto bub_id = bubbles1.structure_ending_at(match_set.walks1[k][j]);
                        if (bub_id != -1) {
                            auto mm_dist = bub_dists1.structure_min_max_dist(bub_id);
                            if (mm_dist.second - mm_dist.first >= min_path_length_spread) {
                                found_branch = true;
                            }
                        }
                    }
                    for (size_t k = 0; k < match_set.walks2.size() && !found_branch; ++k) {
                        auto bub_id = bubbles2.structure_ending_at(match_set.walks2[k][j]);
                        if (bub_id != -1) {
                            auto mm_dist = bub_dists2.structure_min_max_dist(bub_id);
                            if (mm_dist.second - mm_dist.first >= min_path_length_spread) {
                                found_branch = true;
                            }
                        }
                    }
                    
                    if (found_branch) {
                        division_idxs.push_back(j);
                    }
                }
                
                // TODO: repetitive
                if (j + 1 != match_set.walks1.front().size()) {
                    // look for splits facing forwards
                    bool found_branch = false;
                    for (size_t k = 0; k < match_set.walks1.size() && !found_branch; ++k) {
                        auto bub_id = bubbles1.structure_beginning_at(match_set.walks1[k][j]);
                        if (bub_id != -1) {
                            auto mm_dist = bub_dists1.structure_min_max_dist(bub_id);
                            if (mm_dist.second - mm_dist.first >= min_path_length_spread) {
                                found_branch = true;
                            }
                        }
                    }
                    for (size_t k = 0; k < match_set.walks2.size() && !found_branch; ++k) {
                        auto bub_id = bubbles2.structure_beginning_at(match_set.walks2[k][j]);
                        if (bub_id != -1) {
                            auto mm_dist = bub_dists2.structure_min_max_dist(bub_id);
                            if (mm_dist.second - mm_dist.first >= min_path_length_spread) {
                                found_branch = true;
                            }
                        }
                    }
                    
                    if (found_branch) {
                        division_idxs.push_back(j + 1);
                    }
                }
            }
        }
        
        if (!division_idxs.empty()) {
            // create new match sets for the divisions
            size_t end = matches[i].walks1.front().size();
            for (size_t j = division_idxs.size() - 1; j < division_idxs.size(); --j) {
                
                if (masked_matches) {
                    auto it = set_masked_matches.find(i);
                    if (it != set_masked_matches.end()) {
                        for (auto mask : it->second) {
                            masked_matches->emplace(matches.size(), mask.first, mask.second);
                        }
                    }
                }
                
                size_t idx = division_idxs[j];
                matches.emplace_back();
                const auto& match_set = matches[i];
                auto& split_match_set = matches.back();
                
                // copy the partial walks
                split_match_set.walks1.reserve(match_set.walks1.size());
                for (const auto& walk1 : match_set.walks1) {
                    split_match_set.walks1.emplace_back(walk1.begin() + idx, walk1.begin() + end);
                }
                split_match_set.walks2.reserve(match_set.walks2.size());
                for (const auto& walk2 : match_set.walks2) {
                    split_match_set.walks2.emplace_back(walk2.begin() + idx, walk2.begin() + end);
                }
                // copy the annotations
                split_match_set.count1 = match_set.count1;
                split_match_set.count2 = match_set.count2;
                split_match_set.full_length = match_set.full_length;
                
                num_new_split_pairs += split_match_set.walks1.size() * split_match_set.walks2.size();
                
                end = idx;
            }
            
            // shorten the original match set to create the final division
            auto& match_set = matches[i];
            for (auto& walk1 : match_set.walks1) {
                walk1.resize(division_idxs.front());
            }
            for (auto& walk2 : match_set.walks2) {
                walk2.resize(division_idxs.front());
            }
        }
    }
    
    logging::log(logging::Debug, "Add " + std::to_string(matches.size() - num_original_sets) + " match sets containing " + std::to_string(num_new_split_pairs) + " match pairs as a result of splitting.");
}

template<class BGraph, class XMerge>
std::vector<anchor_t> Anchorer::anchor_chain(std::vector<match_set_t>& matches,
                                             const BGraph& graph1,
                                             const BGraph& graph2,
                                             const SentinelTableau& tableau1,
                                             const SentinelTableau& tableau2,
                                             const XMerge& xmerge1,
                                             const XMerge& xmerge2,
                                             bool restrain_memory,
                                             std::unordered_set<std::tuple<size_t, size_t, size_t>>* masked_matches,
                                             double* override_scale) const {
    
    if (split_matches_at_branchpoints) {
        split_branching_matches(matches, graph1, graph2, tableau1, tableau2, masked_matches);
    }
    
    double scale = 1.0;
    if (override_scale) {
        scale = *override_scale;
    }
    else if (chaining_algorithm == SparseAffine && autocalibrate_gap_penalties) {
        // this is only to adjust gap penalties, so don't bother if we're not using them
        logging::log(logging::Verbose, "Calibrating gap penalties.");
        scale = estimate_score_scale(matches, graph1, graph2, tableau1, tableau2, xmerge1, xmerge2, restrain_memory, nullptr, masked_matches);
        logging::log(logging::Debug, "Estimated score scale: " + std::to_string(scale));
    }
    auto anchors = anchor_chain(matches, graph1, graph2, tableau1, tableau2, xmerge1, xmerge2,
                                restrain_memory, chaining_algorithm, false, scale, masked_matches);
    
    static const bool instrument = false;
    if (instrument) {
        instrument_anchor_chain(anchors, scale, graph1, graph2, xmerge1, xmerge2);
    }
    
    return anchors;
}

template<class BGraph, class XMerge>
double Anchorer::estimate_score_scale(std::vector<match_set_t>& matches,
                                      const BGraph& graph1, const BGraph& graph2,
                                      const SentinelTableau& tableau1, const SentinelTableau& tableau2,
                                      const XMerge& xmerge1, const XMerge& xmerge2,
                                      bool restrain_memory,
                                      std::vector<anchor_t>* chain_out,
                                      std::unordered_set<std::tuple<size_t, size_t, size_t>>* masked_matches) const {
    
    // get an anchoring with unscored gaps
    // FIXME: should i handle masked matches here?
    auto anchors = anchor_chain(matches, graph1, graph2, tableau1, tableau2,
                                xmerge1, xmerge2, restrain_memory, Sparse, true, 1.0, masked_matches);
    
    // measure its weight
    double total_weight = 0.0;
    size_t total_length = 0;
    for (const auto& anchor : anchors) {
        total_weight += anchor_weight(anchor);
        total_length += anchor.walk1.size();
    }
    
    // estimate the length of the in-between bits
    auto fill_in_graphs = extract_graphs_between(anchors, graph1, graph2,
                                                 tableau1, tableau2,
                                                 xmerge1, xmerge2);
    
    for (const auto& fill_in_pair : fill_in_graphs) {
        
        size_t fill_in_length = std::numeric_limits<size_t>::max();
        for (auto subgraph_ptr : {&fill_in_pair.first, &fill_in_pair.second}) {
            // compute the minimum distance across either of the stitch graphs
            const auto& subgraph = *subgraph_ptr;
            if (subgraph.subgraph.node_size() == 0) {
                fill_in_length = 0;
            }
            else {
                fill_in_length = std::min<size_t>(fill_in_length, source_sink_minmax(subgraph).first);
            }
        }
        
        total_length += fill_in_length;
    }
    
    if (chain_out) {
        *chain_out = std::move(anchors);
    }
    
    return total_weight / total_length;
}


template<class BGraph, class XMerge>
std::vector<anchor_t> Anchorer::anchor_chain(std::vector<match_set_t>& matches,
                                             const BGraph& graph1,
                                             const BGraph& graph2,
                                             const SentinelTableau& tableau1,
                                             const SentinelTableau& tableau2,
                                             const XMerge& xmerge1,
                                             const XMerge& xmerge2,
                                             bool restrain_memory,
                                             ChainAlgorithm local_chaining_algorithm,
                                             bool suppress_verbose_logging,
                                             double anchor_scale,
                                             std::unordered_set<std::tuple<size_t, size_t, size_t>>* masked_matches) const {
    
    // adjust the number of matches we're going to look at downward if the alignment isn't promising
    size_t adjusted_max_num_match_pairs = std::min<size_t>(round((anchor_scale / score_function->score_scale) * max_num_match_pairs),
                                                           max_num_match_pairs);

    std::vector<anchor_t> anchors;
    if (global_anchoring) {
        anchors = std::move(anchor_chain(matches, graph1, graph2, xmerge1, xmerge2, restrain_memory,
                                         &graph1.next(tableau1.src_id), &graph2.next(tableau2.src_id),
                                         &graph1.previous(tableau1.snk_id), &graph2.previous(tableau2.snk_id),
                                         adjusted_max_num_match_pairs, suppress_verbose_logging, local_chaining_algorithm, anchor_scale,
                                         masked_matches));
    }
    else {
        anchors = std::move(anchor_chain(matches, graph1, graph2, xmerge1, xmerge2, restrain_memory,
                                         nullptr, nullptr, nullptr, nullptr,
                                         adjusted_max_num_match_pairs, suppress_verbose_logging, local_chaining_algorithm, anchor_scale,
                                         masked_matches));
    }
    
    if (do_fill_in_anchoring) {
        fill_in_anchor_chain(anchors, matches, graph1, graph2, tableau1, tableau2, xmerge1, xmerge2, restrain_memory,
                             local_chaining_algorithm, anchor_scale, masked_matches);
    }
    
    return anchors;
}

template<class BGraph, class XMerge>
std::vector<anchor_t> Anchorer::anchor_chain(std::vector<match_set_t>& matches,
                                             const BGraph& graph1,
                                             const BGraph& graph2,
                                             const XMerge& chain_merge1,
                                             const XMerge& chain_merge2,
                                             bool restrain_memory,
                                             const std::vector<uint64_t>* sources1,
                                             const std::vector<uint64_t>* sources2,
                                             const std::vector<uint64_t>* sinks1,
                                             const std::vector<uint64_t>* sinks2,
                                             size_t local_max_num_match_pairs,
                                             bool suppress_verbose_logging,
                                             ChainAlgorithm local_chaining_algorithm,
                                             double anchor_scale,
                                             std::unordered_set<std::tuple<size_t, size_t, size_t>>* masked_matches) const {
    
    size_t total_num_pairs = 0;
    for (const auto& match_set : matches) {
        total_num_pairs += match_set.walks1.size() * match_set.walks2.size();
    }
    
    size_t removed = 0;
    size_t max_match_size = 0;
    size_t pairs_left = local_max_num_match_pairs;
    
    if (total_num_pairs <= local_max_num_match_pairs) {
        // we can just take them all (and update the bookkeeping)
        pairs_left = local_max_num_match_pairs - total_num_pairs;
    }
    else {
        // we need to limit the number of nodes
        
        if (!suppress_verbose_logging) {
            logging::log(logging::Debug, "Selecting a maximum of " + std::to_string(local_max_num_match_pairs) + " matches for anchoring");
        }
        
        // prioritize based on the score
        // note: we prioritize based on the weight of the full length anchor in the case of shortened anchors
        // TODO: adjust this by masked match count?
        auto order = range_vector(matches.size());
        std::stable_sort(order.begin(), order.end(), [&](size_t i, size_t j) {
            return (score_function->anchor_weight(matches[i].count1, matches[i].count2, matches[i].full_length) >
                    score_function->anchor_weight(matches[j].count1, matches[j].count2, matches[j].full_length));
        });
        
        // greedily choose matches as long as we have budget left
        for (size_t i = 0; i < order.size(); ++i) {
            auto& match = matches[order[i]];
            if (score_function->anchor_weight(match.count1, match.count2, match.walks1.front().size()) < 0.0) {
                // these are sorted by score, so nothing else will have positive score
                removed += (order.size() - i);
                break;
            }
            size_t pair_count = match.walks1.size() * match.walks2.size();
            if (pairs_left >= pair_count) {
                max_match_size = std::max(max_match_size, std::max(match.walks1.size(), match.walks2.size()));
                pairs_left -= pair_count;
                std::swap(order[i - removed], order[i]);
            }
            else {
                // remove this match
                ++removed;
            }
        }
        
        // actually put the matches in this order
        auto index = invert(order);
        if (masked_matches) {
            // update the match set indexes for the new ordering (have to recreate because modifying keys)
            std::unordered_set<std::tuple<size_t, size_t, size_t>> reindexed_masked_matches;
            reindexed_masked_matches.reserve(masked_matches->size());
            for (const auto& mask : *masked_matches) {
                reindexed_masked_matches.emplace(index[std::get<0>(mask)], std::get<1>(mask), std::get<2>(mask));
            }
            std::swap(reindexed_masked_matches, *masked_matches);
        }
        reorder(matches, index);
        
        if (debug_anchorer) {
            std::cerr << "moved " << removed << " unique anchor sequences to back limit to " << local_max_num_match_pairs << " total pairs\n";
        }
    }
    
    bool switch_graphs = (graph1.node_size() * chain_merge1.chain_size() > graph2.node_size() * chain_merge2.chain_size());
    std::unique_ptr<std::unordered_set<std::tuple<size_t, size_t, size_t>>> switched_masked_matches(nullptr);
    if (switch_graphs) {
        // we will use less memory if graph1 and graph2 are switched
        for (auto& match_set : matches) {
            std::swap(match_set.walks1, match_set.walks2);
            std::swap(match_set.count1, match_set.count2);
        }
        if (masked_matches) {
            // we also need to swap any masked matches
            switched_masked_matches.reset(new std::unordered_set<std::tuple<size_t, size_t, size_t>>());
            switched_masked_matches->reserve(masked_matches->size());
            for (const auto& mask : *masked_matches) {
                switched_masked_matches->emplace(std::get<0>(mask), std::get<2>(mask), std::get<1>(mask));
            }
            masked_matches = switched_masked_matches.get();
        }
    }
    
    // compute bounds on the size of integer variables
    size_t num_match_sets = matches.size() - removed;
    size_t max_chain_size = std::max(chain_merge1.chain_size(), chain_merge2.chain_size());
    // TODO: are there tighter bounds possible based on the longest paths? would need to factor in post switch distance too...
    size_t max_dist = std::max(graph1.node_size(), graph2.node_size());
    size_t max_diag_diff = graph1.node_size() + graph2.node_size();
    size_t num_anchors = local_max_num_match_pairs - pairs_left;
    
    // choose the switched or unswitched arguments for the internal algorithms
    const auto& graph1_arg = switch_graphs ? graph2 : graph1;
    const auto& graph2_arg = switch_graphs ? graph1 : graph2;
    const auto& chain_merge1_arg = switch_graphs ? chain_merge2 : chain_merge1;
    const auto& chain_merge2_arg = switch_graphs ? chain_merge1 : chain_merge2;
    const auto sources1_arg = switch_graphs ? sources2 : sources1;
    const auto sources2_arg = switch_graphs ? sources1 : sources2;
    const auto sinks1_arg = switch_graphs ? sinks2 : sinks1;
    const auto sinks2_arg = switch_graphs ? sinks1 : sinks2;
    
    // macros to generate the calls for the template functions/
    #define _gen_sparse_affine(UIntSet, UIntMatch, UIntDist, IntShift, UIntAnchor, ShiftMatchVec, DistMatchVec, DistVec, AnchorVec, MBank) \
        if (logging::level >= logging::Debug && !suppress_verbose_logging) { \
            logging::log(logging::Debug, std::string("Integer widths: set ") + #UIntSet + ", match " + #UIntMatch + ", dist " + #UIntDist + ", shift " + #IntShift);\
        }\
        chain = std::move(sparse_affine_chain_dp<UIntSet, UIntMatch, UIntDist, IntShift, UIntAnchor, float, ShiftMatchVec, DistMatchVec, DistVec, AnchorVec, MBank> \
                                                (matches, graph1_arg, graph2_arg, chain_merge1_arg, chain_merge2_arg, \
                                                 gap_open, gap_extend, anchor_scale, num_match_sets, suppress_verbose_logging, \
                                                 sources1_arg, sources2_arg, sinks1_arg, sinks2_arg, masked_matches))
    
    #define _gen_sparse(UIntDist, UIntSet, UIntMatch, UIntAnchor, DisMatchVec, AnchorVec, MBank) \
        chain = std::move(sparse_chain_dp<UIntDist, UIntSet, UIntMatch, UIntAnchor, float, std::vector<std::pair<UIntDist, MBank::match_id_t>>, std::vector<UIntAnchor>, MBank> \
                                         (matches, graph1_arg, chain_merge1_arg, chain_merge2_arg, \
                                          num_match_sets, suppress_verbose_logging, \
                                          sources1_arg, sources2_arg, sinks1_arg, sinks2_arg, \
                                          masked_matches))
    
    #define _gen_exhaustive(UIntSet, UIntMatch) \
        chain = std::move(exhaustive_chain_dp<UIntSet, UIntMatch>(matches, graph1_arg, graph2_arg, chain_merge1_arg, chain_merge2_arg, false, \
                                                                  anchor_scale, num_match_sets, sources1_arg, sources2_arg, sinks1_arg, sinks2_arg, \
                                                                  masked_matches))
    
    // we have to give these aliases so that the macros don't interpret the commas as separate arguments
    using SmallMatchBank = MatchBank<uint32_t, uint16_t, float>;
    using LargeMatchBank = MatchBank<uint64_t, uint32_t, float>;
    using SmallPackedMatchBank = PackedMatchBank<uint32_t, float>;
    using LargePackedMatchBank = PackedMatchBank<uint64_t, float>;
    
    using SmallShiftMatchVector = std::vector<std::pair<int32_t, SmallMatchBank::match_id_t>>;
    using LargeShiftMatchVector = std::vector<std::pair<int64_t, LargeMatchBank::match_id_t>>;
    using SmallPackedShiftMatchVector = VectorPair<SignedPackedVector, PackedVector, int32_t, SmallPackedMatchBank::match_id_t>;
    using LargePackedShiftMatchVector = VectorPair<SignedPackedVector, PackedVector, int64_t, LargePackedMatchBank::match_id_t>;
    
    using SmallDistMatchVector = std::vector<std::pair<uint32_t, SmallMatchBank::match_id_t>>;
    using LargeDistMatchVector = std::vector<std::pair<uint64_t, LargeMatchBank::match_id_t>>;
    using SmallPackedDistMatchVector = VectorPair<PackedVector, PackedVector, uint32_t, SmallPackedMatchBank::match_id_t>;
    using LargePackedDistMatchVector = VectorPair<PackedVector, PackedVector, uint64_t, LargePackedMatchBank::match_id_t>;
    // compute the optimal chain using DP
    std::vector<anchor_t> chain;
    switch (local_chaining_algorithm) {
        case SparseAffine:
            // we limit the number of cases dramatically here to avoid excessive compile time
            // TODO: add them back in?
            if (restrain_memory) {
                // use bit-packed data structures to reduce memory use
                if (num_anchors < std::numeric_limits<uint32_t>::max()
                    && max_diag_diff < std::numeric_limits<int32_t>::max()) {
                    
                    _gen_sparse_affine(uint32_t, uint16_t, uint32_t, int32_t, uint32_t,
                                       SmallPackedShiftMatchVector, SmallPackedDistMatchVector, PackedVector, PackedVector, SmallPackedMatchBank);
                }
                else {
                    _gen_sparse_affine(uint64_t, uint32_t, uint64_t, int64_t, uint64_t,
                                       LargePackedShiftMatchVector, LargePackedDistMatchVector, PackedVector, PackedVector, LargePackedMatchBank);
                }
            }
            else if (num_match_sets < std::numeric_limits<uint32_t>::max()
                     && max_match_size < std::numeric_limits<uint16_t>::max()
                     && max_diag_diff < std::numeric_limits<int32_t>::max()
                     && num_anchors < std::numeric_limits<uint32_t>::max()) {
                
                _gen_sparse_affine(uint32_t, uint16_t, uint32_t, int32_t, uint32_t,
                                   SmallShiftMatchVector, SmallDistMatchVector, std::vector<uint32_t>, std::vector<uint32_t>, SmallMatchBank);
            }
            else {
                _gen_sparse_affine(uint64_t, uint32_t, uint64_t, int64_t, uint64_t,
                                   LargeShiftMatchVector, LargeDistMatchVector, std::vector<uint64_t>, std::vector<uint64_t>, LargeMatchBank);
            }
            break;
            
        case Sparse:
            if (max_dist < std::numeric_limits<uint32_t>::max()
                && num_match_sets < std::numeric_limits<uint32_t>::max()
                && max_match_size < std::numeric_limits<uint16_t>::max()
                && num_anchors < std::numeric_limits<uint32_t>::max()) {
                
                _gen_sparse(uint32_t, uint32_t, uint16_t, uint32_t, SmallDistMatchVector, std::vector<uint32_t>, SmallMatchBank);
            }
            else {
                _gen_sparse(uint64_t, uint64_t, uint32_t, uint64_t, LargeDistMatchVector, std::vector<uint64_t>, LargeMatchBank);
            }
            // TODO: I could add more cases here, but this isn't currently the memory bottleneck
            break;
            
        case Exhaustive:
            _gen_exhaustive(uint64_t, uint32_t);
            // TODO: killed these cases to avoid excessive compile time on GCC
            break;
            
        default:
            throw std::runtime_error("Unrecognized chaining algorithm: " + std::to_string((int) local_chaining_algorithm));
            break;
    }
    
    #undef _gen_sparse_affine
    #undef _gen_sparse
    #undef _gen_exhaustive
    
    if (switch_graphs) {
        // return the matches and chain to the original order
        for (auto& match_set : matches) {
            std::swap(match_set.walks1, match_set.walks2);
            std::swap(match_set.count1, match_set.count2);
        }
        for (auto& anchor : chain) {
            std::swap(anchor.walk1, anchor.walk2);
            std::swap(anchor.count1, anchor.count2);
            std::swap(anchor.idx1, anchor.idx2);
            anchor.gap_before = -anchor.gap_before;
            anchor.gap_after = -anchor.gap_after;
        }
    }
    
    return chain;
}

inline void Anchorer::annotate_scores(std::vector<anchor_t>& anchors) const {
    for (auto& anchor : anchors) {
        anchor.score = anchor_weight(anchor);
    }
}

inline double Anchorer::anchor_weight(const anchor_t& anchor) const {
    return score_function->anchor_weight(anchor.count1, anchor.count2, anchor.walk1.size(), anchor.full_length);
}


template <typename UIntSet, typename UIntMatch, class BGraph, class XMerge>
std::vector<anchor_t> Anchorer::exhaustive_chain_dp(const std::vector<match_set_t>& match_sets,
                                                    const BGraph& graph1,
                                                    const BGraph& graph2,
                                                    const XMerge& chain_merge1,
                                                    const XMerge& chain_merge2,
                                                    bool score_edges,
                                                    double local_scale,
                                                    size_t num_match_sets,
                                                    const std::vector<uint64_t>* sources1,
                                                    const std::vector<uint64_t>* sources2,
                                                    const std::vector<uint64_t>* sinks1,
                                                    const std::vector<uint64_t>* sinks2,
                                                    const std::unordered_set<std::tuple<size_t, size_t, size_t>>* masked_matches) const {
            
    if (debug_anchorer) {
        std::cerr << "beginning exhaustive DP chaining algorithm\n";
    }
    
    assert((sources1 == nullptr) == (sources2 == nullptr));
    assert((sinks1 == nullptr) == (sinks2 == nullptr));
    
    std::vector<std::vector<size_t>> switch_dists1, switch_dists2;
    if (score_edges) {
        switch_dists1 = std::move(post_switch_distances<size_t>(graph1, chain_merge1));
        switch_dists2 = std::move(post_switch_distances<size_t>(graph2, chain_merge2));
    }
    
    // make a graph of match nodes
    AnchorGraph anchor_graph;
    for (size_t i = 0; i < num_match_sets; ++i) {
        
        const auto& match_set = match_sets[i];
        
        double weight = score_function->anchor_weight(match_set.count1, match_set.count2,
                                                      match_set.walks1.front().size(), match_set.full_length);
        
        for (size_t idx1 = 0; idx1 < match_set.walks1.size(); ++idx1) {
            for (size_t idx2 = 0; idx2 < match_set.walks2.size(); ++idx2) {
                if (masked_matches && masked_matches->count(std::tuple<size_t, size_t, size_t>(i, idx1, idx2))) {
                    continue;
                }
                double initial_weight = 0.0;
                double final_weight = 0.0;
                if (sources1) {
                    initial_weight = std::numeric_limits<double>::lowest();
                    for (auto src_id1 : *sources1) {
                        for (auto src_id2 : *sources2) {
                            if ((chain_merge1.reachable(src_id1, match_set.walks1[idx1].front()) ||
                                 src_id1 == match_set.walks1[idx1].front()) &&
                                (chain_merge2.reachable(src_id2, match_set.walks2[idx2].front()) ||
                                 src_id2 == match_set.walks2[idx2].front())) {
                                
                                if (score_edges) {
                                    initial_weight = std::max(initial_weight, edge_weight(src_id1, match_set.walks1[idx1].front(),
                                                                                          src_id2, match_set.walks2[idx2].front(),
                                                                                          local_scale,
                                                                                          chain_merge1, chain_merge2, switch_dists1, switch_dists2));
                                }
                                else  {
                                    initial_weight = 0.0;
                                }
                            }
                        }
                    }
                }
                if (sinks1) {
                    final_weight = std::numeric_limits<double>::lowest();
                    for (auto snk_id1 : *sinks1) {
                        for (auto snk_id2 : *sinks2) {
                            if ((chain_merge1.reachable(match_set.walks1[idx1].back(), snk_id1) ||
                                 match_set.walks1[idx1].back() == snk_id1) &&
                                (chain_merge2.reachable(match_set.walks2[idx2].back(), snk_id2) ||
                                 match_set.walks2[idx2].back() == snk_id2)) {
                                
                                if (score_edges) {
                                    final_weight = std::max(final_weight, edge_weight(match_set.walks1[idx1].back(), snk_id1,
                                                                                      match_set.walks2[idx2].back(), snk_id2,
                                                                                      local_scale,
                                                                                      chain_merge1, chain_merge2, switch_dists1, switch_dists2));
                                }
                                else  {
                                    final_weight = 0.0;
                                }
                            }
                            
                        }
                    }
                }
                anchor_graph.add_node(i, idx1, idx2, weight, initial_weight, final_weight);
            }
        }
    }
    
    // TODO: with no edge costs and positive weights, we can always do DP over the transitive
    // reduction, so we could speed this up by figuring out a better way to reduce transitive edges
    
    // add all possible edges
    for (uint64_t node_id1 = 0; node_id1 < anchor_graph.node_size(); ++node_id1) {
        for (uint64_t node_id2 = 0; node_id2 < anchor_graph.node_size(); ++node_id2) {
            size_t set1, idx11, idx21, set2, idx12, idx22;
            std::tie(set1, idx11, idx21) = anchor_graph.label(node_id1);
            std::tie(set2, idx12, idx22) = anchor_graph.label(node_id2);
            
            auto& match_set1 = match_sets[set1];
            auto& match_set2 = match_sets[set2];
            
            if (chain_merge1.reachable(match_set1.walks1[idx11].back(),
                                       match_set2.walks1[idx12].front()) &&
                chain_merge2.reachable(match_set1.walks2[idx21].back(),
                                       match_set2.walks2[idx22].front())) {
                if (score_edges) {
                    anchor_graph.add_edge(node_id1, node_id2,
                                          edge_weight(match_set1.walks1[idx11].back(), match_set2.walks1[idx12].front(),
                                                      match_set1.walks2[idx21].back(), match_set2.walks2[idx22].front(),
                                                      local_scale,
                                                      chain_merge1, chain_merge2, switch_dists1, switch_dists2));
                }
                else {
                    anchor_graph.add_edge(node_id1, node_id2);
                }
            }
        }
    }
    
    double min_score = 0.0;
    if (sources1 && sinks1 && score_edges) {
        // account for the score of the empty chain
        min_score = std::numeric_limits<double>::lowest();
        for (auto src_id1 : *sources1) {
            for (auto src_id2 : *sources2) {
                for (auto snk_id1 : *sinks1) {
                    for (auto snk_id2 : *sinks2) {
                        min_score = std::max<double>(min_score, edge_weight(src_id1, snk_id1, src_id2, snk_id2,
                                                                            local_scale,
                                                                            chain_merge1, chain_merge2, switch_dists1, switch_dists2));
                    }
                }
            }
        }
    }
    
    // get heaviest path and convert into a chain
    std::vector<anchor_t> chain;
    for (auto node_id : anchor_graph.heaviest_weight_path(min_score)) {
        size_t set, idx1, idx2;
        std::tie(set, idx1, idx2) = anchor_graph.label(node_id);
        chain.emplace_back();
        auto& match_set = match_sets[set];
        auto& chain_node = chain.back();
        chain_node.walk1 = match_set.walks1[idx1];
        chain_node.walk2 = match_set.walks2[idx2];
        chain_node.count1 = match_set.count1;
        chain_node.count2 = match_set.count2;
        chain_node.full_length = match_set.full_length;
        chain_node.match_set = set;
        chain_node.idx1 = idx1;
        chain_node.idx2 = idx2;
    }
    
    annotate_scores(chain);
    
    if (debug_anchorer) {
        std::cerr << "constructed anchor chain of size " << chain.size() << '\n';
    }
    
    return chain;
}

template<typename UIntDist, typename UIntSet, typename UIntMatch, typename UIntAnchor, typename ScoreFloat, class DistMatchVector, class AnchorVector, class MBank, class BGraph, class XMerge>
std::vector<anchor_t> Anchorer::sparse_chain_dp(const std::vector<match_set_t>& match_sets,
                                                const BGraph& graph1,
                                                const XMerge& chain_merge1,
                                                const XMerge& chain_merge2,
                                                size_t num_match_sets,
                                                bool suppress_verbose_logging,
                                                const std::vector<uint64_t>* sources1,
                                                const std::vector<uint64_t>* sources2,
                                                const std::vector<uint64_t>* sinks1,
                                                const std::vector<uint64_t>* sinks2,
                                                const std::unordered_set<std::tuple<size_t, size_t, size_t>>* masked_matches) const {
    
    assert((sources1 == nullptr) == (sources2 == nullptr));
    assert((sinks1 == nullptr) == (sinks2 == nullptr));
    
    if (!suppress_verbose_logging) {
        logging::log(logging::Debug, "Beginning free-gap sparse chaining algorithm");
    }
    
    MBank match_bank(graph1, match_sets, num_match_sets, suppress_verbose_logging, masked_matches);
    using match_id_t = typename MBank::match_id_t;
    
    // a key of (chain index, anchor set, walk1 index, walk2 index)
    using key_t = std::pair<UIntDist, match_id_t>;
    
    // for each chain2, the initial search tree data for chains, records of (key_t, weight)
    std::vector<std::vector<std::pair<key_t, ScoreFloat>>> search_tree_data(chain_merge2.chain_size());
    
    static const ScoreFloat mininf = std::numeric_limits<ScoreFloat>::lowest();
    
    if (debug_anchorer) {
        std::cerr << "gathering anchor endpoint information\n";
    }
    
    for (auto it = match_bank.begin(), end = match_bank.end(); it != end; ++it) {
        
        uint64_t chain2;
        size_t index2;
        std::tie(chain2, index2) = chain_merge2.chain(match_bank.walk2(*it).back());
        
        search_tree_data[chain2].emplace_back(key_t(index2, *it), mininf);
        
        auto& match_set = match_bank.match_set(*it);
        
        ScoreFloat weight = score_function->anchor_weight(match_set.count1, match_set.count2,
                                                          match_set.walks1.front().size(), match_set.full_length);
        
        if (sources1) {
            // not all nodes can be starts for DP, must be reachable by one of the sources
            bool found1 = false, found2 = false;
            for (auto src_id1 : *sources1) {
                if (src_id1 == match_bank.walk1(*it).front() || chain_merge1.reachable(src_id1, match_bank.walk1(*it).front())) {
                    found1 = true;
                    break;
                }
            }
            for (auto src_id2 : *sources2) {
                if (src_id2 == match_bank.walk2(*it).front() || chain_merge2.reachable(src_id2, match_bank.walk2(*it).front())) {
                    found2 = true;
                    break;
                }
            }
            if (!found1 || !found2) {
                weight = mininf;
            }
        }
        
        match_bank.update_dp(*it, weight, match_bank.max());
    }
    
    
    // sort one time to speed up construction across chains
    // TODO: we could use integer sort to exchange O(n log n) with O(V)
    for (auto& chain_search_tree_data : search_tree_data) {
        std::stable_sort(chain_search_tree_data.begin(), chain_search_tree_data.end());
    }
    
    // for each chain1, for each chain2, a tree over seed end chain indexes
    using MaxTree = MaxSearchTree<key_t, ScoreFloat, DistMatchVector, std::vector<ScoreFloat>, AnchorVector>;
    size_t tree_mem_size = 0;
    std::vector<std::vector<MaxTree>> search_trees;
    search_trees.resize(chain_merge1.chain_size());
    
    for (size_t i = 0; i < chain_merge1.chain_size(); ++i) {
        auto& chain_search_trees = search_trees[i];
        chain_search_trees.reserve(chain_merge2.chain_size());
        for (size_t j = 0; j < search_tree_data.size(); ++j) {
            chain_search_trees.emplace_back(search_tree_data[j]);
            if (logging::level >= logging::Debug) {
                tree_mem_size += chain_search_trees.back().memory_size();
            }
        }
    }
    // clear the search tree data
    {
        auto dummy = std::move(search_tree_data);
    }
    
    if (logging::level >= logging::Debug && !suppress_verbose_logging) {
        logging::log(logging::Debug, "Sparse query structures are occupying " + format_memory_usage(tree_mem_size) + " of memory.");
        logging::log(logging::Debug, "Current memory usage is " + format_memory_usage(current_memory_usage()) + ".");
    }
    
    if (debug_anchorer) {
        std::cerr << "computing forward edges\n";
    }
    
    // get the edges to chain neighbors for all of the nodes that have a match start
    std::vector<bool> mask_to, mask_from;
    std::tie(mask_to, mask_from) = generate_forward_edge_masks(graph1, match_sets, num_match_sets);
    auto forward_edges = chain_merge1.chain_forward_edges(&mask_to, &mask_from);
    {
        // clear out this memory
        auto dummy = std::move(mask_to);
        dummy = std::move(mask_from);
    }
    
    if (logging::level >= logging::Debug && !suppress_verbose_logging) {
        size_t forward_size = sizeof(forward_edges) + forward_edges.capacity() * sizeof(typename decltype(forward_edges)::value_type);
        for (const auto& row : forward_edges) {
            forward_size += row.capacity() * sizeof(typename decltype(forward_edges)::value_type::value_type);
        }
        
        logging::log(logging::Debug, "Forward edges are occupying " + format_memory_usage(forward_size) + " of memory.");
        logging::log(logging::Debug, "Current memory usage is " + format_memory_usage(current_memory_usage()) + ".");
    }
    
    if (debug_anchorer) {
        std::cerr << "beginning main DP iteration\n";
    }
    
    size_t iter = 0;
    for (uint64_t node_id : topological_order(graph1)) {
        
        ++iter;
        if (logging::level >= logging::Verbose && iter % 1000000 == 0 && !suppress_verbose_logging) {
            logging::log(logging::Verbose, "Iteration " + std::to_string(iter) + " of " + std::to_string(graph1.node_size()) + " in sparse chaining algorithm.");
        }
        
        if (debug_anchorer) {
            std::cerr << "on node " << node_id << '\n';
        }
        
        uint64_t chain1 = chain_merge1.chain(node_id).first;
        
        for (const auto& match_id : match_bank.ends_on(node_id)) {
            // we've hit the end of this match, so we can enter its DP value into the trees
            // to update other DP values as we move forward
            
            uint64_t chain2;
            size_t index2;
            std::tie(chain2, index2) = chain_merge2.chain(match_bank.walk2(match_id).back());
            
            auto& tree = search_trees[chain1][chain2];
            auto it = tree.find(key_t(index2, match_id));
            ScoreFloat dp_val = match_bank.dp_value(match_id);
            if ((*it).second < dp_val) {
                tree.update(it, dp_val);
            }
        }
        
        if (debug_anchorer) {
            std::cerr << "looking for chain forward edges\n";
        }
                
        // carry the current updates to any anchor starts for which this is the
        // the last node to reach from this chain (then all DP values from anchors
        // that end in this chain on graph1 have been completed)
        for (auto edge : forward_edges[node_id]) {
            
            uint64_t fwd_id, chain1;
            std::tie(fwd_id, chain1) = edge;
            
            if (debug_anchorer) {
                std::cerr << "there is a forward edge to graph1 node " << fwd_id << ", from node " << node_id << " on chain " << chain1 << '\n';
            }
            
            for (const auto& match_id : match_bank.starts_on(fwd_id)) {
                
                // an anchor starts here in graph1
                
                const auto& match_set = match_bank.match_set(match_id);
                
                // the weight of this anchors in this set
                ScoreFloat weight = score_function->anchor_weight(match_set.count1, match_set.count2,
                                                                  match_set.walks1.front().size(), match_set.full_length);
                
                // we will check for previous DP values that can reach this one in graph2 from
                // each of the chains in the chain partition of graph2
                for (uint64_t chain2 = 0; chain2 < chain_merge2.chain_size(); ++chain2) {
                    
                    auto chain_pred2 = chain_merge2.predecessor_index(match_bank.walk2(match_id).front(), chain2);
                    
                    if (chain_pred2 == -1) {
                        // there is no reachable node from this chain
                        continue;
                    }
                    if (debug_anchorer) {
                        std::cerr << "looking for predecessor on chain " << chain2 << " with predecessor at or before index " << chain_pred2 << '\n';
                    }
                    // find the max DP value up to (and including) the predecessor
                    const auto& tree = search_trees[chain1][chain2];
                    
                    auto it = tree.range_max(key_t(0, match_bank.min()),
                                             key_t(chain_pred2 + 1, match_bank.min())); // +1 because past-the-last
                    
                    if (it == tree.end()) {
                        // there aren't any ends of anchors before the predecessor in this chain
                        if (debug_anchorer) {
                            std::cerr << "there are no predecessors on this chain\n";
                        }
                        continue;
                    }
                    
                    // the weight of this anchor plus all previous anchors
                    ScoreFloat dp_weight = (*it).second + weight;
                    match_bank.update_dp(match_id, dp_weight, (*it).first.second);
                }
            }
        }
    }
    
    auto final_term = [&](const match_id_t& match_id) -> ScoreFloat {
        if (sinks1) {
            uint64_t last_id1 = match_bank.walk1(match_id).back();
            uint64_t last_id2 = match_bank.walk2(match_id).back();
            for (auto snk_id1 : *sinks1) {
                for (auto snk_id2 : *sinks2) {
                    if ((snk_id1 == last_id1 || chain_merge1.reachable(last_id1, snk_id1)) &&
                        (snk_id2 == last_id2 || chain_merge2.reachable(last_id2, snk_id2))) {
                        return 0.0;
                    }
                }
            }
            return mininf;
        }
        else {
            return 0.0;
        }   
    };
    
    return traceback_sparse_dp<ScoreFloat>(match_sets, match_bank, final_term, 0.0, suppress_verbose_logging);
}



template<typename UIntDist, class BGraph, class XMerge>
std::vector<std::vector<UIntDist>> Anchorer::post_switch_distances(const BGraph& graph, const XMerge& xmerge) const {
        
    std::vector<std::vector<UIntDist>> dists(graph.node_size(),
                                             std::vector<UIntDist>(xmerge.chain_size(), -1));
    
    // note: this DP is different from the paper because i have the predecessor on the same
    // path being the previous node rather than the node itself
    for (auto node_id : topological_order(graph)) {
        auto& row = dists[node_id];
        for (uint64_t p = 0; p < xmerge.chain_size(); ++p) {
            for (auto prev_id : graph.previous(node_id)) {
                auto pred = xmerge.predecessor_index(node_id, p);
                if (xmerge.index_on(prev_id, p) == pred) {
                    // switching paths you here immediately, no distance
                    row[p] = 0;
                    break;
                }
                else if (xmerge.predecessor_index(prev_id, p) == pred) {
                    // travel through this predecessor after switching onto it from path p
                    // note: this will overwrite the default -1 because of overflow
                    row[p] = std::min<UIntDist>(row[p], dists[prev_id][p] + graph.label_size(prev_id));
                }
            }
        }
    }
    
    return dists;
}

template<class BGraph>
std::pair<std::vector<bool>, std::vector<bool>> Anchorer::generate_forward_edge_masks(const BGraph& graph1, const std::vector<match_set_t>& match_sets, size_t num_match_sets) const {
    
    static const bool debug = false;
    
    std::pair<std::vector<bool>, std::vector<bool>> return_val;
    return_val.first.resize(graph1.node_size(), false);
    return_val.second.resize(graph1.node_size(), false);
    auto& have_match_start = return_val.first;
    auto& follow_match_end = return_val.second;
    
    // label all the nodes that are the start of some match
    for (size_t i = 0; i < num_match_sets; ++i) {
        for (const auto& walk1 : match_sets[i].walks1) {
            have_match_start[walk1.front()] = true;
        }
    }
    // label all the nodes that follow the end of some match
    for (size_t i = 0; i < num_match_sets; ++i) {
        for (const auto& walk1 : match_sets[i].walks1) {
            follow_match_end[walk1.back()] = true;
        }
    }
    if (debug) {
        std::cerr << "init match ends\n";
        for (size_t i = 0; i < follow_match_end.size(); ++i) {
            std::cerr << i << '\t' << follow_match_end[i] << '\n';
        }
    }
    std::vector<uint64_t> queue;
    for (uint64_t node_id1 = 0; node_id1 < graph1.node_size(); ++node_id1) {
        if (follow_match_end[node_id1]) {
            // DFS of unlabeled nodes
            queue.push_back(node_id1);
            while (!queue.empty()) {
                auto here = queue.back();
                queue.pop_back();
                for (auto next_id : graph1.next(here)) {
                    if (!follow_match_end[next_id]) {
                        if (debug) {
                            std::cerr << "walk to " << next_id << " from " << here << '\n';
                        }
                        follow_match_end[next_id] = true;
                        queue.push_back(next_id);
                    }
                }
            }
        }
    }
    
    if (debug) {
        std::cerr << "masks:\n";
        for (size_t i = 0; i < have_match_start.size(); ++i) {
            std::cerr << i << '\t' << have_match_start[i] << '\t' << follow_match_end[i] << '\n';
        }
    }
    
    return return_val;
}

template<typename UIntSet, typename UIntMatch, typename UIntDist, typename IntShift, typename UIntAnchor, typename ScoreFloat,
         class ShiftMatchVector, class DistMatchVector, class DistVector, class AnchorVector, class MBank, class BGraph, class XMerge, size_t NumPW>
std::vector<anchor_t> Anchorer::sparse_affine_chain_dp(const std::vector<match_set_t>& match_sets,
                                                       const BGraph& graph1,
                                                       const BGraph& graph2,
                                                       const XMerge& xmerge1,
                                                       const XMerge& xmerge2,
                                                       const std::array<double, NumPW>& gap_open,
                                                       const std::array<double, NumPW>& gap_extend,
                                                       double local_scale,
                                                       size_t num_match_sets,
                                                       bool suppress_verbose_logging,
                                                       const std::vector<uint64_t>* sources1,
                                                       const std::vector<uint64_t>* sources2,
                                                       const std::vector<uint64_t>* sinks1,
                                                       const std::vector<uint64_t>* sinks2,
                                                       const std::unordered_set<std::tuple<size_t, size_t, size_t>>* masked_matches) const {
    
    assert((sources1 == nullptr) == (sources2 == nullptr));
    assert((sinks1 == nullptr) == (sinks2 == nullptr));
    
    if (!suppress_verbose_logging) {
        logging::log(logging::Debug, "Beginning affine-gap sparse chaining algorithm with local scale " + std::to_string(local_scale));
    }
    
    if (debug_anchorer) {
        std::cerr << "chaining match sets:\n";
        for (UIntSet i = 0; i < num_match_sets; ++i) {
            std::cerr << "set " << i << ":\n";
            std::cerr << "\twalks on 1:\n";
            for (auto& walk : match_sets[i].walks1) {
                std::cerr << "\t\t";
                for (auto n : walk) {
                    std::cerr << n << ' ';
                }
                std::cerr << '\n';
            }
            std::cerr << "\twalks on 2:\n";
            for (auto& walk : match_sets[i].walks2) {
                std::cerr << "\t\t";
                for (auto n : walk) {
                    std::cerr << n << ' ';
                }
                std::cerr << '\n';
            }
        }
    }
    
    
    MBank match_bank(graph1, match_sets, num_match_sets, suppress_verbose_logging, masked_matches);
    using match_id_t = typename MBank::match_id_t;
    
    // (offset diff, match set, walk1, walk2)
    using key_t = std::pair<IntShift, match_id_t>;
    using gf_key_t = std::pair<UIntDist, match_id_t>;
    
    static const ScoreFloat mininf = std::numeric_limits<ScoreFloat>::lowest();
    
    // the D arrays from chandra & jain 2023
    auto switch_dists1 = post_switch_distances<UIntDist>(graph1, xmerge1);
    auto switch_dists2 = post_switch_distances<UIntDist>(graph2, xmerge2);
    
    // the gap contribution from a source anchor
    auto basic_source_shift = [&](uint64_t src_id1, uint64_t src_id2, uint64_t path1, uint64_t path2) -> IntShift {
        return xmerge1.index_on(src_id1, path1) - xmerge2.index_on(src_id2, path2);
    };
    auto source_shift = [&](const match_id_t& match_id, uint64_t path1, uint64_t path2) -> IntShift {
        return basic_source_shift(match_bank.walk1(match_id).back(), match_bank.walk2(match_id).back(), path1, path2);
    };
    // add the source coordinates to grab as backpointers
    auto get_key = [&](const match_id_t& match_id, uint64_t path1, uint64_t path2) -> key_t {
        return key_t(source_shift(match_id, path1, path2), match_id);
    };
    // the gap contribution from the destination anchor
    auto basic_query_shift = [&](uint64_t query_id1, uint64_t query_id2, uint64_t path1, uint64_t path2) -> IntShift {
        return (xmerge1.predecessor_index(query_id1, path1) - xmerge2.predecessor_index(query_id2, path2)
                + switch_dists1[query_id1][path1] - switch_dists2[query_id2][path2]);
    };
    auto query_shift = [&](const match_id_t& match_id, uint64_t path1, uint64_t path2) -> IntShift {
        return basic_query_shift(match_bank.walk1(match_id).front(), match_bank.walk2(match_id).front(), path1, path2);
    };
    // the offset on path from graph2
    auto get_key_offset = [&](const match_id_t& match_id, uint64_t path2) -> UIntDist {
        return xmerge2.index_on(match_bank.walk2(match_id).back(), path2);
    };
    // the "effective offset" on path of graph2 (true offsets on the chain before this can reach it)
    auto get_query_offset = [&](const match_id_t& match_id, uint64_t path2) -> UIntDist {
        // note: we rely on -1's overflowing to 0
        return xmerge2.predecessor_index(match_bank.walk2(match_id).front(), path2) + 1;
    };
    auto get_gap_free_key = [&](const match_id_t& match_id, uint64_t path2) {
        return gf_key_t(get_key_offset(match_id, path2), match_id);
    };
    // directly measure a gap between two node pairs
    auto score_gap = [&](IntShift gap) -> ScoreFloat {
        // score the gap
        ScoreFloat score = mininf;
        if (gap == 0) {
            score = 0.0;
        }
        else if (gap != std::numeric_limits<IntShift>::max()) {
            for (int pw = 0; pw < NumPW; ++pw) {
                score = std::max<ScoreFloat>(score, -local_scale * (gap_open[pw] + gap_extend[pw] * std::abs(gap)));
            }
        }
        return score;
    };
    auto measure_gap = [&](uint64_t prev_id1, uint64_t prev_id2, uint64_t curr_id1, uint64_t curr_id2) -> IntShift {
        
        IntShift gap = std::numeric_limits<IntShift>::max();
        
        if ((prev_id1 == curr_id1 || xmerge1.reachable(prev_id1, curr_id1)) &&
            (prev_id2 == curr_id2 || xmerge2.reachable(prev_id2, curr_id2))) {
            // iterate over combos of source ids and their paths
            for (auto p1 : xmerge1.chains_on(prev_id1)) {
                for (auto p2 : xmerge2.chains_on(prev_id2)) {
                    IntShift gap_here = basic_source_shift(prev_id1, prev_id2, p1, p2) - basic_query_shift(curr_id1, curr_id2, p1, p2);
                    if (std::abs(gap_here) < std::abs(gap)) {
                        gap = gap_here;
                    }
                }
            }
        }
        return gap;
    };
    auto measure_gap_nn = [&](uint64_t prev_id1, uint64_t prev_id2, uint64_t curr_id1, uint64_t curr_id2) -> std::pair<IntShift, ScoreFloat> {
        
        std::pair<IntShift, ScoreFloat> return_val;
        return_val.first = measure_gap(prev_id1, prev_id2, curr_id1, curr_id2);
        return_val.second = score_gap(return_val.first);
        return return_val;
    };

    // measure from set to node pair
    auto measure_gap_sn = [&](const std::vector<uint64_t>& prev1,
                              const std::vector<uint64_t>& prev2,
                              uint64_t curr_id1, uint64_t curr_id2) -> std::pair<IntShift, ScoreFloat> {
        
        std::pair<IntShift, ScoreFloat> return_val(std::numeric_limits<IntShift>::max(), mininf);
        for (uint64_t prev_id1 : prev1) {
            for (uint64_t prev_id2 : prev2) {
                IntShift gap_here = measure_gap(prev_id1, prev_id2, curr_id1, curr_id2);
                if (std::abs(gap_here) < return_val.first) {
                    return_val.first = gap_here;
                }
            }
        }
        return_val.second = score_gap(return_val.first);
        return return_val;
    };
    // measure from node pair to set
    auto measure_gap_ns = [&](uint64_t prev_id1, uint64_t prev_id2,
                              const std::vector<uint64_t>& curr1,
                              const std::vector<uint64_t>& curr2) -> std::pair<IntShift, ScoreFloat> {
        
        std::pair<IntShift, ScoreFloat> return_val(std::numeric_limits<IntShift>::max(), mininf);
        for (uint64_t curr_id1 : curr1) {
            for (uint64_t curr_id2 : curr2) {
                IntShift gap_here = measure_gap(prev_id1, prev_id2, curr_id1, curr_id2);
                if (std::abs(gap_here) < return_val.first) {
                    return_val.first = gap_here;
                }
            }
        }
        return_val.second = score_gap(return_val.first);
        return return_val;
    };
    // measure from set to set
    auto measure_gap_ss = [&](const std::vector<uint64_t>& prev1,
                              const std::vector<uint64_t>& prev2,
                              const std::vector<uint64_t>& curr1,
                              const std::vector<uint64_t>& curr2) -> std::pair<IntShift, ScoreFloat> {
        
        std::pair<IntShift, ScoreFloat> return_val(std::numeric_limits<IntShift>::max(), mininf);
        for (uint64_t curr_id1 : curr1) {
            for (uint64_t curr_id2 : curr2) {
                for (uint64_t prev_id1 : prev1) {
                    for (uint64_t prev_id2 : prev2) {
                        IntShift gap_here = measure_gap(prev_id1, prev_id2, curr_id1, curr_id2);
                        if (std::abs(gap_here) < return_val.first) {
                            return_val.first = gap_here;
                        }
                    }
                }
            }
        }
        return_val.second = score_gap(return_val.first);
        return return_val;
    };
        
    // for each path1, for each path2, list of (search index, value)
    std::vector<std::vector<std::vector<std::tuple<key_t, UIntDist, ScoreFloat>>>> search_tree_data;
    search_tree_data.resize(xmerge1.chain_size(), std::vector<std::vector<std::tuple<key_t, UIntDist, ScoreFloat>>>(xmerge2.chain_size()));
    
    // for each path1, for each path2, for each shift value, list of (offset, value)
    std::vector<std::vector<std::deque<std::vector<std::pair<gf_key_t, ScoreFloat>>>>> gap_free_search_tree_data;
    gap_free_search_tree_data.resize(xmerge1.chain_size(), std::vector<std::deque<std::vector<std::pair<gf_key_t, ScoreFloat>>>>(xmerge2.chain_size()));
    // for each path1, for each path2, shift value represented by the start of the corresponding deque
    std::vector<std::vector<IntShift>> min_shift(xmerge1.chain_size(), std::vector<IntShift>(xmerge2.chain_size()));
    
    if (debug_anchorer) {
        std::cerr << "doing bookkeeping for sparse affine algorithm\n";
    }
    
    uint64_t num_pairs = 0;
    
    // do the bookkeeping
    for (auto it = match_bank.begin(), end = match_bank.end(); it != end; ++it) {
        ++num_pairs;
        
        if (debug_anchorer) {
            std::cerr << "initializing data for match " << std::get<0>(match_bank.get_match_indexes(*it)) << ", " << std::get<1>(match_bank.get_match_indexes(*it)) << ", " << std::get<2>(match_bank.get_match_indexes(*it)) << '\n';
        }
        
        const auto& match_set = match_bank.match_set(*it);
        
        // initialize the DP structure with a single-anchor chain at each position
        ScoreFloat weight = score_function->anchor_weight(match_set.count1, match_set.count2,
                                                      match_set.walks1.front().size(), match_set.full_length);
        
        if (sources1) {
            
            // TODO: could i get rid of these loops by querying the structure somehow?
            
            ScoreFloat lead_indel_score = measure_gap_sn(*sources1, *sources2,
                                                         match_bank.walk1(*it).front(), match_bank.walk2(*it).front()).second;
            if (lead_indel_score == mininf) {
                // this anchor was not reachable
                weight = mininf;
            }
            else {
                weight += lead_indel_score;
            }
        }
        
        match_bank.update_dp(*it, weight, match_bank.max());
        
        for (auto p1 : xmerge1.chains_on(match_bank.walk1(*it).back())) {
            for (auto p2 : xmerge2.chains_on(match_bank.walk2(*it).back())) {
                // initialize the sparse value data
                search_tree_data[p1][p2].emplace_back(get_key(*it, p1, p2), get_key_offset(*it, p2), mininf);
                
                
                // initialize the gap free value in its appropriate bank
                auto shift = source_shift(*it, p1, p2);
                if (debug_anchorer) {
                    std::cerr << "adding match to gap free data for shift " << shift << " on chain combo " << p1 << ", " << p2 << '\n';
                }
                auto& gf_data = gap_free_search_tree_data[p1][p2];
                auto& gf_data_min = min_shift[p1][p2];
                if (gf_data.empty()) {
                    gf_data_min = shift;
                    gf_data.emplace_front();
                    gf_data.front().emplace_back(get_gap_free_key(*it, p2), mininf);
                }
                else {
                    while (gf_data_min > shift) {
                        gf_data.emplace_front();
                        --gf_data_min;
                    }
                    while (gf_data_min + IntShift(gf_data.size()) <= shift) {
                        gf_data.emplace_back();
                    }
                    
                    gf_data[shift - gf_data_min].emplace_back(get_gap_free_key(*it, p2), mininf);
                }
            }
        }
    }
    
    if (logging::level >= logging::Debug && !suppress_verbose_logging) {
        // measure fine grain memory usage from local structs
        
        // search tree data
        size_t search_data_size = sizeof(search_tree_data) + search_tree_data.capacity() * sizeof(typename decltype(search_tree_data)::value_type);
        for (const auto& row : search_tree_data) {
            search_data_size += row.capacity() * sizeof(typename decltype(search_tree_data)::value_type::value_type);
            for (const auto& match_list : row) {
                search_data_size += match_list.capacity() * sizeof(typename decltype(search_tree_data)::value_type::value_type::value_type);
            }
        }
        // FIXME: include the gap free search tree data
        
        
        size_t switch_dists_size = sizeof(switch_dists1) + sizeof(switch_dists2) + sizeof(typename decltype(switch_dists1)::value_type) * (switch_dists1.capacity() + switch_dists2.capacity());
        for (const auto& row : switch_dists1) {
            switch_dists_size += row.capacity() * sizeof(typename decltype(switch_dists1)::value_type::value_type);
        }
        for (const auto& row : switch_dists2) {
            switch_dists_size += row.capacity() * sizeof(typename decltype(switch_dists2)::value_type::value_type);
        }
        
        logging::log(logging::Debug, "Initialized search tree data is occupying " + format_memory_usage(search_data_size) + ".");
        logging::log(logging::Debug, "Post switch distances are occupying " + format_memory_usage(switch_dists_size) + ".");
    }
    
    if (!suppress_verbose_logging) {
        logging::log(logging::Verbose, "Chaining " + std::to_string(num_pairs) + " matches.");
    }
    
    if (debug_anchorer) {
        std::cerr << "constructing search trees\n";
    }
    
    if (!suppress_verbose_logging) {
        logging::log(logging::Debug, "Initializing sparse query data structures");
    }
    static const bool verbose_memory = true;
    
    // grids of search trees over (path1, path2) where scores correspond to different distance scenaries
    // odds        d1 > d2
    // evens       d1 < d2
    using OrthoMaxTree = OrthogonalMaxSearchTree<key_t, UIntDist, ScoreFloat, UIntAnchor, ShiftMatchVector, DistVector, std::vector<ScoreFloat>, AnchorVector>;
    std::array<std::vector<std::vector<OrthoMaxTree>>, 2 * NumPW> search_trees;
    for (size_t pw = 0; pw < 2 * NumPW; ++pw) {
        auto& pw_trees = search_trees[pw];
        pw_trees.resize(xmerge1.chain_size());
        for (uint64_t p1 = 0; p1 < xmerge1.chain_size(); ++p1) {
            auto& tree_row = pw_trees[p1];
            tree_row.reserve(xmerge2.chain_size());
            for (uint64_t p2 = 0; p2 < xmerge2.chain_size(); ++p2) {
                tree_row.emplace_back(search_tree_data[p1][p2]);
                if (verbose_memory && !suppress_verbose_logging) {
                    logging::log(logging::Debug, "Search tree for path combination (" + std::to_string(p1) + ", " + std::to_string(p2) + ") in piece-wise component " + std::to_string(pw) + " is occupying " + format_memory_usage(tree_row.back().memory_size()) + ".");
                    logging::log(logging::Debug, "Current memory usage is " + format_memory_usage(current_memory_usage()) + ".");
                }
            }
        }
    }
    if (logging::level >= logging::Debug && !suppress_verbose_logging) {
        size_t gapped_tree_mem_size = 0;
        for (size_t pw = 0; pw < 2 * NumPW; ++pw) {
            const auto& pw_trees = search_trees[pw];
            for (uint64_t p1 = 0; p1 < xmerge1.chain_size(); ++p1) {
                const auto& tree_row = pw_trees[p1];
                for (uint64_t p2 = 0; p2 < xmerge2.chain_size(); ++p2) {
                    gapped_tree_mem_size += tree_row[p2].memory_size();
                }
            }
        }
        
        logging::log(logging::Debug, "Gapped sparse query structures are occupying " + format_memory_usage(gapped_tree_mem_size) + " of memory.");
        logging::log(logging::Debug, "Current memory usage is " + format_memory_usage(current_memory_usage()) + ".");
    }
    
    // for each path1, for each path2, for each shift value, a search tree
    using MaxTree = MaxSearchTree<gf_key_t, ScoreFloat, DistMatchVector, std::vector<ScoreFloat>, AnchorVector>;
    std::vector<std::vector<std::vector<MaxTree>>> gap_free_search_trees;
    gap_free_search_trees.resize(xmerge1.chain_size(), std::vector<std::vector<MaxTree>>(xmerge2.chain_size()));
    
    for (uint64_t p1 = 0; p1 < xmerge1.chain_size(); ++p1) {
        for (uint64_t p2 = 0; p2 < xmerge2.chain_size(); ++p2) {
            auto& tree_bank = gap_free_search_trees[p1][p2];
            auto& tree_data_vecs = gap_free_search_tree_data[p1][p2];
            tree_bank.reserve(tree_data_vecs.size());
            for (size_t i = 0; i < tree_data_vecs.size(); ++i) {
                tree_bank.emplace_back(tree_data_vecs[i]);
            }
            if (verbose_memory && !suppress_verbose_logging) {
                size_t mem_size = 0;
                size_t num_nonempty = 0;
                for (const auto& tree : tree_bank) {
                    mem_size += tree.memory_size();
                    num_nonempty += (tree.empty() ? 0 : 1);
                }
                logging::log(logging::Debug, "Ungapped sparse query structures for path combination (" + std::to_string(p1) + ", " + std::to_string(p2) + ") are occupying " + format_memory_usage(mem_size) + " of memory with " + std::to_string(num_nonempty) + " nonempty trees out of " + std::to_string(tree_bank.size()) + " total.");
                logging::log(logging::Debug, "Current memory usage is " + format_memory_usage(current_memory_usage()) + ".");
            }
        }
    }
    
    
    if (logging::level >= logging::Debug && !suppress_verbose_logging) {
        
        size_t ungapped_tree_mem_size = 0;
        for (uint64_t p1 = 0; p1 < xmerge1.chain_size(); ++p1) {
            for (uint64_t p2 = 0; p2 < xmerge2.chain_size(); ++p2) {
                for (const auto& tree : gap_free_search_trees[p1][p2]) {
                    ungapped_tree_mem_size += tree.memory_size();
                }
            }
        }
        logging::log(logging::Debug, "Ungapped sparse query structures are occupying " + format_memory_usage(ungapped_tree_mem_size) + " of memory.");
        logging::log(logging::Debug, "Current memory usage is " + format_memory_usage(current_memory_usage()) + ".");
    }
    
    // clear the search tree data
    {
        decltype(search_tree_data) dummy = std::move(search_tree_data);
        decltype(gap_free_search_tree_data) dummy2 = std::move(gap_free_search_tree_data);
    }
    
    if (debug_anchorer) {
        std::cerr << "computing forward edges\n";
    }
    
    // get the edges to chain neighbors
    std::vector<bool> mask_to, mask_from;
    std::tie(mask_to, mask_from) = generate_forward_edge_masks(graph1, match_sets, num_match_sets);
    auto forward_edges = xmerge1.chain_forward_edges(&mask_to, &mask_from);
    {
        // clear out this memory
        auto dummy = std::move(mask_to);
        dummy = std::move(mask_from);
    }
    
    if (logging::level >= logging::Debug && !suppress_verbose_logging) {
        size_t forward_size = sizeof(forward_edges) + forward_edges.capacity() * sizeof(typename decltype(forward_edges)::value_type);
        for (const auto& row : forward_edges) {
            forward_size += row.capacity() * sizeof(typename decltype(forward_edges)::value_type::value_type);
        }
        
        logging::log(logging::Debug, "Forward edges are occupying " + format_memory_usage(forward_size) + " of memory.");
        logging::log(logging::Debug, "Current memory usage is " + format_memory_usage(current_memory_usage()) + ".");
    }
    
    if (debug_anchorer) {
        std::cerr << "beginning main DP iteration\n";
    }
    
    if (!suppress_verbose_logging) {
        logging::log(logging::Debug, "Beginning sparse dynamic programming");
    }
    
    size_t iter = 0;
    for (uint64_t node_id : topological_order(graph1)) {
        
        ++iter;
        if (logging::level >= logging::Verbose && iter % 250000 == 0 && !suppress_verbose_logging) {
            logging::log(logging::Verbose, "Iteration " + std::to_string(iter) + " of " + std::to_string(graph1.node_size()) + " in sparse chaining algorithm.");
        }
        
        if (debug_anchorer) {
            std::cerr << "on node " << node_id << '\n';
        }
        
        for (const auto& match_id : match_bank.ends_on(node_id)) {
            
            if (debug_anchorer) {
                std::cerr << "match " << std::get<0>(match_bank.get_match_indexes(match_id)) << ", " << std::get<1>(match_bank.get_match_indexes(match_id)) << ", " << std::get<2>(match_bank.get_match_indexes(match_id)) << " is on this node\n";
            }
            
            ScoreFloat dp_val = match_bank.dp_value(match_id);
            
            for (auto p1 : xmerge1.chains_on(match_bank.walk1(match_id).back())) {
                for (auto p2 : xmerge2.chains_on(match_bank.walk2(match_id).back())) {
                    
                    auto key1 = get_key(match_id, p1, p2);
                    auto key2 = get_key_offset(match_id, p2);
                    
                    if (debug_anchorer) {
                        std::cerr << "extending for path combo " << p1 << ',' << p2 << ", giving shift key " << key1.first << " and offset key " << key2 << '\n';
                    }
                    auto shift = source_shift(match_id, p1, p2);
                    {
                        // update the within diagonal search tree
                        auto& tree = gap_free_search_trees[p1][p2][shift - min_shift[p1][p2]];
                        auto it = tree.find(get_gap_free_key(match_id, p2));
                        tree.update(it, dp_val);
                    }
                    for (size_t pw = 0; pw < 2 * NumPW; ++pw) {
                        // save the anchor-independent portion of the score in the search tree
                        ScoreFloat value;
                        if (pw % 2 == 1) {
                            // d1 > d2
                            value = dp_val + local_scale * gap_extend[pw / 2] * shift;
                        }
                        else {
                            // d1 < d2
                            value = dp_val - local_scale * gap_extend[pw / 2] * shift;
                        }
                        auto& tree = search_trees[pw][p1][p2];
                        auto it = tree.find(key1, key2);
                        if (value > std::get<2>(*it)) {
                            // TODO: shouldn't this condition always be met, since keys are unique to a match pair?
                            tree.update(it, value);
                        }
                    }
                }
            }
        }
        
        
        if (debug_anchorer) {
            std::cerr << "looking for forward edges\n";
        }
        
        for (auto edge : forward_edges[node_id]) {
            
            uint64_t fwd_id, chain1;
            std::tie(fwd_id, chain1) = edge;
            
            if (debug_anchorer) {
                std::cerr << "there is a forward edge to " << fwd_id << ", from node " << node_id << " on chain " << chain1 << '\n';
            }
            
            for (const auto& match_id : match_bank.starts_on(fwd_id)) {
                // an anchor starts here in graph1
                                
                const auto& match_set = match_bank.match_set(match_id);
                
                // the weight of this anchors in this set
                ScoreFloat weight = score_function->anchor_weight(match_set.count1, match_set.count2,
                                                                  match_set.walks1.front().size(), match_set.full_length);
                
                for (uint64_t chain2 = 0; chain2 < xmerge2.chain_size(); ++chain2) {
                    // note: we have to check all of the chains because the best distance measure might
                    // not originate from a path that contains the head of the path
                    
                    IntShift query = query_shift(match_id, chain1, chain2);
                    UIntDist offset = get_query_offset(match_id, chain2);
                    if (debug_anchorer) {
                        std::cerr << "query shift is " << query << " and offset is " << offset << " on chain combo " << chain1 << "," << chain2 << '\n';
                    }
                    if (query >= min_shift[chain1][chain2] && query - min_shift[chain1][chain2] < gap_free_search_trees[chain1][chain2].size()) {
                        // check within the same diagonal
                        const auto& tree = gap_free_search_trees[chain1][chain2][query - min_shift[chain1][chain2]];
                        auto it = tree.range_max(gf_key_t(0, match_bank.min()), gf_key_t(offset, match_bank.min()));
                        if (it != tree.end()) {
                            ScoreFloat value = (*it).second + weight;
                            match_bank.update_dp(match_id, value, (*it).first.second);
                        }
                    }
                    for (size_t pw = 0; pw < 2 * NumPW; ++pw) {
                        // combine the anchor-dependent and anchor-independent portions of the score
                        const auto& tree = search_trees[pw][chain1][chain2];

                        if (pw % 2 == 1) {
                            // d1 > d2, search leftward of the query value
                            auto it = tree.range_max(key_t(std::numeric_limits<IntShift>::min(), match_bank.min()),
                                                     key_t(query, match_bank.min()),
                                                     0, offset);
                            if (it != tree.end()) {
                                ScoreFloat value = std::get<2>(*it) + weight - local_scale * (gap_open[pw / 2] + gap_extend[pw / 2] * query);
                                match_bank.update_dp(match_id, value, std::get<0>(*it).second);
                            }
                        }
                        else {
                            auto it = tree.range_max(key_t(query + 1, match_bank.min()),
                                                     key_t(std::numeric_limits<IntShift>::max(), match_bank.max()),
                                                     0, offset);
                            if (it != tree.end()) {
                                ScoreFloat value = std::get<2>(*it) + weight - local_scale * (gap_open[pw / 2] - gap_extend[pw / 2] * query);
                                match_bank.update_dp(match_id, value, std::get<0>(*it).second);
                            }
                        }
                    }
                }
            }
        }
    }
    
    ScoreFloat min_score = 0.0;
    if (sources1 && sinks1) {
        // we have to do better than the empty chain, so let's figure out its score
        
        min_score = measure_gap_ss(*sources1, *sources2, *sinks1, *sinks2).second;
    }
    
    if (debug_anchorer) {
        std::cerr << "traceback must exceed " << min_score << " to be counted\n";
    }
    
    // the score for the final indel
    auto final_indel_score = [&](const match_id_t& match_id) -> ScoreFloat {
        
        if (!sinks1) {
            return 0.0;
        }
        
        return measure_gap_ns(match_bank.walk1(match_id).back(), match_bank.walk2(match_id).back(), *sinks1, *sinks2).second;
    };
        
    auto traceback = traceback_sparse_dp<ScoreFloat>(match_sets, match_bank, final_indel_score, min_score, suppress_verbose_logging);
    
    // annotate the gap length and score between the anchors
    for (size_t i = 0; i < traceback.size(); ++i) {
        auto& anchor = traceback[i];
        if (i == 0) {
            if (sources1) {
                auto gap = measure_gap_sn(*sources1, *sources2, anchor.walk1.front(), anchor.walk2.front());
                anchor.gap_before = gap.first;
                anchor.gap_score_before = gap.second;
            }
        }
        else {
            auto& prev_anchor = traceback[i - 1];
            auto gap = measure_gap_nn(prev_anchor.walk1.back(), prev_anchor.walk2.back(),
                                      anchor.walk1.front(), anchor.walk2.front());
            prev_anchor.gap_after = gap.first;
            prev_anchor.gap_score_after = gap.second;
            anchor.gap_before = gap.first;
            anchor.gap_score_before = gap.second;
        }
        if (i + 1 == traceback.size()) {
            if (sinks1) {
                auto gap = measure_gap_ns(anchor.walk1.back(), anchor.walk2.back(), *sinks1, *sinks2);
                anchor.gap_after = gap.first;
                anchor.gap_score_after = gap.second;
            }
        }
    }
    
    return traceback;
}

template<typename ScoreFloat, class FinalFunc, class MBank>
std::vector<anchor_t> Anchorer::traceback_sparse_dp(const std::vector<match_set_t>& match_sets,
                                                    const MBank& match_bank,
                                                    const FinalFunc& final_function, ScoreFloat min_score,
                                                    bool suppress_verbose_logging) const {
    
    if (debug_anchorer) {
        std::cerr << "finding optimum\n";
    }
    
    // find the optimum dynamic programming values
    ScoreFloat opt_value = std::numeric_limits<ScoreFloat>::lowest();
    auto opt_match = match_bank.max();
    for (auto it = match_bank.begin(), end = match_bank.end(); it != end; ++it) {
        ScoreFloat dp_val = match_bank.dp_value(*it);
        ScoreFloat final_term = final_function(*it);
        if (final_term == std::numeric_limits<ScoreFloat>::lowest()) {
            dp_val = final_term;
        }
        else {
            dp_val += final_term;
        }
        if (dp_val > opt_value && dp_val > min_score) {
            opt_value = dp_val;
            opt_match = *it;
        }
    }
    
    if (debug_anchorer) {
        std::cerr << "doing traceback from opt score " << opt_value << "\n";
    }
    
    // traceback into a chain
    std::vector<anchor_t> anchors;
    auto here = opt_match;
    while (here != match_bank.max()) {
        
        auto indexes = match_bank.get_match_indexes(here);
        if (debug_anchorer) {
            std::cerr << "following traceback to set " << std::get<0>(indexes) << ", walk pair " << std::get<1>(indexes) << " " << std::get<2>(indexes) << '\n';
        }
        
        // grab the anchors that we used from their set
        anchors.emplace_back();
        auto& anchor = anchors.back();
        
        auto& match_set = match_sets[std::get<0>(indexes)];
        anchor.walk1 = match_set.walks1[std::get<1>(indexes)];
        anchor.count1 = match_set.count1;
        anchor.walk2 = match_set.walks2[std::get<2>(indexes)];
        anchor.count2 = match_set.count2;
        anchor.full_length = match_set.full_length;
        anchor.match_set = std::get<0>(indexes);
        anchor.idx1 = std::get<1>(indexes);
        anchor.idx2 = std::get<2>(indexes);
        
        // follow the backpointer from the DP structure
        here = match_bank.backpointer(here);
    }
    
    // take out of reverse order
    std::reverse(anchors.begin(), anchors.end());
    
    annotate_scores(anchors);
    
    if (debug_anchorer) {
        std::cerr << "completed sparse chaining\n";
    }
    
    if (!suppress_verbose_logging) {
        logging::log(logging::Debug, "Optimal chain consists of " + std::to_string(anchors.size()) + " matches with score " + (opt_value == std::numeric_limits<ScoreFloat>::lowest() ? std::string("-inf") : std::to_string(opt_value)));
    }
    
    return anchors;
}


template<class XMerge>
double Anchorer::edge_weight(uint64_t from_id1, uint64_t to_id1, uint64_t from_id2, uint64_t to_id2,
                             double local_scale,
                             const XMerge& xmerge1, const XMerge& xmerge2,
                             const std::vector<std::vector<size_t>>& switch_dists1,
                             const std::vector<std::vector<size_t>>& switch_dists2) const {
    
    
    double weight = std::numeric_limits<double>::lowest();
    for (auto chain1 : xmerge1.chains_on(from_id1)) {
        for (auto chain2 : xmerge2.chains_on(from_id2)) {
            int64_t dist1 = (xmerge1.predecessor_index(to_id1, chain1)
                             - xmerge1.index_on(from_id1, chain1) + switch_dists1[to_id1][chain1]);
            int64_t dist2 = (xmerge2.predecessor_index(to_id2, chain2)
                             - xmerge2.index_on(from_id2, chain2) + switch_dists2[to_id2][chain2]);
            
            int64_t gap = abs(dist1 - dist2);
            
            if (gap == 0) {
                // this is the best possible score
                weight = 0.0;
            }
            else {
                for (size_t i = 0; i < gap_open.size(); ++i) {
                    double comp_weight = -local_scale * (gap_open[i] + gap_extend[i] * gap);
                    weight = std::max(weight, comp_weight);
                }
            }
        }
    }
    return weight;
}

template<class BGraph, class XMerge>
void Anchorer::instrument_anchor_chain(const std::vector<anchor_t>& chain, double local_scale,
                                       const BGraph& graph1, const BGraph& graph2,
                                       const XMerge& xmerge1, const XMerge& xmerge2) const {
    
    std::vector<std::vector<size_t>> switch_dists1, switch_dists2;
    if (chaining_algorithm == SparseAffine) {
        switch_dists1 = std::move(post_switch_distances<size_t>(graph1, xmerge1));
        switch_dists2 = std::move(post_switch_distances<size_t>(graph2, xmerge2));
    }
    
    for (size_t i = 0; i < chain.size(); ++i) {
        const auto& anchor = chain[i];
        std::cerr << '@' << '\t' << i << '\t' << anchor.walk1.size() << '\t' << anchor.count1 << '\t' << anchor.count2 << '\t' << (anchor.count1 * anchor.count2) << '\t' << anchor.full_length << '\t' << score_function->anchor_weight(anchor.count1, anchor.count2, anchor.walk1.size(), anchor.full_length) << '\t';
        if (chaining_algorithm != SparseAffine || i == 0) {
            std::cerr << 0.0;
        }
        else {
            std::cerr << edge_weight(chain[i-1].walk1.back(), anchor.walk1.front(),
                                     chain[i-1].walk2.back(), anchor.walk2.front(),
                                     local_scale,
                                     xmerge1, xmerge2, switch_dists1, switch_dists2);
        }
        std::cerr << '\t' << anchor.walk1.front() << '\t' << anchor.walk1.back() << '\t' << anchor.walk2.front() << '\t' << anchor.walk2.back();
        size_t size1, size2;
        if (i == 0) {
            size1 = anchor.walk1.front();
            size2 = anchor.walk2.front();
            std::cerr << '\t' << anchor.walk1.front() << '\t' << anchor.walk2.front();
        }
        else {
            size1 = (anchor.walk1.front() - chain[i - 1].walk1.back());
            size2 = (anchor.walk2.front() - chain[i - 1].walk2.back());
        }
        std::vector<double> sizes{(double) size1, (double) size2};
        std::cerr << '\t' << size1 << '\t' << size2 << '\t' << generalized_mean(sizes.begin(), sizes.end(), -0.5);
        std::cerr << '\t' << anchor.gap_before << '\t'<< anchor.gap_score_before << '\t' << anchor.gap_after << '\t' << anchor.gap_score_after << '\n';
    }
}


}

#endif /* centrolign_anchorer_hpp */
