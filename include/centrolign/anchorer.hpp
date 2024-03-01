#ifndef centrolign_anchorer_hpp
#define centrolign_anchorer_hpp

#include <vector>
#include <cstdint>
#include <cmath>
#include <iostream>
#include <limits>
#include <array>

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

namespace centrolign {

struct anchor_t;

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
    
    static std::pair<int64_t, int64_t> source_sink_minmax(const SubGraphInfo& extraction);
    
    static std::vector<size_t> get_logging_indexes(size_t size);
    
private:
    
    template<class BGraph1, class BGraph2, class XMerge1, class XMerge2>
    static std::pair<SubGraphInfo, SubGraphInfo>
    do_extraction(uint64_t from1, uint64_t to1, uint64_t from2, uint64_t to2,
                  const BGraph1& graph1, const BGraph2& graph2,
                  const XMerge1& chain_merge1, const XMerge2& chain_merge2);
    
    template<class BGraph>
    static void do_project(const BGraph& graph, SubGraphInfo& subgraph,
                           const StepIndex& step_index);
};


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
    double score = 0.0;
    size_t match_set = -1;
    size_t idx1 = -1;
    size_t idx2 = -1;
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
                                       std::unordered_set<std::tuple<size_t, size_t, size_t>>* masked_matches = nullptr) const;
    
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
    using dp_entry_t = std::tuple<double, size_t, size_t, size_t>;
    
    // assumes that the graphs have already been given unique sentinels
    // note: these will never show up in anchors because they can't match
    
    template<class BGraph, class XMerge>
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
    
    template<class BGraph, class XMerge>
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
    
    template<size_t NumPW, class BGraph, class XMerge>
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
    template<class BGraph, class XMerge>
    std::vector<std::vector<size_t>> post_switch_distances(const BGraph& graph, const XMerge& xmerge) const;
    
    template<class FinalFunc>
    std::vector<anchor_t> traceback_sparse_dp(const std::vector<match_set_t>& match_sets,
                                              const std::vector<std::vector<std::vector<dp_entry_t>>>& dp,
                                              const FinalFunc& final_function, double min_score,
                                              bool suppress_verbose_logging,
                                              const std::unordered_set<std::tuple<size_t, size_t, size_t>>* masked_matches) const;
    
    
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
                                const XMerge& xmerge1, const XMerge& xmerge2) const;
    
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
    
    size_t next_log_idx = 0;
    std::vector<size_t> logging_indexes;
    if (logging::level >= logging::Debug) {
        logging_indexes = get_logging_indexes(anchor_chain.size());
        
        logging::log(logging::Debug, "Extracting graphs in between chain of " + std::to_string(anchor_chain.size()) + " anchors");
    }
    
    std::vector<std::pair<SubGraphInfo, SubGraphInfo>> stitch_pairs;
    stitch_pairs.reserve(anchor_chain.size() + 1);
    
    if (anchor_chain.empty()) {
        stitch_pairs.emplace_back(do_extraction(tableau1.src_id, tableau1.snk_id, tableau2.src_id, tableau2.snk_id,
                                                graph1, graph2, xmerge1, xmerge2));
    }
    else {
        stitch_pairs.emplace_back(do_extraction(tableau1.src_id, anchor_chain.front().walk1.front(),
                                                tableau2.src_id, anchor_chain.front().walk2.front(),
                                                graph1, graph2, xmerge1, xmerge2));
        
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
        
        stitch_pairs.emplace_back(do_extraction(anchor_chain.back().walk1.back(), tableau1.snk_id,
                                                anchor_chain.back().walk2.back(), tableau2.snk_id,
                                                graph1, graph2, xmerge1, xmerge2));
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
                                    ChainAlgorithm local_chaining_algorithm,
                                    double anchor_scale,
                                    const std::unordered_set<std::tuple<size_t, size_t, size_t>>* masked_matches) const {
    
    
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
                        size_t orig_idx2 = walks.first[idx2];
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
                    // the counts get retained from the original match sets
                    divvied[stitch_idx].back().count1 = match_set.count1;
                    divvied[stitch_idx].back().count2 = match_set.count2;
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



template<class BGraph, class XMerge>
std::vector<anchor_t> Anchorer::anchor_chain(std::vector<match_set_t>& matches,
                                             const BGraph& graph1,
                                             const BGraph& graph2,
                                             const SentinelTableau& tableau1,
                                             const SentinelTableau& tableau2,
                                             const XMerge& xmerge1,
                                             const XMerge& xmerge2,
                                             std::unordered_set<std::tuple<size_t, size_t, size_t>>* masked_matches) const {
    
    double scale = 1.0;
    if (chaining_algorithm == SparseAffine && autocalibrate_gap_penalties) {
        // this is only to adjust gap penalties, so don't bother if we're not using them
        logging::log(logging::Verbose, "Calibrating gap penalties");
        scale = estimate_score_scale(matches, graph1, graph2, tableau1, tableau2, xmerge1, xmerge2);
        logging::log(logging::Debug, "Estimated score scale: " + std::to_string(scale));
    }
    auto anchors = anchor_chain(matches, graph1, graph2, tableau1, tableau2, xmerge1, xmerge2, chaining_algorithm, false, scale, masked_matches);
    
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
                                      const XMerge& xmerge1, const XMerge& xmerge2) const {
    
    // get an anchoring with unscored gaps
    // FIXME: should i handle masked matches here?
    auto anchors = anchor_chain(matches, graph1, graph2, tableau1, tableau2, xmerge1, xmerge2, Sparse, true, 1.0, nullptr);
    
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
                                             ChainAlgorithm local_chaining_algorithm,
                                             bool suppress_verbose_logging,
                                             double anchor_scale,
                                             std::unordered_set<std::tuple<size_t, size_t, size_t>>* masked_matches) const {
    
    std::vector<anchor_t> anchors;
    if (global_anchoring) {
        anchors = std::move(anchor_chain(matches, graph1, graph2, xmerge1, xmerge2,
                                         &graph1.next(tableau1.src_id), &graph2.next(tableau2.src_id),
                                         &graph1.previous(tableau1.snk_id), &graph2.previous(tableau2.snk_id),
                                         max_num_match_pairs, suppress_verbose_logging, local_chaining_algorithm, anchor_scale,
                                         masked_matches));
    }
    else {
        anchors = std::move(anchor_chain(matches, graph1, graph2, xmerge1, xmerge2,
                                         nullptr, nullptr, nullptr, nullptr,
                                         max_num_match_pairs, suppress_verbose_logging, local_chaining_algorithm, anchor_scale,
                                         masked_matches));
    }
    
    if (do_fill_in_anchoring) {
        fill_in_anchor_chain(anchors, matches, graph1, graph2, tableau1, tableau2, xmerge1, xmerge2,
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
    if (total_num_pairs > local_max_num_match_pairs) {
        // we need to limit the number of nodes
        
        if (!suppress_verbose_logging) {
            logging::log(logging::Debug, "Selecting a maximum of " + std::to_string(local_max_num_match_pairs) + " matches for anchoring");
        }
        
        // prioritize based on the score
        // TODO: adjust this by masked match count?
        auto order = range_vector(matches.size());
        std::stable_sort(order.begin(), order.end(), [&](size_t i, size_t j) {
            return (score_function->anchor_weight(matches[i].count1, matches[i].count2, matches[i].walks1.front().size()) >
                    score_function->anchor_weight(matches[j].count1, matches[j].count2, matches[j].walks1.front().size()));
        });
        
        // greedily choose matches as long as we have budget left
        size_t pairs_left = local_max_num_match_pairs;
        for (size_t i = 0; i < order.size(); ++i) {
            auto& match = matches[order[i]];
            if (score_function->anchor_weight(match.count1, match.count2, match.walks1.front().size()) < 0.0) {
                // these are sorted by score, so nothing else will have positive score
                removed += (order.size() - i);
                break;
            }
            size_t pair_count = match.walks1.size() * match.walks2.size();
            if (pairs_left >= pair_count) {
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
    
    size_t num_match_sets = matches.size() - removed;
    
    // compute the optimal chain using DP
    std::vector<anchor_t> chain;
    switch (local_chaining_algorithm) {
        case SparseAffine:
            chain = std::move(sparse_affine_chain_dp(matches, graph1, graph2, chain_merge1, chain_merge2,
                                                     gap_open, gap_extend, anchor_scale, num_match_sets, suppress_verbose_logging,
                                                     sources1, sources2, sinks1, sinks2, masked_matches));
            break;
            
        case Sparse:
            chain = std::move(sparse_chain_dp(matches, graph1, chain_merge1, chain_merge2, num_match_sets,
                                              suppress_verbose_logging, sources1, sources2, sinks1, sinks2,
                                              masked_matches));
            break;
            
        case Exhaustive:
            chain = std::move(exhaustive_chain_dp(matches, graph1, graph2, chain_merge1, chain_merge2, false,
                                                  anchor_scale, num_match_sets, sources1, sources2, sinks1, sinks2,
                                                  masked_matches));
            break;
            
        default:
            throw std::runtime_error("Unrecognized chaining algorithm: " + std::to_string((int) local_chaining_algorithm));
            break;
    }
    return chain;
}

inline void Anchorer::annotate_scores(std::vector<anchor_t>& anchors) const {
    for (auto& anchor : anchors) {
        anchor.score = anchor_weight(anchor);
    }
}

inline double Anchorer::anchor_weight(const anchor_t& anchor) const {
    return score_function->anchor_weight(anchor.count1, anchor.count2, anchor.walk1.size());
}


template <class BGraph, class XMerge>
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
        switch_dists1 = std::move(post_switch_distances(graph1, chain_merge1));
        switch_dists2 = std::move(post_switch_distances(graph2, chain_merge2));
    }
    
    // make a graph of match nodes
    AnchorGraph anchor_graph;
    for (size_t i = 0; i < num_match_sets; ++i) {
        
        const auto& match_set = match_sets[i];
        
        double weight = score_function->anchor_weight(match_set.count1, match_set.count2,
                                                      match_set.walks1.front().size());
        
        for (size_t idx1 = 0; idx1 < match_set.walks1.size(); ++idx1) {
            for (size_t idx2 = 0; idx2 < match_set.walks2.size(); ++idx2) {
                if (masked_matches && masked_matches->count(std::make_tuple(i, idx1, idx2))) {
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

template<class BGraph, class XMerge>
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
    
    // a key of (chain index, anchor set, walk1 index, walk2 index)
    using key_t = std::tuple<size_t, size_t, size_t, size_t>;
    
    // for each chain2, the initial search tree data for chains, records of (key_t, weight)
    std::vector<std::vector<std::pair<key_t, double>>> search_tree_data(chain_merge2.chain_size());
    
    // for each node, the list of (set, walk1) that start/end on it
    std::vector<std::vector<std::pair<size_t, size_t>>> starts(graph1.node_size()), ends(graph1.node_size());
    
    // for each set, for each walk1, for each walk2, the max weight and (set, walk1, walk2) for the traceback
    static const double mininf = std::numeric_limits<double>::lowest();
    std::vector<std::vector<std::vector<dp_entry_t>>> dp(num_match_sets);
    
    if (debug_anchorer) {
        std::cerr << "gathering anchor endpoint information\n";
    }
    
    // do the bookkeeping
    for (int64_t i = 0; i < dp.size(); ++i) {
        // get the starts and ends of anchors on graph 1
        auto& match_set = match_sets[i];
        for (size_t j = 0; j < match_set.walks1.size(); ++j) {
            // get the starts and ends of anchor on graph 1
            auto& walk = match_set.walks1[j];
            starts[walk.front()].emplace_back(i, j);
            ends[walk.back()].emplace_back(i, j);
            
            // get the chain index of anchor ends on graph 2
            for (size_t k = 0; k < match_set.walks2.size(); ++k) {
                if (masked_matches && masked_matches->count(std::make_tuple(i, j, k))) {
                    continue;
                }
                uint64_t chain2;
                size_t index;
                std::tie(chain2, index) = chain_merge2.chain(match_set.walks2[k].back());
                
                search_tree_data[chain2].emplace_back(key_t(index, i, j, k), mininf);
            }
        }
        
        // initialize the DP structure with a single-anchor chain at each position
        double weight = score_function->anchor_weight(match_set.count1, match_set.count2,
                                                      match_set.walks1.front().size());
        
        
        dp[i].resize(match_set.walks1.size(),
                     std::vector<dp_entry_t>(match_set.walks2.size(),
                                             dp_entry_t(weight, -1, -1, -1)));
        
        if (sources1) {
            // not all nodes can be starts for DP, must be reachable by one of the sources
            for (size_t j = 0; j < match_set.walks1.size(); ++j) {
                bool found1 = false;
                for (auto src_id1 : *sources1) {
                    if (src_id1 == match_set.walks1[j].front() || chain_merge1.reachable(src_id1, match_set.walks1[j].front())) {
                        found1 = true;
                    }
                }
                for (size_t k = 0; k < match_set.walks2.size(); ++k) {
                    bool found2 = false;
                    if (found1) {
                        for (auto src_id2 : *sources2) {
                            if (src_id2 == match_set.walks2[k].front() || chain_merge2.reachable(src_id2, match_set.walks2[k].front())) {
                                found2 = true;
                            }
                        }
                    }
                    if (!found2) {
                        std::get<0>(dp[i][j][k]) = mininf;
                    }
                }
            }
        }
    }
    
    if (debug_anchorer) {
        std::cerr << "initial DP state:\n";
        for (size_t i = 0; i < match_sets.size(); ++i) {
            const auto& table = dp[i];
            for (size_t j = 0; j < table.size(); ++j) {
                const auto& row = table[j];
                for (size_t k = 0; k < row.size(); ++k) {
                    std::cerr << i << ' ' << j << ' ' << k << ' ' << std::get<0>(row[k]) << '\n';
                }
            }
        }
        
        std::cerr << "initializing search trees\n";
    }
    
    // sort one time to speed up construction across chains
    // TODO: we could use integer sort to exchange O(n log n) with O(V)
    for (auto& chain_search_tree_data : search_tree_data) {
        std::stable_sort(chain_search_tree_data.begin(), chain_search_tree_data.end());
    }
    
    if (debug_anchorer) {
        for (size_t i = 0; i < search_tree_data.size(); ++i) {
            std::cerr << "search tree keys for chain " << i << " of graph2:\n";
            for (auto& key_val_pair : search_tree_data[i]) {
                std::cerr << '\t' << std::get<0>(key_val_pair.first) << ' ' << std::get<1>(key_val_pair.first) << ' ' << std::get<2>(key_val_pair.first) << ' ' << std::get<3>(key_val_pair.first) << '\n';
            }
        }
    }
    
    // for each chain1, for each chain2, a tree over seed end chain indexes
    std::vector<std::vector<MaxSearchTree<key_t, double>>> search_trees;
    search_trees.resize(chain_merge1.chain_size());
    
    for (size_t i = 0; i < chain_merge1.chain_size(); ++i) {
        auto& chain_search_trees = search_trees[i];
        chain_search_trees.reserve(chain_merge2.chain_size());
        for (size_t j = 0; j < search_tree_data.size(); ++j) {
            chain_search_trees.emplace_back(search_tree_data[j]);
        }
    }
    // clear the search tree data
    {
        auto dummy = std::move(search_tree_data);
    }
    
    if (debug_anchorer) {
        std::cerr << "computing forward edges\n";
    }
    
    // get the edges to chain neighbors
    // TODO: we could thin these out based on which forward edges are actually needed based on the
    // anchor positions
    auto forward_edges = chain_merge1.chain_forward_edges();
    
    
    if (debug_anchorer) {
        std::cerr << "beginning main DP iteration\n";
    }
    
    size_t iter = 0;
    for (uint64_t node_id : topological_order(graph1)) {
        
        ++iter;
        if (logging::level >= logging::Verbose && iter % 1000000 == 0 && !suppress_verbose_logging) {
            logging::log(logging::Verbose, "Iteration " + std::to_string(iter) + " of " + std::to_string(graph1.node_size()) + " in sparse chaining algorithm");
        }
        
        if (debug_anchorer) {
            std::cerr << "on node " << node_id << '\n';
        }
        
        for (const auto& end : ends[node_id]) {
            // we've hit the end of this match, so we can enter its DP value into the trees
            // to update other DP values as we move forward
            if (debug_anchorer) {
                std::cerr << "node is end of graph 1 anchor " << end.first << ',' << end.second << '\n';
            }
            
            auto& match_set = match_sets[end.first];
            uint64_t chain1 = chain_merge1.chain(node_id).first;
            auto& dp_row = dp[end.first][end.second];
            
            // add a value in the appropriate tree for each occurrences in graph2
            for (size_t i = 0; i < match_set.walks2.size(); ++i) {
                if (masked_matches && masked_matches->count(std::make_tuple(end.first, end.second, i))) {
                    continue;
                }
                const auto& walk2 = match_set.walks2[i];
                uint64_t chain2;
                size_t index;
                std::tie(chain2, index) = chain_merge2.chain(walk2.back());
                
                auto& tree = search_trees[chain1][chain2];
                auto it = tree.find(key_t(index, end.first, end.second, i));
                if (it->second < std::get<0>(dp_row[i])) {
                    tree.update(it, std::get<0>(dp_row[i]));
                    if (debug_anchorer) {
                        std::cerr << "recording increased DP value of " << it->second << " on chain " << chain2 << " with key " << std::get<0>(it->first) << ' ' << std::get<1>(it->first) << ' ' << std::get<2>(it->first) << ' ' << std::get<3>(it->first) << '\n';
                    }
                }
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
            
            for (const auto& start : starts[fwd_id]) {
                
                // an anchor starts here in graph1
                
                if (debug_anchorer) {
                    std::cerr << "anchor " << start.first << ',' << start.second << " starts on " << fwd_id << '\n';
                }
                
                const auto& match_set = match_sets[start.first];
                auto& dp_row = dp[start.first][start.second];
                
                // the weight of this anchors in this set
                double weight = score_function->anchor_weight(match_set.count1, match_set.count2,
                                                              match_set.walks1.front().size());
                
                // we will consider all occurrences of this anchor in graph2
                for (size_t j = 0; j < match_set.walks2.size(); ++j) {
                    
                    if (masked_matches && masked_matches->count(std::make_tuple(start.first, start.second, j))) {
                        continue;
                    }
                    
                    auto& dp_entry = dp_row[j];
                    
                    if (debug_anchorer) {
                        std::cerr << "walk2 index " << j << " of " << match_set.walks2.size() << " starts on node " << match_set.walks2[j].front() << " and has current DP entry:\n";
                        std::cerr << std::get<0>(dp_entry) << ' ' << std::get<1>(dp_entry) << ' ' << std::get<2>(dp_entry) << ' ' << std::get<3>(dp_entry) << '\n';
                    }
                    // we will check for previous DP values that can reach this one in graph2 from
                    // each of the chains in the chain partition of graph2
                    const auto& chain_preds2 = chain_merge2.predecessor_indexes(match_set.walks2[j].front());
                    for (uint64_t chain2 = 0; chain2 < chain_preds2.size(); ++chain2) {
                        
                        if (chain_preds2[chain2] == -1) {
                            // there is no reachable node from this chain
                            continue;
                        }
                        if (debug_anchorer) {
                            std::cerr << "looking for predecessor on chain " << chain2 << " with predecessor at or before index " << chain_preds2[chain2] << '\n';
                        }
                        // find the max DP value up to (and including) the predecessor
                        const auto& tree = search_trees[chain1][chain2];
                        
                        auto it = tree.range_max(key_t(0, 0, 0, 0),
                                                 key_t(chain_preds2[chain2] + 1, 0, 0, 0)); // +1 because past-the-last
                        
                        if (it == tree.end()) {
                            // there aren't any ends of anchors before the predecessor in this chain
                            if (debug_anchorer) {
                                std::cerr << "there are no predecessors on this chain\n";
                            }
                            continue;
                        }
                        
                        // the weight of this anchor plus all previous anchors
                        double dp_weight = it->second + weight;
                        if (dp_weight > std::get<0>(dp_entry)) {
                            // update the DP values and the backpointer
                            if (debug_anchorer) {
                                std::cerr << "got increased DP weight of " << dp_weight << " along chain " << chain2 << " from " << std::get<1>(it->first) << ' ' << std::get<2>(it->first) << ' ' << std::get<3>(it->first) << '\n';
                            }
                            dp_entry = dp_entry_t(dp_weight, std::get<1>(it->first), std::get<2>(it->first), std::get<3>(it->first));
                        }
                        else if (debug_anchorer) {
                            std::cerr << "DP weight of " << dp_weight << " from predecessor " << std::get<1>(it->first) << ' ' << std::get<2>(it->first) << ' ' << std::get<3>(it->first) << " does not beat current DP value of " << std::get<0>(dp_entry) << '\n';
                        }
                    }
                }
            }
        }
    }
    
    auto final_term = [&](size_t i, size_t j, size_t k) -> double {
        if (sinks1) {
            const auto& match_set = match_sets[i];
            uint64_t last_id1 = match_set.walks1[j].back();
            uint64_t last_id2 = match_set.walks2[k].back();
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
    
    return traceback_sparse_dp(match_sets, dp, final_term, 0.0, suppress_verbose_logging, masked_matches);
}



template<class BGraph, class XMerge>
std::vector<std::vector<size_t>> Anchorer::post_switch_distances(const BGraph& graph, const XMerge& xmerge) const {
    
    // TODO: use label size instead of assuming node length 1?
    
    std::vector<std::vector<size_t>> dists(graph.node_size(),
                                           std::vector<size_t>(xmerge.chain_size(), -1));
    
    // note: this DP is different from the paper because i have the predecessor on the same
    // path being the previous node rather than the node itself
    for (auto node_id : topological_order(graph)) {
        auto& row = dists[node_id];
        const auto& preds = xmerge.predecessor_indexes(node_id);
        for (uint64_t p = 0; p < xmerge.chain_size(); ++p) {
            for (auto prev_id : graph.previous(node_id)) {
                if (xmerge.index_on(prev_id, p) == preds[p]) {
                    // switching paths you here immediately, no distance
                    row[p] = 0;
                    break;
                }
                else if (xmerge.predecessor_indexes(prev_id)[p] == preds[p]) {
                    // travel through this predecessor after switching onto it from path p
                    // note: this will overwrite the default -1
                    row[p] = std::min(row[p], dists[prev_id][p] + 1);
                }
            }
        }
    }
    
    return dists;
}

template<size_t NumPW, class BGraph, class XMerge>
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
        for (size_t i = 0; i < match_sets.size(); ++i) {
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
    
    // (offset diff, match set, walk1, walk2)
    using key_t = std::tuple<int64_t, size_t, size_t, size_t>;
    
    // the D arrays from chandra & jain 2023
    auto switch_dists1 = post_switch_distances(graph1, xmerge1);
    auto switch_dists2 = post_switch_distances(graph2, xmerge2);
    
    // the gap contribution from a source anchor
    auto basic_source_shift = [&](uint64_t src_id1, uint64_t src_id2, uint64_t path1, uint64_t path2) -> int64_t {
        return xmerge1.index_on(src_id1, path1) - xmerge2.index_on(src_id2, path2);
    };
    auto source_shift = [&](size_t i, size_t j, size_t k, uint64_t path1, uint64_t path2) -> int64_t {
        const auto& match_set = match_sets[i];
        return basic_source_shift(match_set.walks1[j].back(), match_set.walks2[k].back(), path1, path2);
    };
    // add the source coordinates to grab as backpointers
    auto get_key = [&](size_t i, size_t j, size_t k, uint64_t path1, uint64_t path2) -> key_t {
        return key_t(source_shift(i, j, k, path1, path2), i, j, k);
    };
    // the gap contribution from the destination anchor
    auto basic_query_shift = [&](uint64_t query_id1, uint64_t query_id2, uint64_t path1, uint64_t path2) -> int64_t {
        return (xmerge1.predecessor_indexes(query_id1)[path1] - xmerge2.predecessor_indexes(query_id2)[path2]
                + switch_dists1[query_id1][path1] - switch_dists2[query_id2][path2]);
    };
    auto query_shift = [&](size_t i, size_t j, size_t k, uint64_t path1, uint64_t path2) -> int64_t {
        const auto& match_set = match_sets[i];
        return basic_query_shift(match_set.walks1[j].front(), match_set.walks2[k].front(), path1, path2);
    };
    // the offset on path from graph2
    auto get_key_offset = [&](size_t i, size_t k, uint64_t path2) -> size_t {
        return xmerge2.index_on(match_sets[i].walks2[k].back(), path2);
    };
    // the "effective offset" on path of graph2 (true offsets on the chain before this can reach it)
    auto get_query_offset = [&](size_t i, size_t k, uint64_t path2) -> size_t {
        // note: we rely on -1's overflowing to 0
        return xmerge2.predecessor_indexes(match_sets[i].walks2[k].front())[path2] + 1;
    };
    
    // for each node, the list of (set, walk1) that start/end on it
    std::vector<std::vector<std::pair<size_t, size_t>>> starts(graph1.node_size()), ends(graph1.node_size());
    
    // for each set, for each walk1, for each walk2, the max weight and (set, walk1, walk2) for the traceback
    static const double mininf = std::numeric_limits<double>::lowest();
    std::vector<std::vector<std::vector<dp_entry_t>>> dp(num_match_sets);
        
    // for each path1, for each path2, list of (search index, value)
    std::vector<std::vector<std::vector<std::tuple<key_t, size_t, double>>>> search_tree_data;
    search_tree_data.resize(xmerge1.chain_size(), std::vector<std::vector<std::tuple<key_t, size_t, double>>>(xmerge2.chain_size()));
    
    if (debug_anchorer) {
        std::cerr << "doing bookkeeping for sparse affine algorithm\n";
    }
    
    uint64_t num_pairs = 0;
    
    // do the bookkeeping
    for (int64_t i = 0; i < dp.size(); ++i) {
        // get the starts and ends of anchors on graph 1
        auto& match_set = match_sets[i];
        num_pairs += match_set.walks1.size() * match_set.walks2.size();
        for (size_t j = 0; j < match_set.walks1.size(); ++j) {
            // get the starts and ends of anchor on graph 1
            auto& walk1 = match_set.walks1[j];
            starts[walk1.front()].emplace_back(i, j);
            ends[walk1.back()].emplace_back(i, j);
            
            auto path1s = xmerge1.chains_on(walk1.back());
            for (size_t k = 0; k < match_set.walks2.size(); ++k) {
                if (masked_matches && masked_matches->count(std::make_tuple(i, j, k))) {
                    continue;
                }
                for (auto p1 : path1s) {
                    for (auto p2 : xmerge2.chains_on(match_set.walks2[k].back())) {
                                                
                        // initialize the sparse value data
                        search_tree_data[p1][p2].emplace_back(get_key(i, j, k, p1, p2),
                                                              get_key_offset(i, k, p2),
                                                              mininf);
                    }
                }
            }
        }
        
        // initialize the DP structure with a single-anchor chain at each position
        double weight = score_function->anchor_weight(match_set.count1, match_set.count2,
                                                      match_set.walks1.front().size());
        
        dp[i].resize(match_set.walks1.size(),
                     std::vector<dp_entry_t>(match_set.walks2.size(),
                                             dp_entry_t(weight, -1, -1, -1)));
        if (sources1) {
            // we are doing global anchoring, so we also have to pay for the lead indel
            for (size_t j = 0; j < match_set.walks1.size(); ++j) {
                auto& walk1 = match_set.walks1[j];
                for (size_t k = 0; k < match_set.walks2.size(); ++k) {
                    auto& walk2 = match_set.walks2[k];
                    
                    // TODO: could i get rid of these loops by querying the structure somehow?
                    
                    int64_t lead_indel = std::numeric_limits<int64_t>::max();
                    
                    // iterate over combos of source ids and their paths
                    for (auto src_id1 : *sources1) {
                        if (walk1.front() != src_id1 && !xmerge1.reachable(src_id1, walk1.front())) {
                            continue;
                        }
                        for (auto src_id2 : *sources2) {
                            if (walk2.front() != src_id2 && !xmerge2.reachable(src_id2, walk2.front())) {
                                continue;
                            }
                            for (auto p1 : xmerge1.chains_on(src_id1)) {
                                for (auto p2 : xmerge2.chains_on(src_id2)) {
                                    lead_indel = std::min<int64_t>(lead_indel, std::abs(basic_source_shift(src_id1, src_id2, p1, p2)
                                                                                        - query_shift(i, j, k, p1, p2)));
                                    
                                }
                            }
                        }
                    }
                    
                    if (lead_indel == std::numeric_limits<int64_t>::max()) {
                        // this anchor was not reachable
                        std::get<0>(dp[i][j][k]) = mininf;
                    }
                    else if (lead_indel != 0) {
                        double lead_indel_weight = mininf;
                        for (int pw = 0; pw < NumPW; ++pw) {
                            lead_indel_weight = std::max(lead_indel_weight, -local_scale * (gap_open[pw] + gap_extend[pw] * lead_indel));
                        }
                        std::get<0>(dp[i][j][k]) += lead_indel_weight;
                    }
                }
            }
        }
    }
    
    if (debug_anchorer) {
        std::cerr << "initial DP state:\n";
        for (size_t i = 0; i < match_sets.size(); ++i) {
            const auto& table = dp[i];
            for (size_t j = 0; j < table.size(); ++j) {
                const auto& row = table[j];
                for (size_t k = 0; k < row.size(); ++k) {
                    std::cerr << i << ' ' << j << ' ' << k << ' ' << std::get<0>(row[k]) << " (" << score_function->anchor_weight(match_sets[i].count1, match_sets[i].count2, match_sets[i].walks1.front().size()) << ")\n";
                }
            }
        }
    }
    
    if (!suppress_verbose_logging) {
        logging::log(logging::Verbose, "Chaining " + std::to_string(num_pairs) + " matches");
    }
    
    if (debug_anchorer) {
        std::cerr << "constructing search trees\n";
    }
    
    if (!suppress_verbose_logging) {
        logging::log(logging::Debug, "Initializing sparse query data structures");
    }
    
    // grids of search trees over (path1, path2) where scores correspond to different distance scenaries
    // 0           d1 = d2
    // odds        d1 > d2
    // evens > 0   d1 < d2
    std::array<std::vector<std::vector<OrthogonalMaxSearchTree<key_t, size_t, double>>>, 2 * NumPW + 1> search_trees;
    for (uint64_t pw = 0; pw < 2 * NumPW + 1; ++pw) {
        auto& pw_trees = search_trees[pw];
        pw_trees.resize(xmerge1.chain_size());
        for (uint64_t p1 = 0; p1 < xmerge1.chain_size(); ++p1) {
            auto& tree_row = pw_trees[p1];
            tree_row.reserve(xmerge2.chain_size());
            for (uint64_t p2 = 0; p2 < xmerge2.chain_size(); ++p2) {
                tree_row.emplace_back(search_tree_data[p1][p2]);
            }
        }
    }
    
    // clear the search tree data
    {
        auto dummy = std::move(search_tree_data);
    }
    
    if (debug_anchorer) {
        std::cerr << "computing forward edges\n";
    }
    
    // get the edges to chain neighbors
    // TODO: we could thin these out based on which forward edges are actually needed based on the
    // anchor positions
    auto forward_edges = xmerge1.chain_forward_edges();
    
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
            logging::log(logging::Verbose, "Iteration " + std::to_string(iter) + " of " + std::to_string(graph1.node_size()) + " in sparse chaining algorithm");
        }
        
        if (debug_anchorer) {
            std::cerr << "on node " << node_id << '\n';
        }
        
        for (const auto& end : ends[node_id]) {
            // we've hit the end of this match, so we can enter its DP value into the trees
            // to update other DP values as we move forward
            if (debug_anchorer) {
                std::cerr << "node is end of graph 1 anchor " << end.first << ',' << end.second << '\n';
            }
            
            auto& match_set = match_sets[end.first];
            auto& dp_row = dp[end.first][end.second];
            
            for (size_t k = 0; k < match_set.walks2.size(); ++k) {
                if (masked_matches && masked_matches->count(std::make_tuple(end.first, end.second, k))) {
                    continue;
                }
                const auto& dp_val = std::get<0>(dp[end.first][end.second][k]);
                if (debug_anchorer) {
                    std::cerr << "considering matched graph 2 anchor " << k << " with DP value " << dp_val << '\n';
                }
                for (auto p2 : xmerge2.chains_on(match_set.walks2[k].back())) {
                    auto key2 = get_key_offset(end.first, k, p2);
                    for (auto p1 : xmerge1.chains_on(match_set.walks1[end.second].back())) {
                        auto key1 = get_key(end.first, end.second, k, p1, p2);
                        if (debug_anchorer) {
                            std::cerr << "extending for path combo " << p1 << ',' << p2 << ", giving shift key " << std::get<0>(key1) << " and offset key " << key2 << '\n';
                        }
                        for (size_t pw = 0; pw < 2 * NumPW + 1; ++pw) {
                            // save the anchor-independent portion of the score in the search tree
                            double value;
                            if (pw == 0) {
                                // d1 == d2
                                value = dp_val;
                            }
                            else if (pw % 2 == 1) {
                                // d1 > d2
                                value = dp_val + local_scale * gap_extend[pw / 2] * source_shift(end.first, end.second, k, p1, p2);
                            }
                            else {
                                // d1 < d2
                                value = dp_val - local_scale * gap_extend[pw / 2 - 1] * source_shift(end.first, end.second, k, p1, p2);
                            }
                            auto& tree = search_trees[pw][p1][p2];
                            auto it = tree.find(key1, key2);
                            if (value > std::get<2>(*it)) {
                                // TODO: shouldn't this condition always be met, since keys are unique to a match pair?
                                if (debug_anchorer) {
                                    std::cerr << "register " << value << " for key " << std::get<0>(key1) << ',' << std::get<1>(key1) << ',' << std::get<2>(key1) << ',' << std::get<3>(key1) << ", offset " << key2 << " in piecewise component " << pw << "\n";
                                }
                                tree.update(it, value);
                            }
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
            
            for (const auto& start : starts[fwd_id]) {
                
                // an anchor starts here in graph1
                
                if (debug_anchorer) {
                    std::cerr << "graph1 anchor " << start.first << ',' << start.second << " starts on " << fwd_id << '\n';
                }
                
                const auto& match_set = match_sets[start.first];
                auto& dp_row = dp[start.first][start.second];
                
                // the weight of this anchors in this set
                double weight = score_function->anchor_weight(match_set.count1, match_set.count2,
                                                              match_set.walks1.front().size());
                
                // we will consider all occurrences of this anchor in graph2
                for (size_t k = 0; k < match_set.walks2.size(); ++k) {
                    if (masked_matches && masked_matches->count(std::make_tuple(start.first, start.second, k))) {
                        continue;
                    }
                    auto& dp_entry = dp_row[k];
                    if (debug_anchorer) {
                        std::cerr << "checking graph2 anchor " << k << " with DP value " << std::get<0>(dp_entry) << '\n';
                    }
                    for (uint64_t chain2 = 0; chain2 < xmerge2.chain_size(); ++chain2) {
                        // note: we have to check all of the chains because the best distance measure might
                        // not originate from a path that contains the head of the path
                        
                        int64_t query = query_shift(start.first, start.second, k, chain1, chain2);
                        size_t offset = get_query_offset(start.first, k, chain2);
                        if (debug_anchorer) {
                            std::cerr << "query shift is " << query << " and offset is " << offset << " on chain combo " << chain1 << "," << chain2 << '\n';
                        }
                        for (size_t pw = 0; pw < 2 * NumPW + 1; ++pw) {
                            // combine the anchor-dependent and anchor-independent portions of the score
                            const auto& tree = search_trees[pw][chain1][chain2];
                            if (pw == 0) {
                                // d1 = d2, search only at this query value
                                auto it = tree.range_max(key_t(query, 0, 0, 0),
                                                         key_t(query + 1, 0, 0, 0),
                                                         0, offset);
                                // note: nodes can query themselves here, but only before their true value is included in the
                                // search trees (it is just the placeholder min inf)
                                if (it != tree.end()) {
                                    double value = std::get<2>(*it) + weight;
                                    if (debug_anchorer) {
                                        auto key1 = std::get<0>(*it);
                                        std::cerr << "piecewise component " << pw << " got source hit " << std::get<0>(key1) << ',' << std::get<1>(key1) << ',' << std::get<2>(key1) << ',' << std::get<3>(key1) << ": " << std::get<2>(*it) << ", extends to score " << value << '\n';
                                    }
                                    if (value > std::get<0>(dp_entry)) {
                                        if (debug_anchorer) {
                                            std::cerr << "this hit is the new opt\n";
                                        }
                                        dp_entry = dp_entry_t(value, std::get<1>(std::get<0>(*it)),
                                                              std::get<2>(std::get<0>(*it)), std::get<3>(std::get<0>(*it)));
                                    }
                                }
                            }
                            else if (pw % 2 == 1) {
                                // d1 > d2, search leftward of the query value
                                auto it = tree.range_max(key_t(std::numeric_limits<int64_t>::min(), 0, 0, 0),
                                                         key_t(query, 0, 0, 0),
                                                         0, offset);
                                if (it != tree.end()) {
                                    double value = std::get<2>(*it) + weight - local_scale * (gap_open[pw / 2] + gap_extend[pw / 2] * query);
                                    if (debug_anchorer) {
                                        auto key1 = std::get<0>(*it);
                                        std::cerr << "piecewise component " << pw << " got source hit " << std::get<0>(key1) << ',' << std::get<1>(key1) << ',' << std::get<2>(key1) << ',' << std::get<3>(key1) << ": " << std::get<2>(*it) << ", extends to score " << value << '\n';
                                    }
                                    if (value > std::get<0>(dp_entry)) {
                                        if (debug_anchorer) {
                                            std::cerr << "this hit is the new opt\n";
                                        }
                                        dp_entry = dp_entry_t(value, std::get<1>(std::get<0>(*it)),
                                                              std::get<2>(std::get<0>(*it)), std::get<3>(std::get<0>(*it)));
                                    }
                                }
                            }
                            else {
                                // d1 < d2, search right of the query value
                                auto it = tree.range_max(key_t(query + 1, 0, 0, 0),
                                                         key_t(std::numeric_limits<int64_t>::max(),
                                                               std::numeric_limits<size_t>::max(),
                                                               std::numeric_limits<size_t>::max(),
                                                               std::numeric_limits<size_t>::max()),
                                                         0, offset);
                                if (it != tree.end()) {
                                    double value = std::get<2>(*it) + weight - local_scale * (gap_open[pw / 2 - 1] - gap_extend[pw / 2 - 1] * query);
                                    if (debug_anchorer) {
                                        auto key1 = std::get<0>(*it);
                                        std::cerr << "piecewise component " << pw << " got source hit " << std::get<0>(key1) << ',' << std::get<1>(key1) << ',' << std::get<2>(key1) << ',' << std::get<3>(key1) << ": " << std::get<2>(*it) << ", extends to score " << value << '\n';
                                    }
                                    if (value > std::get<0>(dp_entry)) {
                                        if (debug_anchorer) {
                                            std::cerr << "this hit is the new opt\n";
                                        }
                                        dp_entry = dp_entry_t(value, std::get<1>(std::get<0>(*it)),
                                                              std::get<2>(std::get<0>(*it)), std::get<3>(std::get<0>(*it)));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    double min_score = 0.0;
    if (sources1 && sinks1) {
        // we have to do better than the empty chain, so let's figure out its score
        
        min_score = mininf;
        
        // find minimum length indel
        int64_t min_indel = std::numeric_limits<int64_t>::max();
        for (auto src_id1 : *sources1) {
            for (auto src_id2 : *sources2) {
                for (auto p1 : xmerge1.chains_on(src_id1)) {
                    for (auto p2 : xmerge2.chains_on(src_id2)) {
                        int64_t src_shift = basic_source_shift(src_id1, src_id2, p1, p2);
                        for (auto snk_id1 : *sinks1) {
                            for (auto snk_id2 : *sinks2) {
                                int64_t shift = std::abs(src_shift - basic_query_shift(snk_id1, snk_id2, p1, p2));
                                min_indel = std::min(min_indel, shift);
                            }
                        }
                    }
                }
            }
        }
        
        // score it
        if (min_indel == 0) {
            min_score = 0.0;
        }
        else {
            for (int pw = 0; pw < NumPW; ++pw) {
                min_score = std::max(min_score, -local_scale * (gap_open[pw] + gap_extend[pw] * min_indel));
            }
        }
    }
    
    if (debug_anchorer) {
        std::cerr << "traceback must exceed " << min_score << " to be counted\n";
    }
    
    // the score for the final indel
    auto final_indel_score = [&](size_t i, size_t j, size_t k) -> double {
        if (!sinks1) {
            return 0.0;
        }
        
        const auto& match_set = match_sets[i];
        const auto& walk1 = match_set.walks1[j];
        const auto& walk2 = match_set.walks2[k];
        
        // find the minimum length indel
        int64_t final_indel = std::numeric_limits<int64_t>::max();
        for (auto p1 : xmerge1.chains_on(walk1.back())) {
            for (auto p2 : xmerge2.chains_on(walk2.back())) {
                
                int64_t src_shift = basic_source_shift(walk1.back(), walk2.back(), p1, p2);
                for (auto snk_id1 : *sinks1) {
                    if (walk1.back() != snk_id1 && !xmerge1.reachable(walk1.back(), snk_id1)) {
                        continue;
                    }
                    for (auto snk_id2 : *sinks2) {
                        if (walk2.back() != snk_id2 && !xmerge2.reachable(walk2.back(), snk_id2)) {
                            continue;
                        }
                        int64_t shift = std::abs(src_shift - basic_query_shift(snk_id1, snk_id2, p1, p2));
                        final_indel = std::min<int64_t>(final_indel, shift);
                        
                    }
                }
            }
        }
        
        // score it
        if (final_indel == std::numeric_limits<int64_t>::max()) {
            return mininf;
        }
        else if (final_indel == 0) {
            return 0.0;
        }
        else {
            double score = mininf;
            for (int pw = 0; pw < NumPW; ++pw) {
                score = std::max(score, -local_scale * (gap_open[pw] + gap_extend[pw] * final_indel));
            }
            return score;
        }
    };
        
    return traceback_sparse_dp(match_sets, dp, final_indel_score, min_score, suppress_verbose_logging, masked_matches);
}

template<class FinalFunc>
std::vector<anchor_t> Anchorer::traceback_sparse_dp(const std::vector<match_set_t>& match_sets,
                                                    const std::vector<std::vector<std::vector<dp_entry_t>>>& dp,
                                                    const FinalFunc& final_function, double min_score,
                                                    bool suppress_verbose_logging,
                                                    const std::unordered_set<std::tuple<size_t, size_t, size_t>>* masked_matches) const {
    
    if (debug_anchorer) {
        std::cerr << "finding optimum\n";
    }
    
    // find the optimum dynamic programming values
    dp_entry_t opt(std::numeric_limits<double>::lowest(), -1, -1, -1);
    for (size_t set = 0; set < dp.size(); ++set) {
        const auto& set_dp = dp[set];
        for (size_t i = 0; i < set_dp.size(); ++i) {
            const auto& dp_row = set_dp[i];
            for (size_t j = 0; j < dp_row.size(); ++j) {
                if (masked_matches && masked_matches->count(std::make_tuple(set, i, j))) {
                    continue;
                }
                
                double final_term = final_function(set, i, j);
                if (debug_anchorer && final_term == std::numeric_limits<double>::lowest()) {
                    std::cerr << "no valid final indel term for " << set << " " << i << " " << j << '\n';
                }
                
                if (final_term != std::numeric_limits<double>::lowest()) {
                    double score = std::get<0>(dp_row[j]) + final_term;
                    if (debug_anchorer) {
                        std::cerr << "traceback at " << set << " " << i << " " << j << ": " << score << " (" << std::get<0>(dp_row[j]) << " dp, " << final_term << " final)\n";
                    }
                    if (score > std::get<0>(opt) && score > min_score) {
                        opt = dp_entry_t(score, set, i, j);
                    }
                }
            }
        }
    }
    
    if (debug_anchorer) {
        std::cerr << "doing traceback from opt score " << std::get<0>(opt) << ": " << std::get<1>(opt) << " " << std::get<2>(opt) << " " << std::get<3>(opt) << "\n";
    }
    
    // traceback into a chain
    std::vector<anchor_t> anchors;
    auto here = opt;
    while (std::get<1>(here) != -1) {
        
        if (debug_anchorer) {
            std::cerr << "following traceback to set " << std::get<1>(here) << ", walk pair " << std::get<2>(here) << " " << std::get<3>(here) << " with value " << std::get<0>(here) << '\n';
        }
        
        // grab the anchors that we used from their set
        auto& match_set = match_sets[std::get<1>(here)];
        anchors.emplace_back();
        auto& anchor = anchors.back();
        anchor.walk1 = match_set.walks1[std::get<2>(here)];
        anchor.count1 = match_set.count1;
        anchor.walk2 = match_set.walks2[std::get<3>(here)];
        anchor.count2 = match_set.count2;
        anchor.match_set = std::get<1>(here);
        anchor.idx1 = std::get<2>(here);
        anchor.idx2 = std::get<3>(here);
        
        // follow the backpointer from the DP structure
        here = dp[std::get<1>(here)][std::get<2>(here)][std::get<3>(here)];
    }
    
    // take out of reverse order
    std::reverse(anchors.begin(), anchors.end());
    
    annotate_scores(anchors);
    
    if (debug_anchorer) {
        std::cerr << "completed sparse chaining\n";
    }
    
    if (!suppress_verbose_logging) {
        logging::log(logging::Debug, "Optimal chain consists of " + std::to_string(anchors.size()) + " matches with score " + std::to_string(std::get<0>(opt)));
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
            int64_t dist1 = (xmerge1.predecessor_indexes(to_id1)[chain1]
                             - xmerge1.index_on(from_id1, chain1) + switch_dists1[to_id1][chain1]);
            int64_t dist2 = (xmerge2.predecessor_indexes(to_id2)[chain2]
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
        switch_dists1 = std::move(post_switch_distances(graph1, xmerge1));
        switch_dists2 = std::move(post_switch_distances(graph2, xmerge2));
    }
    
    for (size_t i = 0; i < chain.size(); ++i) {
        std::cerr << '@' << '\t' << i << '\t' << chain[i].walk1.size() << '\t' << chain[i].count1 << '\t' << chain[i].count2 << '\t' << (chain[i].count1 * chain[i].count2) << '\t' << score_function->anchor_weight(chain[i].count1, chain[i].count2, chain[i].walk1.size()) << '\t';
        if (chaining_algorithm != SparseAffine || i == 0) {
            std::cerr << 0.0;
        }
        else {
            std::cerr << edge_weight(chain[i-1].walk1.back(), chain[i].walk1.front(),
                                     chain[i-1].walk2.back(), chain[i].walk2.front(),
                                     local_scale,
                                     xmerge1, xmerge2, switch_dists1, switch_dists2);
        }
        std::cerr << '\t' << chain[i].walk1.front() << '\t' << chain[i].walk1.back() << '\t' << chain[i].walk2.front() << '\t' << chain[i].walk2.back();
        size_t size1, size2;
        if (i == 0) {
            size1 = chain[i].walk1.front();
            size2 = chain[i].walk2.front();
            std::cerr << '\t' << chain[i].walk1.front() << '\t' << chain[i].walk2.front();
        }
        else {
            size1 = (chain[i].walk1.front() - chain[i - 1].walk1.back());
            size2 = (chain[i].walk2.front() - chain[i - 1].walk2.back());
        }
        std::vector<double> sizes{(double) size1, (double) size2};
        std::cerr << '\t' << size1 << '\t' << size2 << '\t' << generalized_mean(sizes.begin(), sizes.end(), -0.5) << '\n';
    }
}


}

#endif /* centrolign_anchorer_hpp */
