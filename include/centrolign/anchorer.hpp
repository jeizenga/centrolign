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
};

/*
 * Data structure finding anchors between two graphs
 */
class Anchorer {
public:
    Anchorer() = default;
    ~Anchorer() = default;
    
    enum ChainAlgorithm { Exhaustive = 0, Sparse = 1, SparseAffine = 2 };
    
    // compute a heaviest weight anchoring of a set of matches
    // anchoring is "local" in that it can start at any match and end at
    // any match
    // may reorder matches
    template<class BGraph, class XMerge>
    std::vector<anchor_t> anchor_chain(std::vector<match_set_t>& matches,
                                       const BGraph& graph1,
                                       const BGraph& graph2,
                                       const XMerge& chain_merge1,
                                       const XMerge& chain_merge2,
                                       double anchor_scale = 1.0) const;
    // same as previous but with a non-default chaining algorithm
    template<class BGraph, class XMerge>
    std::vector<anchor_t> anchor_chain(std::vector<match_set_t>& matches,
                                       const BGraph& graph1,
                                       const BGraph& graph2,
                                       const XMerge& chain_merge1,
                                       const XMerge& chain_merge2,
                                       ChainAlgorithm chaining_algorithm_override) const;
    
    // compute a heaviest weight anchoring of a set of matches
    // anchoring is "global" in that the first much must be reachable by sources,
    // the last match must be able to reach sinks, and (if scoring indels) a penalty
    // is incurred from the transition from source to match and from match to sink
    // may reorder matches
    template<class BGraph, class XMerge>
    std::vector<anchor_t> global_anchor_chain(std::vector<match_set_t>& matches,
                                              const BGraph& graph1,
                                              const BGraph& graph2,
                                              const XMerge& chain_merge1,
                                              const XMerge& chain_merge2,
                                              const std::vector<uint64_t>& sources1,
                                              const std::vector<uint64_t>& sources2,
                                              const std::vector<uint64_t>& sinks1,
                                              const std::vector<uint64_t>& sinks2,
                                              double anchor_scale = 1.0,
                                              size_t max_num_match_pairs_override = -1) const;
    
    // same as previous but with a non-default chaining algorithm
    template<class BGraph, class XMerge>
    std::vector<anchor_t> global_anchor_chain(std::vector<match_set_t>& matches,
                                              const BGraph& graph1,
                                              const BGraph& graph2,
                                              const XMerge& chain_merge1,
                                              const XMerge& chain_merge2,
                                              const std::vector<uint64_t>& sources1,
                                              const std::vector<uint64_t>& sources2,
                                              const std::vector<uint64_t>& sinks1,
                                              const std::vector<uint64_t>& sinks2,
                                              ChainAlgorithm chaining_algorithm_override,
                                              size_t max_num_match_pairs_override = -1) const;
    
    /*
     * Configurable parameters
     */
    
    // select which algorithm to use
    ChainAlgorithm chaining_algorithm = Sparse;
    // anchor weight is proportional to length
    bool length_scale = true;
    // power to raise the pair count to in the weight function
    double pair_count_power = 1.0;
    // pair count at which count penalty becomes negative
    // TODO: should this be lower to induce more penalty?
    double count_penalty_threshold = 32.0;
    
    double global_scale = 1.0;
    // affine gap parameters
    std::array<double, 3> gap_open{1.0, 100.0, 4000.0};
    std::array<double, 3> gap_extend{2.0, 0.5, 0.001};
    // the max number of match pairs we will use for anchoring
    size_t max_num_match_pairs = 1000000;

protected:
    
    
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
                                       size_t max_num_match_pairs_override,
                                       bool suppress_verbose_logging,
                                       ChainAlgorithm chaining_algorithm_override,
                                       double anchor_scale) const;

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
                                              const std::vector<uint64_t>* sinks2 = nullptr) const;
    
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
                                          const std::vector<uint64_t>* sinks2 = nullptr) const;
    
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
                                                 const std::vector<uint64_t>* sinks2 = nullptr) const;
    
    // the shortest distance along any chain to each node after jumping from a given chain
    template<class BGraph, class XMerge>
    std::vector<std::vector<size_t>> post_switch_distances(const BGraph& graph, const XMerge& xmerge) const;
    
    template<class FinalFunc>
    std::vector<anchor_t> traceback_sparse_dp(const std::vector<match_set_t>& match_sets,
                                              const std::vector<std::vector<std::vector<dp_entry_t>>>& dp,
                                              const FinalFunc& final_function, double min_score,
                                              bool suppress_verbose_logging) const;
    
    inline double anchor_weight(size_t count1, size_t count2, size_t length) const;
    
    // affine edge score, assumes that both pairs of nodes are reachable
    template<class XMerge>
    double edge_weight(uint64_t from_id1, uint64_t to_id1, uint64_t from_id2, uint64_t to_id2,
                       double local_scale,
                       const XMerge& xmerge1, const XMerge& xmerge2,
                       const std::vector<std::vector<size_t>>& switch_dists1,
                       const std::vector<std::vector<size_t>>& switch_dists2) const;
    
public:
    
    // TODO: feels a bit inelegant exposing this here
    inline double anchor_weight(const anchor_t& anchor) const;
    
    template<class BGraph, class XMerge>
    void instrument_anchor_chain(const std::vector<anchor_t>& chain, double local_scale,
                                 const BGraph& graph1, const BGraph& graph2,
                                 const XMerge& xmerge1, const XMerge& xmerge2) const;
    
};






/*
 * Template implementations
 */


template<class BGraph, class XMerge>
std::vector<anchor_t> Anchorer::anchor_chain(std::vector<match_set_t>& matches,
                                             const BGraph& graph1,
                                             const BGraph& graph2,
                                             const XMerge& chain_merge1,
                                             const XMerge& chain_merge2,
                                             double anchor_scale) const {
    
    return anchor_chain(matches, graph1, graph2, chain_merge1, chain_merge2,
                        nullptr, nullptr, nullptr, nullptr, -1, false, chaining_algorithm, anchor_scale);
}


template<class BGraph, class XMerge>
std::vector<anchor_t> Anchorer::anchor_chain(std::vector<match_set_t>& matches,
                                             const BGraph& graph1,
                                             const BGraph& graph2,
                                             const XMerge& chain_merge1,
                                             const XMerge& chain_merge2,
                                             ChainAlgorithm chaining_algorithm_override) const {
    
    return anchor_chain(matches, graph1, graph2, chain_merge1, chain_merge2,
                        nullptr, nullptr, nullptr, nullptr,
                        -1, logging::level != logging::Debug, // this is for the calibration path, so we only log in the most verbose state
                        chaining_algorithm_override, global_scale);
}


template<class BGraph, class XMerge>
std::vector<anchor_t> Anchorer::global_anchor_chain(std::vector<match_set_t>& matches,
                                                    const BGraph& graph1,
                                                    const BGraph& graph2,
                                                    const XMerge& chain_merge1,
                                                    const XMerge& chain_merge2,
                                                    const std::vector<uint64_t>& sources1,
                                                    const std::vector<uint64_t>& sources2,
                                                    const std::vector<uint64_t>& sinks1,
                                                    const std::vector<uint64_t>& sinks2,
                                                    double anchor_scale,
                                                    size_t max_num_match_pairs_override) const {
    
    return anchor_chain(matches, graph1, graph2, chain_merge1, chain_merge2,
                        &sources1, &sources2, &sinks1, &sinks2,
                        max_num_match_pairs_override, true, chaining_algorithm, anchor_scale);
}

template<class BGraph, class XMerge>
std::vector<anchor_t> Anchorer::global_anchor_chain(std::vector<match_set_t>& matches,
                                                    const BGraph& graph1,
                                                    const BGraph& graph2,
                                                    const XMerge& chain_merge1,
                                                    const XMerge& chain_merge2,
                                                    const std::vector<uint64_t>& sources1,
                                                    const std::vector<uint64_t>& sources2,
                                                    const std::vector<uint64_t>& sinks1,
                                                    const std::vector<uint64_t>& sinks2,
                                                    ChainAlgorithm chaining_algorithm_override,
                                                    size_t max_num_match_pairs_override) const {
    
    return anchor_chain(matches, graph1, graph2, chain_merge1, chain_merge2,
                        &sources1, &sources2, &sinks1, &sinks2,
                        max_num_match_pairs_override, true, chaining_algorithm_override, global_scale);
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
                                             size_t max_num_match_pairs_override,
                                             bool suppress_verbose_logging,
                                             ChainAlgorithm chaining_algorithm_override,
                                             double anchor_scale) const {
    
    size_t max_num_match_pairs_local = max_num_match_pairs_override == -1 ? max_num_match_pairs : max_num_match_pairs_override;
    double local_scale = anchor_scale / global_scale;
    
    
    size_t total_num_pairs = 0;
    for (const auto& match_set : matches) {
        total_num_pairs += match_set.walks1.size() * match_set.walks2.size();
    }
    
    size_t removed = 0;
    if (total_num_pairs > max_num_match_pairs_local) {
        // we need to limit the number of nodes
        
        if (!suppress_verbose_logging) {
            logging::log(logging::Debug, "Selecting a maximum of " + std::to_string(max_num_match_pairs_local) + " matches for anchoring");
        }
        
        // prioritize based on the minimum count
        // TODO: is this a good criterion to use?
        std::stable_sort(matches.begin(), matches.end(), [](const match_set_t& a, const match_set_t& b) {
            //return std::min(a.count1, a.count2) < std::min(b.count1, b.count2);
            return a.count1 * a.count2 < b.count1 * b.count2;
        });
        
        // greedily choose matches as long as we have budget left
        size_t pairs_left = max_num_match_pairs_local;
        for (size_t i = 0; i < matches.size(); ++i) {
            auto& match = matches[i];
            size_t pair_count = match.walks1.size() * match.walks2.size();
            if (pairs_left >= pair_count) {
                pairs_left -= pair_count;
                std::swap(matches[i - removed], matches[i]);
            }
            else {
                // remove this match
                ++removed;
            }
        }
        
        if (debug_anchorer) {
            std::cerr << "moved " << removed << " unique anchor sequences to back limit to " << max_num_match_pairs_local << " total pairs\n";
        }
    }
    
    size_t num_match_sets = matches.size() - removed;
    
    // compute the optimal chain using DP
    std::vector<anchor_t> chain;
    switch (chaining_algorithm_override) {
        case SparseAffine:
            chain = std::move(sparse_affine_chain_dp(matches, graph1, graph2, chain_merge1, chain_merge2,
                                                     gap_open, gap_extend, local_scale, num_match_sets, suppress_verbose_logging,
                                                     sources1, sources2, sinks1, sinks2));
            break;
            
        case Sparse:
            chain = std::move(sparse_chain_dp(matches, graph1, chain_merge1, chain_merge2, num_match_sets,
                                              suppress_verbose_logging, sources1, sources2, sinks1, sinks2));
            break;
            
        case Exhaustive:
            chain = std::move(exhaustive_chain_dp(matches, graph1, graph2, chain_merge1, chain_merge2, false,
                                                  local_scale, num_match_sets, sources1, sources2, sinks1, sinks2));
            break;
            
        default:
            throw std::runtime_error("Unrecognized chaining algorithm: " + std::to_string((int) chaining_algorithm));
            break;
    }
    return chain;
}

inline double Anchorer::anchor_weight(const anchor_t& anchor) const {
    return anchor_weight(anchor.count1, anchor.count2, anchor.walk1.size());
}

inline double Anchorer::anchor_weight(size_t count1, size_t count2, size_t length) const {
    return (length_scale ? (double) length : 1.0) - pair_count_power * (log(count1 * count2 / count_penalty_threshold));
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
                                                    const std::vector<uint64_t>* sinks2) const {
            
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
        
        double weight = anchor_weight(match_set.count1, match_set.count2,
                                      match_set.walks1.front().size());
        
        for (size_t idx1 = 0; idx1 < match_set.walks1.size(); ++idx1) {
            for (size_t idx2 = 0; idx2 < match_set.walks2.size(); ++idx2) {
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
    }
    
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
                                                const std::vector<uint64_t>* sinks2) const {
    
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
                uint64_t chain2;
                size_t index;
                std::tie(chain2, index) = chain_merge2.chain(match_set.walks2[k].back());
                
                search_tree_data[chain2].emplace_back(key_t(index, i, j, k), mininf);
            }
        }
        
        // initialize the DP structure with a single-anchor chain at each position
        double weight = anchor_weight(match_set.count1, match_set.count2,
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
                double weight = anchor_weight(match_set.count1, match_set.count2,
                                              match_set.walks1.front().size());
                
                // we will consider all occurrences of this anchor in graph2
                for (size_t j = 0; j < match_set.walks2.size(); ++j) {
                    
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
    
    return traceback_sparse_dp(match_sets, dp, final_term, 0.0, suppress_verbose_logging);
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
                                                       const std::vector<uint64_t>* sinks2) const {
    
    assert((sources1 == nullptr) == (sources2 == nullptr));
    assert((sinks1 == nullptr) == (sinks2 == nullptr));
    
    if (!suppress_verbose_logging) {
        logging::log(logging::Debug, "Beginning affine-gap sparse chaining algorithm with local scale " + std::to_string(local_scale));
        
//        std::cerr << "modified score params:\n";
//        for (size_t i = 0; i < NumPW; ++i) {
//            std::cerr << local_scale * gap_open[i] << '\t' << local_scale * gap_extend[i] << '\n';
//        }
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
        double weight = anchor_weight(match_set.count1, match_set.count2,
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
                    else {
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
                    std::cerr << i << ' ' << j << ' ' << k << ' ' << std::get<0>(row[k]) << '\n';
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
                double weight = anchor_weight(match_set.count1, match_set.count2,
                                              match_set.walks1.front().size());
                
                // we will consider all occurrences of this anchor in graph2
                for (size_t k = 0; k < match_set.walks2.size(); ++k) {
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
        
    return traceback_sparse_dp(match_sets, dp, final_indel_score, min_score, suppress_verbose_logging);
}

template<class FinalFunc>
std::vector<anchor_t> Anchorer::traceback_sparse_dp(const std::vector<match_set_t>& match_sets,
                                                    const std::vector<std::vector<std::vector<dp_entry_t>>>& dp,
                                                    const FinalFunc& final_function, double min_score,
                                                    bool suppress_verbose_logging) const {
    
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
        
        // follow the backpointer from the DP structure
        here = dp[std::get<1>(here)][std::get<2>(here)][std::get<3>(here)];
    }
    
    // take out of reverse order
    std::reverse(anchors.begin(), anchors.end());
    
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
        std::cerr << '@' << '\t' << i << '\t' << chain[i].walk1.size() << '\t' << chain[i].count1 << '\t' << chain[i].count2 << '\t' << (chain[i].count1 * chain[i].count2) << '\t' << anchor_weight(chain[i].count1, chain[i].count2, chain[i].walk1.size()) << '\t';
        if (chaining_algorithm != SparseAffine || i == 0) {
            std::cerr << 0.0;
        }
        else {
            std::cerr << edge_weight(chain[i-1].walk1.back(), chain[i].walk1.front(),
                                     chain[i-1].walk2.back(), chain[i].walk2.front(),
                                     local_scale,
                                     xmerge1, xmerge2, switch_dists1, switch_dists2);
        }
        std::cerr << '\n';
    }
}


}

#endif /* centrolign_anchorer_hpp */
