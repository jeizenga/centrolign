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
    anchor_t() = default;
    ~anchor_t() = default;
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
    
    // compute a heaviest weight anchoring of a set of matches (consumes
    // the matches)
    template<class BGraph, class XMerge>
    std::vector<anchor_t> anchor_chain(std::vector<match_set_t>& matches,
                                       const BGraph& graph1,
                                       const BGraph& graph2,
                                       const XMerge& chain_merge1,
                                       const XMerge& chain_merge2) const;
    
    /*
     * Configurable parameters
     */
    
    // select which algorithm to use
    enum ChainAlgorithm { Exhaustive = 0, Sparse = 1, SparseAffine = 2 };
    ChainAlgorithm chaining_algorithm = Sparse;
    
    // power to raise the pair count to in the weight function
    double pair_count_power = 1.0;
    // anchor weight is proportional to length
    bool length_scale = false;
    // affine gap parameters
    // TODO: reframe count power in log-space, figure out the length indel i want to detect
    std::array<double, 2> gap_open{0.1, 10.0};
    std::array<double, 2> gap_extend{0.01, 0.0001};
    
protected:

    static const bool debug_anchorer;
    
    
    // a pair of two of the occurrence of a match
    struct AnchorNode {
        AnchorNode(size_t set, size_t idx1, size_t idx2, double weight) :
            set(set), idx1(idx1), idx2(idx2), weight(weight), in_degree(0) { }
        AnchorNode() = default;
        ~AnchorNode() = default;
        size_t set = 0;
        size_t idx1 = 0;
        size_t idx2 = 0;
        double weight = 0.0;
        size_t in_degree = 0; // for the sake of satisfying the topological_order interface
        std::vector<size_t> edges;
        std::vector<double> edge_weights;
    };
    
    // mutual reachability graph over anchor pairs
    class AnchorGraph {
    public:
        AnchorGraph() = default;
        ~AnchorGraph() = default;
        
        uint64_t add_node(size_t set, size_t idx1, size_t idx2, double weight);
        void add_edge(uint64_t from_id, uint64_t to_id, double weight = 0.0);
        
        std::vector<uint64_t> heaviest_weight_path() const;
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
    std::vector<anchor_t> exhaustive_chain_dp(std::vector<match_set_t>& match_sets,
                                              const BGraph& graph1,
                                              const BGraph& graph2,
                                              const XMerge& chain_merge1,
                                              const XMerge& chain_merge2,
                                              bool score_edges) const;
    
    template<class BGraph, class XMerge>
    std::vector<anchor_t> sparse_chain_dp(std::vector<match_set_t>& match_sets,
                                          const BGraph& graph1,
                                          const XMerge& chain_merge1,
                                          const XMerge& chain_merge2) const;
    
    template<size_t NumPW, class BGraph, class XMerge>
    std::vector<anchor_t> sparse_affine_chain_dp(std::vector<match_set_t>& match_sets,
                                                 const BGraph& graph1,
                                                 const BGraph& graph2,
                                                 const XMerge& xmerge1,
                                                 const XMerge& xmerge2,
                                                 const std::array<double, NumPW>& gap_open,
                                                 const std::array<double, NumPW>& gap_extend) const;
    
    // the shortest distance along any chain to each node after jumping from a given chain
    template<class BGraph, class XMerge>
    std::vector<std::vector<size_t>> post_switch_distances(const BGraph& graph, const XMerge& xmerge) const;
    
    std::vector<anchor_t> traceback_sparse_dp(std::vector<match_set_t>& match_sets,
                                              const std::vector<std::vector<std::vector<dp_entry_t>>>& dp) const;
    
    inline double anchor_weight(size_t count1, size_t count2, size_t length) const;
    
    // affine edge score, assumes that both pairs of nodes are reachable
    template<class XMerge>
    double edge_weight(uint64_t from_id1, uint64_t to_id1, uint64_t from_id2, uint64_t to_id2,
                       const XMerge& xmerge1, const XMerge& xmerge2,
                       const std::vector<std::vector<size_t>>& switch_dists1,
                       const std::vector<std::vector<size_t>>& switch_dists2) const;
    
};






/*
 * Template implementations
 */

template<class BGraph, class XMerge>
std::vector<anchor_t> Anchorer::anchor_chain(std::vector<match_set_t>& matches,
                                             const BGraph& graph1,
                                             const BGraph& graph2,
                                             const XMerge& chain_merge1,
                                             const XMerge& chain_merge2) const {
        
    // compute the optimal chain using DP
    std::vector<anchor_t> chain;
    switch (chaining_algorithm) {
        case SparseAffine:
            chain = std::move(sparse_affine_chain_dp(matches, graph1, graph2, chain_merge1, chain_merge2,
                                                     gap_open, gap_extend));
            break;
            
        case Sparse:
            chain = std::move(sparse_chain_dp(matches, graph1, chain_merge1, chain_merge2));
            break;
            
        case Exhaustive:
            chain = std::move(exhaustive_chain_dp(matches, graph1, graph2, chain_merge1, chain_merge2, false));
            break;
            
        default:
            throw std::runtime_error("Unrecognized chaining algorithm: " + std::to_string((int) chaining_algorithm));
            break;
    }
    return chain;
}

inline double Anchorer::anchor_weight(size_t count1, size_t count2, size_t length) const {
    double weight = 1.0 / pow(count1 * count2, pair_count_power);
    if (length_scale) {
        weight *= length;
    }
    return weight;
}

template <class BGraph, class XMerge>
std::vector<anchor_t> Anchorer::exhaustive_chain_dp(std::vector<match_set_t>& match_sets,
                                                    const BGraph& graph1,
                                                    const BGraph& graph2,
                                                    const XMerge& chain_merge1,
                                                    const XMerge& chain_merge2,
                                                    bool score_edges) const {
    
    std::vector<std::vector<size_t>> switch_dists1, switch_dists2;
    if (score_edges) {
        switch_dists1 = std::move(post_switch_distances(graph1, chain_merge1));
        switch_dists2 = std::move(post_switch_distances(graph2, chain_merge2));
    }
    
    // make a graph of match nodes
    AnchorGraph anchor_graph;
    for (size_t i = 0; i < match_sets.size(); ++i) {
        
        const auto& match_set = match_sets[i];
        
        double weight = anchor_weight(match_set.walks1.size(), match_set.walks2.size(),
                                      match_set.walks1.front().size());
        
        for (size_t idx1 = 0; idx1 < match_set.walks1.size(); ++idx1) {
            for (size_t idx2 = 0; idx2 < match_set.walks2.size(); ++idx2) {
                anchor_graph.add_node(i, idx1, idx2, weight);
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
                                                      chain_merge1, chain_merge2, switch_dists1, switch_dists2));
                }
                else {
                    anchor_graph.add_edge(node_id1, node_id2);
                }
            }
        }
    }
    
    // get heaviest path and convert into a chain
    std::vector<anchor_t> chain;
    for (auto node_id : anchor_graph.heaviest_weight_path()) {
        size_t set, idx1, idx2;
        std::tie(set, idx1, idx2) = anchor_graph.label(node_id);
        chain.emplace_back();
        auto& match_set = match_sets[set];
        auto& chain_node = chain.back();
        chain_node.walk1 = std::move(match_set.walks1[idx1]);
        chain_node.walk2 = std::move(match_set.walks2[idx2]);
        chain_node.count1 = match_set.walks1.size();
        chain_node.count2 = match_set.walks2.size();
    }
    return chain;
}

template<class BGraph, class XMerge>
std::vector<anchor_t> Anchorer::sparse_chain_dp(std::vector<match_set_t>& match_sets,
                                                const BGraph& graph1,
                                                const XMerge& chain_merge1,
                                                const XMerge& chain_merge2) const {
    
    logging::log(logging::Debug, "Beginning sparse chaining algorithm");
    
    // a key of (chain index, anchor set, walk1 index, walk2 index)
    using key_t = std::tuple<size_t, size_t, size_t, size_t>;
    
    // for each chain2, the initial search tree data for chains, records of (key_t, weight)
    std::vector<std::vector<std::pair<key_t, double>>> search_tree_data(chain_merge2.chain_size());
    
    // for each node, the list of (set, walk1) that start/end on it
    std::vector<std::vector<std::pair<size_t, size_t>>> starts(graph1.node_size()), ends(graph1.node_size());
    
    // for each set, for each walk1, for each walk2, the max weight and (set, walk1, walk2) for the traceback
    static const double mininf = std::numeric_limits<double>::lowest();
    std::vector<std::vector<std::vector<dp_entry_t>>> dp(match_sets.size());
    
    if (debug_anchorer) {
        std::cerr << "gathering anchor endpoint information\n";
    }
    
    // do the bookkeeping
    for (int64_t i = 0; i < match_sets.size(); ++i) {
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
        double weight = anchor_weight(match_set.walks1.size(), match_set.walks2.size(),
                                      match_set.walks1.front().size());
        dp[i].resize(match_set.walks1.size(),
                     std::vector<dp_entry_t>(match_set.walks2.size(),
                                             dp_entry_t(weight, -1, -1, -1)));
    }
    
    if (debug_anchorer) {
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
        if (logging::level >= logging::Debug && iter % 1000000 == 0) {
            logging::log(logging::Debug, "Entering iteration " + std::to_string(iter) + " of " + std::to_string(graph1.node_size()) + " in sparse chaining algorithm");
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
                std::cerr << "there is a forward edge to " << fwd_id << ", from node " << node_id << " on chain " << chain1 << '\n';
            }
            
            for (const auto& start : starts[fwd_id]) {
                
                // an anchor starts here in graph1
                
                if (debug_anchorer) {
                    std::cerr << "anchor " << start.first << ',' << start.second << " starts on " << fwd_id << '\n';
                }
                
                const auto& match_set = match_sets[start.first];
                auto& dp_row = dp[start.first][start.second];
                
                // the weight of this anchors in this set
                double weight = anchor_weight(match_set.walks1.size(), match_set.walks2.size(),
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
                            std::cerr << "looking for predecessor on chain " << chain2 << " with predecessor at index " << chain_preds2[chain2] << '\n';
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
    
    return traceback_sparse_dp(match_sets, dp);
}



template<class BGraph, class XMerge>
std::vector<std::vector<size_t>> Anchorer::post_switch_distances(const BGraph& graph, const XMerge& xmerge) const {
    
    std::vector<std::vector<size_t>> dists(graph.node_size(),
                                           std::vector<size_t>(xmerge.chain_size(), -1));
    
    for (auto node_id : topological_order(graph)) {
        auto& row = dists[node_id];
        for (uint64_t p = 0; p < xmerge.chain_size(); ++p) {
            if (xmerge.index_on(node_id, p) != -1) {
                row[p] = 0;
            }
            else {
                const auto& preds = xmerge.predecessor_indexes(node_id);
                for (auto prev_id : graph.previous(node_id)) {
                    if (xmerge.predecessor_indexes(prev_id)[p] == preds[p]) {
                        row[p] = std::min(row[p], dists[prev_id][p] + graph.label_size(prev_id));
                    }
                }
            }
        }
    }
    
    return dists;
}

template<size_t NumPW, class BGraph, class XMerge>
std::vector<anchor_t> Anchorer::sparse_affine_chain_dp(std::vector<match_set_t>& match_sets,
                                                       const BGraph& graph1,
                                                       const BGraph& graph2,
                                                       const XMerge& xmerge1,
                                                       const XMerge& xmerge2,
                                                       const std::array<double, NumPW>& gap_open,
                                                       const std::array<double, NumPW>& gap_extend) const {
    
    
    logging::log(logging::Debug, "Beginning sparse affine chaining algorithm");
    
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
    auto source_shift = [&](size_t i, size_t j, size_t k, uint64_t path1, uint64_t path2) -> int64_t {
        const auto& match_set = match_sets[i];
        return xmerge1.index_on(match_set.walks1[j].back(), path1) - xmerge2.index_on(match_set.walks2[k].back(), path2);
    };
    // add the source coordinates to grab as backpointers
    auto get_key = [&](size_t i, size_t j, size_t k, uint64_t path1, uint64_t path2) -> key_t {
        return key_t(source_shift(i, j, k, path1, path2), i, j, k);
    };
    // the offset on path from graph2
    auto get_key_offset = [&](size_t i, size_t k, uint64_t path2) -> size_t {
        return xmerge2.index_on(match_sets[i].walks2[k].back(), path2);
    };
    // the "effective offset" on path of graph2 (true offsets before this can reach it)
    auto get_query_offset = [&](size_t i, size_t k, uint64_t path2) -> size_t {
        return xmerge2.predecessor_indexes(match_sets[i].walks2[k].front())[path2] + 1;
    };
    
    // the gap contribution from the destination anchor
    auto query_shift = [&](size_t i, size_t j, size_t k, uint64_t path1, uint64_t path2) {
        const auto& match_set = match_sets[i];
        uint64_t n1 = match_set.walks1[j].front();
        uint64_t n2 = match_set.walks2[k].front();
        return (xmerge1.predecessor_indexes(n1)[path1] - xmerge2.predecessor_indexes(n2)[path2]
                + switch_dists1[n1][path1] - switch_dists2[n2][path2]);
    };
    
    // for each node, the list of (set, walk1) that start/end on it
    std::vector<std::vector<std::pair<size_t, size_t>>> starts(graph1.node_size()), ends(graph1.node_size());
    
    // for each set, for each walk1, for each walk2, the max weight and (set, walk1, walk2) for the traceback
    static const double mininf = std::numeric_limits<double>::lowest();
    std::vector<std::vector<std::vector<dp_entry_t>>> dp(match_sets.size());
        
    // for each path1, for each path2, list of (search index, value)
    std::vector<std::vector<std::vector<std::tuple<key_t, size_t, double>>>> search_tree_data;
    search_tree_data.resize(xmerge1.chain_size(), std::vector<std::vector<std::tuple<key_t, size_t, double>>>(xmerge2.chain_size()));
    
    if (debug_anchorer) {
        std::cerr << "doing bookkeeping for sparse affine algorithm\n";
    }
    
    // do the bookkeeping
    for (int64_t i = 0; i < match_sets.size(); ++i) {
        // get the starts and ends of anchors on graph 1
        auto& match_set = match_sets[i];
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
        double weight = anchor_weight(match_set.walks1.size(), match_set.walks2.size(),
                                      match_set.walks1.front().size());
        dp[i].resize(match_set.walks1.size(),
                     std::vector<dp_entry_t>(match_set.walks2.size(),
                                             dp_entry_t(weight, -1, -1, -1)));
    }
    
    if (debug_anchorer) {
        std::cerr << "constructing search trees\n";
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
    
    size_t iter = 0;
    for (uint64_t node_id : topological_order(graph1)) {
        
        ++iter;
        if (logging::level >= logging::Debug && iter % 100000 == 0) {
            logging::log(logging::Debug, "Entering iteration " + std::to_string(iter) + " of " + std::to_string(graph1.node_size()) + " in sparse chaining algorithm");
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
                                value = dp_val + gap_extend[pw / 2] * source_shift(end.first, end.second, k, p1, p2);
                            }
                            else {
                                // d1 < d2
                                value = dp_val - gap_extend[pw / 2 - 1] * source_shift(end.first, end.second, k, p1, p2);
                            }
                            auto& tree = search_trees[pw][p1][p2];
                            auto it = tree.find(key1, key2);
                            if (value > std::get<2>(*it)) {
                                if (debug_anchorer) {
                                    std::cerr << "found new opt " << value << " for key " << std::get<0>(key1) << ',' << std::get<1>(key1) << ',' << std::get<2>(key1) << ',' << std::get<3>(key1) << ", offset " << key2 << " in piecewise component " << pw << "\n";
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
                double weight = anchor_weight(match_set.walks1.size(), match_set.walks2.size(),
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
                            std::cerr << "query shift is " << query << " on chain " << chain2 << " and offset is " << offset << '\n';
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
                                    double value = std::get<2>(*it) + weight - gap_open[pw / 2] - gap_extend[pw / 2] * query;
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
                                    double value = std::get<2>(*it) + weight - gap_open[pw / 2 - 1] + gap_extend[pw / 2 - 1] * query;
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
    
    return traceback_sparse_dp(match_sets, dp);
}

template<class XMerge>
double Anchorer::edge_weight(uint64_t from_id1, uint64_t to_id1, uint64_t from_id2, uint64_t to_id2,
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
                    weight = std::max(weight, -gap_open[i] - gap_extend[i] * gap);
                }
            }
        }
    }
    return weight;
}


}

#endif /* centrolign_anchorer_hpp */
