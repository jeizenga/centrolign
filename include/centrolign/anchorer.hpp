#ifndef centrolign_anchorer_hpp
#define centrolign_anchorer_hpp

#include <vector>
#include <cstdint>
#include <cmath>
#include <iostream>
#include <limits>

#include "centrolign/chain_merge.hpp"
#include "centrolign/modify_graph.hpp"
#include "centrolign/gesa.hpp"
#include "centrolign/max_search_tree.hpp"
#include "centrolign/topological_order.hpp"

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
    
    // compute a heaviest weight anchoring of
    template<class BGraph>
    std::vector<anchor_t> anchor_chain(const BGraph& graph1,
                                       const BGraph& graph2,
                                       const ChainMerge& chain_merge1,
                                       const ChainMerge& chain_merge2) const;
    
    /*
     * Configurable parameters
     */
    
    // the max count in either of the two graphs
    size_t max_count = 50;
    // the maximum number of occurrences of matches we will consider
    size_t max_num_match_pairs = 10000;
    // power to raise the pair count to in the weight function
    double pair_count_power = 1.0;
    // anchor weight is proportional to length
    bool length_scale = false;
    // use the sparse chaining algorithm
    bool sparse_chaining = false;
    
protected:

    static const bool debug_anchorer;
    static const bool basic_logging;
    
    // a set of walks of the same sequence in two graphs
    struct anchor_set_t {
        anchor_set_t() = default;
        ~anchor_set_t() = default;
        std::vector<std::vector<uint64_t>> walks1;
        std::vector<std::vector<uint64_t>> walks2;
    };
    
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
    };
    
    // mutual reachability graph over anchor pairs
    class AnchorGraph {
    public:
        AnchorGraph() = default;
        ~AnchorGraph() = default;
        
        uint64_t add_node(size_t set, size_t idx1, size_t idx2, double weight);
        void add_edge(uint64_t from_id, uint64_t to_id);
        
        std::vector<uint64_t> heaviest_weight_path() const;
        // get (set, idx1, idx2)
        std::tuple<size_t, size_t, size_t> label(uint64_t node_id) const;
        
        size_t node_size() const;
        const std::vector<size_t> next(uint64_t node_id) const;
        size_t next_size(uint64_t node_id) const;
        size_t previous_size(uint64_t node_id) const;
        
    private:
        
        std::vector<AnchorNode> nodes;
    };
    
    // assumes that the graphs have already been given unique sentinels
    // note: these will never show up in anchors because they can't match
    template<class BGraph>
    std::vector<anchor_set_t> find_matches(const BGraph& graph1,
                                           const BGraph& graph2) const;
    
    std::vector<anchor_t> exhaustive_chain_dp(std::vector<anchor_set_t>& anchor_sets,
                                              const ChainMerge& chain_merge1,
                                              const ChainMerge& chain_merge2) const;
    
    template<class BGraph>
    std::vector<anchor_t> sparse_chain_dp(std::vector<anchor_set_t>& anchor_sets,
                                          const BGraph& graph1,
                                          const ChainMerge& chain_merge1,
                                          const ChainMerge& chain_merge2) const;
    
    inline double anchor_weight(size_t count1, size_t count2, size_t length) const;
    
};

/*
 * Template implementations
 */


template<class BGraph>
std::vector<anchor_t> Anchorer::anchor_chain(const BGraph& graph1,
                                             const BGraph& graph2,
                                             const ChainMerge& chain_merge1,
                                             const ChainMerge& chain_merge2) const {
    
    // get the matches
    std::vector<anchor_set_t> anchor_sets = find_matches(graph1, graph2);
    
    // compute the optimal chain using DP
    std::vector<anchor_t> chain;
    if (sparse_chaining) {
        chain = sparse_chain_dp(anchor_sets, graph1, chain_merge1, chain_merge2);
    }
    else {
        chain = exhaustive_chain_dp(anchor_sets, chain_merge1, chain_merge2);
    }
    return chain;
}

template<class BGraph>
std::vector<Anchorer::anchor_set_t> Anchorer::find_matches(const BGraph& graph1,
                                                           const BGraph& graph2) const {
    
    std::vector<const BGraph*> graph_ptrs{&graph1, &graph2};
    
    GESA gesa(graph_ptrs);
    
    if (basic_logging) {
        std::cerr << "finding minimal rare matches\n";
    }
    
    // records of (min count on either graph, total pairs, length, node)
    std::vector<std::tuple<size_t, size_t, size_t, GESANode>> matches;
    size_t total_num_pairs = 0;
    for (const auto& match : gesa.minimal_rare_matches(max_count)) {
        
        const auto& counts = gesa.component_counts(match.first);
        
        if (debug_anchorer) {
            auto walked = gesa.walk_matches(match.first, match.second);
            const auto& walk_graph = walked.front().first == 0 ? graph1 : graph2;
            std::string seq;
            for (auto node_id : walked.front().second) {
                char base = walk_graph.label(node_id);
                if (base <= 4) {
                    base = decode_base(base);
                }
                seq.push_back(base);
            }
            std::cerr << "found match node " << match.first.begin << ',' << match.first.end << " with length " << match.second << ", counts " << counts[0] << " and " << counts[1] << " and sequence " << seq << '\n';
        }
        
        size_t num_pairs = counts[0] * counts[1];
        matches.emplace_back(std::min(counts[0], counts[1]),
                             num_pairs, match.second, match.first);
        total_num_pairs += num_pairs;
    }
    
    if (basic_logging || debug_anchorer) {
        std::cerr << "completed querying matches, found " << matches.size() << " unique anchor sequences with max count " << max_count << ", giving " << total_num_pairs << " total anchor pairings\n";
    }
    
    if (total_num_pairs > max_num_match_pairs) {
        // we need to limit the number of nodes
        
        // prioritize based on the minimum count
        // TODO: is this a good criterion to use?
        std::stable_sort(matches.begin(), matches.end());
        
        // greedily choose matches as long as we have budget left
        size_t removed = 0;
        size_t pairs_left = max_num_match_pairs;
        for (size_t i = 0; i < matches.size(); ++i) {
            auto& match = matches[i];
            if (pairs_left >= std::get<1>(match)) {
                pairs_left -= std::get<1>(match);
                matches[i - removed] = std::move(match);
            }
            else {
                ++removed;
            }
        }
        matches.resize(matches.size() - removed);
        
        if (basic_logging || debug_anchorer) {
            std::cerr << "removed " << removed << " unique anchor sequences to limit to " << max_num_match_pairs << " total pairs\n";
        }
    }
    
    // walk out the matches into paths
    std::vector<anchor_set_t> anchors;
    for (const auto& match : matches) {
        anchors.emplace_back();
        auto& anchor = anchors.back();
        for (auto& walked : gesa.walk_matches(std::get<3>(match), std::get<2>(match))) {
            // add
            if (walked.first == 0) {
                anchor.walks1.emplace_back(std::move(walked.second));
            }
            else {
                anchor.walks2.emplace_back(std::move(walked.second));
            }
        }
    }
    
    return anchors;
}

inline double Anchorer::anchor_weight(size_t count1, size_t count2, size_t length) const {
    double weight = 1.0 / pow(count1 * count2, pair_count_power);
    if (length_scale) {
        weight *= length;
    }
    return weight;
}



template<class BGraph>
std::vector<anchor_t> Anchorer::sparse_chain_dp(std::vector<anchor_set_t>& anchor_sets,
                                                const BGraph& graph1,
                                                const ChainMerge& chain_merge1,
                                                const ChainMerge& chain_merge2) const {
    
    if (debug_anchorer || basic_logging) {
        std::cerr << "beginning sparse chaining algorithm\n";
    }
    
    // a key of (chain index, anchor set, walk1 index, walk2 index)
    using key_t = std::tuple<size_t, size_t, size_t, size_t>;
    // a dynamic programming value and backpointer of (anchor set, walk1 index, walk2 index)
    using dp_entry_t = std::tuple<double, size_t, size_t, size_t>;
    
    
    // for each chain2, the initial search tree data for chains, records of (key_t, weight)
    std::vector<std::vector<std::pair<key_t, double>>> search_tree_data(chain_merge2.chain_size());
    
    // for each node, the list of (set, walk1) that start/end on it
    std::vector<std::vector<std::pair<size_t, size_t>>> starts(graph1.node_size()), ends(graph1.node_size());
    
    // for each set, for each walk1, for each walk2, the max weight and (set, walk1, walk2) for the traceback
    static const double mininf = std::numeric_limits<double>::lowest();
    std::vector<std::vector<std::vector<dp_entry_t>>> dp(anchor_sets.size());
    
    if (debug_anchorer) {
        std::cerr << "gathering anchor endpoint information\n";
    }
    
    // do the bookkeeping
    for (int64_t i = 0; i < anchor_sets.size(); ++i) {
        // get the starts and ends of anchors on graph 1
        auto& anchor_set = anchor_sets[i];
        for (size_t j = 0; j < anchor_set.walks1.size(); ++j) {
            // get the starts and ends of anchor on graph 1
            auto& walk = anchor_set.walks1[j];
            starts[walk.front()].emplace_back(i, j);
            ends[walk.back()].emplace_back(i, j);
            
            // get the chain index of anchor ends on graph 2
            for (size_t k = 0; k < anchor_set.walks2.size(); ++k) {
                uint64_t chain2 = chain_merge2.chain(anchor_set.walks2[k].back());
                size_t index = chain_merge2.chain_index(anchor_set.walks2[k].back());
                
                search_tree_data[chain2].emplace_back(key_t(index, i, j, k), mininf);
            }
        }
        
        // initialize the DP structure with a single-anchor chain at each position
        double weight = anchor_weight(anchor_set.walks1.size(), anchor_set.walks2.size(),
                                      anchor_set.walks1.front().size());
        dp[i].resize(anchor_set.walks1.size(),
                     std::vector<dp_entry_t>(anchor_set.walks2.size(),
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
    
    // for each chain1, for each chain2, a tree over seed end chain indexes (and a null start)
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
    
//    auto dump_tree = [&](const MaxSearchTree<key_t, double>& tree) {
//        for (auto key_val_pair : tree) {
//            std::cerr << std::get<0>(key_val_pair.first) << ' ' << std::get<1>(key_val_pair.first) << ' ' << std::get<2>(key_val_pair.first) << ' ' << std::get<3>(key_val_pair.first) << " -> " << key_val_pair.second << '\n';
//        }
//    };
    
    if (debug_anchorer) {
        std::cerr << "computing forward edges\n";
    }
    
    // get the edges to chain neighbors
    auto forward_edges = chain_merge1.chain_forward_edges();
    
    
    if (debug_anchorer) {
        std::cerr << "beginning main DP iteration\n";
    }
    
    size_t iter = 0;
    for (uint64_t node_id : topological_order(graph1)) {
        
        ++iter;
        if ((basic_logging || debug_anchorer) && iter % 1000000 == 0) {
            std::cerr << "entering iteration " << iter << " of " << graph1.node_size() << " in sparse chaining algorithm\n";
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
            
            auto& anchor_set = anchor_sets[end.first];
            uint64_t chain1 = chain_merge1.chain(node_id);
            auto& dp_row = dp[end.first][end.second];
            
            // add a value in the appropriate tree for each occurrences in graph2
            for (size_t i = 0; i < anchor_set.walks2.size(); ++i) {
                const auto& walk2 = anchor_set.walks2[i];
                uint64_t chain2 = chain_merge2.chain(walk2.back());
                size_t index = chain_merge2.chain_index(walk2.back());
                
                auto& tree = search_trees[chain1][chain2];
                // TODO: do we really need equal_range if we make everything unique
                auto it = tree.equal_range(key_t(index, end.first, end.second, i)).first;
                if (it->second < std::get<0>(dp_row[i])) {
                    tree.update(it, std::get<0>(dp_row[i]));
                    if (debug_anchorer) {
                        std::cerr << "recording increased DP value of " << it->second << " on chain " << chain2 << " with key " << std::get<0>(it->first) << ' ' << std::get<1>(it->first) << ' ' << std::get<2>(it->first) << ' ' << std::get<3>(it->first) << '\n';
                    }
                }
//                if (debug_anchorer) {
//                    std::cerr << "tree state of " << chain1 << ' ' << chain2 << ":\n";
//                    dump_tree(tree);
//                }
            }
        }
        
        if (debug_anchorer) {
            std::cerr << "looking for chain forward edges\n";
        }
        
        uint64_t chain1 = chain_merge1.chain(node_id);
        
        // carry the current updates to any anchor starts for which this is the
        // the last node to reach from this chain (then all DP values from anchors
        // that end in this chain on graph1 have been completed)
        for (uint64_t fwd_id : forward_edges[node_id]) {
            
            if (debug_anchorer) {
                std::cerr << "there is a forward edge to " << fwd_id << ", from node " << node_id << " on chain " << chain1 << '\n';
            }
            
            for (const auto& start : starts[fwd_id]) {
                
                // an anchor starts here in graph1
                
                if (debug_anchorer) {
                    std::cerr << "anchor " << start.first << ',' << start.second << " starts on " << fwd_id << '\n';
                }
                
                const auto& anchor_set = anchor_sets[start.first];
                auto& dp_row = dp[start.first][start.second];
                
                // the weight of this anchors in this set
                double weight = anchor_weight(anchor_set.walks1.size(), anchor_set.walks2.size(),
                                              anchor_set.walks1.front().size());
                
                // we will consider all occurrences of this anchor in graph2
                for (size_t j = 0; j < anchor_set.walks2.size(); ++j) {
                    
                    auto& dp_entry = dp_row[j];
                    
                    if (debug_anchorer) {
                        std::cerr << "walk2 index " << j << " of " << anchor_set.walks2.size() << " starts on node " << anchor_set.walks2[j].front() << " and has current DP entry:\n";
                        std::cerr << std::get<0>(dp_entry) << ' ' << std::get<1>(dp_entry) << ' ' << std::get<2>(dp_entry) << ' ' << std::get<3>(dp_entry) << '\n';
                    }
                    // we will check for previous DP values that can reach this one in graph2 from
                    // each of the chains in the chain partition of graph2
                    const auto& chain_preds2 = chain_merge2.predecessor_indexes(anchor_set.walks2[j].front());
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
//                        if (debug_anchorer) {
//                            std::cerr << "tree state of " << chain1 << ' ' << chain2 << ":\n";
//                            dump_tree(tree);
//                        }
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
    
    if (debug_anchorer) {
        std::cerr << "finding optimum\n";
    }
    
    // find the optimum dynamic programming values
    dp_entry_t here(mininf, -1, -1, -1);
    for (size_t set = 0; set < anchor_sets.size(); ++set) {
        const auto& set_dp = dp[set];
        for (size_t i = 0; i < set_dp.size(); ++i) {
            const auto& dp_row = set_dp[i];
            for (size_t j = 0; j < dp_row.size(); ++j) {
                if (std::get<0>(dp_row[j]) > std::get<0>(here)) {
                    here = dp_entry_t(std::get<0>(dp_row[j]), set, i, j);
                }
            }
        }
    }
    
    if (debug_anchorer) {
        std::cerr << "doing traceback\n";
    }
    
    // traceback into a chain
    std::vector<anchor_t> anchors;
    while (std::get<1>(here) != -1) {
        
        if (debug_anchorer) {
            std::cerr << "following traceback to set " << std::get<1>(here) << ", walk pair " << std::get<2>(here) << " " << std::get<3>(here) << '\n';
        }
        
        // grab the anchors that we used from their set
        auto& anchor_set = anchor_sets[std::get<1>(here)];
        anchors.emplace_back();
        auto& anchor = anchors.back();
        anchor.walk1 = std::move(anchor_set.walks1[std::get<2>(here)]);
        anchor.count1 = anchor_set.walks1.size();
        anchor.walk2 = std::move(anchor_set.walks2[std::get<3>(here)]);
        anchor.count2 = anchor_set.walks2.size();
        
        // follow the backpointer from the DP structure
        here = dp[std::get<1>(here)][std::get<2>(here)][std::get<3>(here)];
    }
    
    // take out of reverse order
    std::reverse(anchors.begin(), anchors.end());
    
    if (debug_anchorer) {
        std::cerr << "completed sparse chaining\n";
    }
    
    return anchors;
}

}

#endif /* centrolign_anchorer_hpp */
