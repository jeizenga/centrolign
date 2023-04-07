#ifndef centrolign_anchorer_hpp
#define centrolign_anchorer_hpp

#include <vector>
#include <cstdint>
#include <cmath>
#include <iostream>

#include "centrolign/chain_merge.hpp"
#include "centrolign/modify_graph.hpp"
#include "centrolign/gesa.hpp"

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
                                       const ChainMerge& chain_merge2);
    
    /*
     * Configurable parameters
     */
    
    // the max count in either of the two graphs
    size_t max_count = 50;
    // the maximum number of occurrences of matches we will consider
    size_t max_num_match_pairs = 10000;
    // anchor weight is proportional to square root of occurrences
    bool root_scale = false;
    // anchor weight is proportional to length
    bool length_scale = false;
    
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
        std::tuple<size_t, size_t, size_t> anchor(uint64_t node_id) const;
        
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
    
};

/*
 * Template implementations
 */

template<class BGraph>
std::vector<Anchorer::anchor_set_t> Anchorer::find_matches(const BGraph& graph1,
                                                           const BGraph& graph2) const {
    
    std::vector<const BGraph*> graph_ptrs{&graph1, &graph2};
    
    GESA gesa(graph_ptrs);
    
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
                char base = walk_graph.base(node_id);
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


template<class BGraph>
std::vector<anchor_t> Anchorer::anchor_chain(const BGraph& graph1,
                                             const BGraph& graph2,
                                             const ChainMerge& chain_merge1,
                                             const ChainMerge& chain_merge2) {
    
    // get the matches
    std::vector<anchor_set_t> anchor_sets = find_matches(graph1, graph2);
    
    // make a graph of match nodes
    AnchorGraph anchor_graph;
    for (size_t i = 0; i < anchor_sets.size(); ++i) {
        
        const auto& anchor_set = anchor_sets[i];
        
        double weight = 1.0 / double(anchor_set.walks1.size() * anchor_set.walks2.size());
        if (root_scale) {
            weight = sqrt(weight);
        }
        if (length_scale) {
            weight *= anchor_set.walks1.front().size();
        }
        
        for (size_t idx1 = 0; idx1 < anchor_set.walks1.size(); ++idx1) {
            for (size_t idx2 = 0; idx2 < anchor_set.walks2.size(); ++idx2) {
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
            std::tie(set1, idx11, idx21) = anchor_graph.anchor(node_id1);
            std::tie(set2, idx12, idx22) = anchor_graph.anchor(node_id2);
            
            auto& anchor_set1 = anchor_sets[set1];
            auto& anchor_set2 = anchor_sets[set2];
            
            if (chain_merge1.reachable(anchor_set1.walks1[idx11].back(),
                                       anchor_set2.walks1[idx12].front()) &&
                chain_merge2.reachable(anchor_set1.walks2[idx21].back(),
                                       anchor_set2.walks2[idx22].front())) {
                
                anchor_graph.add_edge(node_id1, node_id2);
            }
        }
    }
    
    // get heaviest path and convert into a chain
    std::vector<anchor_t> chain;
    for (auto node_id : anchor_graph.heaviest_weight_path()) {
        size_t set, idx1, idx2;
        std::tie(set, idx1, idx2) = anchor_graph.anchor(node_id);
        chain.emplace_back();
        auto& anchor_set = anchor_sets[set];
        auto& chain_node = chain.back();
        chain_node.walk1 = std::move(anchor_set.walks1[idx1]);
        chain_node.walk2 = std::move(anchor_set.walks2[idx2]);
        chain_node.count1 = anchor_set.walks1.size();
        chain_node.count2 = anchor_set.walks2.size();
    }
    return chain;
}

}

#endif /* centrolign_anchorer_hpp */
