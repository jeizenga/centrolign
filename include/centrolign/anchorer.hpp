#ifndef centrolign_anchorer_hpp
#define centrolign_anchorer_hpp

#include <vector>
#include <cstdint>
#include <iostream>

#include "centrolign/modify_graph.hpp"
#include "centrolign/gesa.hpp"

namespace centrolign {

struct anchor_t {
    anchor_t() = default;
    ~anchor_t() = default;
    std::vector<std::vector<uint64_t>> walks1;
    std::vector<std::vector<uint64_t>> walks2;
};

/*
 * Data structure finding anchors between two graphs
 */
class Anchorer {
public:
    Anchorer() = default;
    ~Anchorer() = default;
    
    // assumes that the graphs have already been given unique sentinels
    template<class BGraph>
    std::vector<anchor_t> find_anchors(const BGraph& graph1,
                                       const BGraph& graph2) const;
    
    
    /*
     * Configurable parameters
     */
    
    // the max count in either of the two graphs
    size_t max_count = 50;
    // the maximum number of occurrences of matches we will consider
    size_t max_num_match_pairs = 10000;
    
private:
    
};

/*
 * Template implementations
 */

template<class BGraph>
std::vector<anchor_t> Anchorer::find_anchors(const BGraph& graph1,
                                             const BGraph& graph2) const {
    
    std::vector<const BGraph*> graph_ptrs{&graph1, &graph2};
    
    GESA gesa(graph_ptrs);
    
    // records of (min count on either graph, total pairs, length, node)
    std::vector<std::tuple<size_t, size_t, size_t, GESANode>> matches;
    size_t total_num_pairs = 0;
    for (const auto& match : gesa.minimal_rare_matches(max_count)) {
        const auto& counts = gesa.component_counts(match.first);
        size_t num_pairs = counts[0] * counts[1];
        matches.emplace_back(std::min(counts[0], counts[1]),
                             num_pairs, match.second, match.first);
        total_num_pairs += num_pairs;
    }
    
    if (total_num_pairs > max_num_match_pairs) {
        // we need to limit the number of nodes
        
        // prioritize based on the minimum count
        // TODO: is this a good criterion to use?
        std::sort(matches.begin(), matches.end());
        
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
    }
    
    // walk out the matches into paths
    std::vector<anchor_t> anchors;
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

}

#endif /* centrolign_anchorer_hpp */
