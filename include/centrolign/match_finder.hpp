#ifndef centrolign_match_finder_hpp
#define centrolign_match_finder_hpp

#include <deque>
#include <limits>
#include <iostream>
#include <sstream>

#include "centrolign/graph.hpp"
#include "centrolign/modify_graph.hpp"

namespace centrolign {


// a set of walks of the same sequence in two graphs
struct match_set_t {
    match_set_t() = default;
    ~match_set_t() = default;
    std::vector<std::vector<uint64_t>> walks1;
    std::vector<std::vector<uint64_t>> walks2;
};

/*
 * Object that finds matches between two graphs
 */
class MatchFinder {
public:
    
    MatchFinder() = default;
    ~MatchFinder() = default;
    
    // returns minimal rare matches
    template<class BGraph>
    std::vector<match_set_t> find_matches(const BGraph& graph1, const BGraph& graph2) const;
    
    // returns minimal rare matches, with a back translation applied
    template<class BGraph>
    std::vector<match_set_t> find_matches(const BGraph& graph1, const BGraph& graph2,
                                          const std::vector<uint64_t>& back_translation1,
                                          const std::vector<uint64_t>& back_translation2) const;
    /*
     * Configurable parameters
     */
    
    // the max count in either of the two graphs
    size_t max_count = 50;
    // the maximum number of occurrences of matches we will consider
    size_t max_num_match_pairs = 10000;
    // the size at which to throw a PathGraphSizeError
    size_t size_limit = 32000000;
    
private:
     
    static const bool debug_match_finder = false;
    
    template<class BGraph>
    std::vector<match_set_t> find_matches(const BGraph& graph1, const BGraph& graph2,
                                          const std::vector<uint64_t>* back_translation1,
                                          const std::vector<uint64_t>* back_translation2) const;
};




/*
 * Template and inline implementations
 */

template<class BGraph>
std::vector<match_set_t> MatchFinder::find_matches(const BGraph& graph1, const BGraph& graph2) const {
    return find_matches(graph1, graph2, nullptr, nullptr);
}

template<class BGraph>
std::vector<match_set_t> MatchFinder::find_matches(const BGraph& graph1, const BGraph& graph2,
                                                   const std::vector<uint64_t>& back_translation1,
                                                   const std::vector<uint64_t>& back_translation2) const {
    return find_matches(graph1, graph2, &back_translation1, &back_translation2);
}

template<class BGraph>
std::vector<match_set_t> MatchFinder::find_matches(const BGraph& graph1, const BGraph& graph2,
                                                   const std::vector<uint64_t>* back_translation1,
                                                   const std::vector<uint64_t>* back_translation2) const {
    
    std::vector<const BGraph*> graph_ptrs{&graph1, &graph2};
    std::vector<const std::vector<uint64_t>*> trans_ptrs{back_translation1, back_translation2};
    
    GESA gesa(graph_ptrs, trans_ptrs, size_limit);
    
    logging::log(logging::Debug, "Finding minimal rare matches");
    
    // records of (min count on either graph, total pairs, length, node)
    std::vector<std::tuple<size_t, size_t, size_t, GESANode>> matches;
    size_t total_num_pairs = 0;
    for (const auto& match : gesa.minimal_rare_matches(max_count)) {
        
        const auto& counts = std::get<2>(match);
        
        if (debug_match_finder) {
            auto walked = gesa.walk_matches(std::get<0>(match), std::get<1>(match));
            const auto& walk_graph = walked.front().first == 0 ? graph1 : graph2;
            std::string seq;
            for (auto node_id : walked.front().second) {
                char base = walk_graph.label(node_id);
                if (base <= 4) {
                    base = decode_base(base);
                }
                seq.push_back(base);
            }
            std::cerr << "found match node " << std::get<0>(match).begin << ',' << std::get<0>(match).end << " with length " << std::get<1>(match) << ", counts " << counts[0] << " and " << counts[1] << " and sequence " << seq << '\n';
        }
        
        size_t num_pairs = counts[0] * counts[1];
        matches.emplace_back(std::min(counts[0], counts[1]),
                             num_pairs, std::get<1>(match), std::get<0>(match));
        total_num_pairs += num_pairs;
    }
    
    if (logging::level >= logging::Debug) {
        logging::log(logging::Debug, "Completed querying matches, found " + std::to_string(matches.size()) + " unique anchor sequences with max count " + std::to_string(max_count) + ", giving " + std::to_string(total_num_pairs) + " total anchor pairings");
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
        
        if (debug_match_finder) {
            std::cerr << "removed " << removed << " unique anchor sequences to limit to " << max_num_match_pairs << " total pairs\n";
        }
    }
    
    // walk out the matches into paths
    std::vector<match_set_t> match_sets;
    for (const auto& match : matches) {
        match_sets.emplace_back();
        auto& match_set = match_sets.back();
        for (auto& walked : gesa.walk_matches(std::get<3>(match), std::get<2>(match))) {
            // add them to the match set
            if (walked.first == 0) {
                match_set.walks1.emplace_back(std::move(walked.second));
            }
            else {
                match_set.walks2.emplace_back(std::move(walked.second));
            }
        }
    }
    
    return match_sets;
}

}

#endif /* centrolign_match_finder_hpp */
