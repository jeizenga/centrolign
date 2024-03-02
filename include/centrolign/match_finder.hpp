#ifndef centrolign_match_finder_hpp
#define centrolign_match_finder_hpp

#include <deque>
#include <limits>
#include <iostream>
#include <sstream>

#include "centrolign/graph.hpp"
#include "centrolign/modify_graph.hpp"
#include "centrolign/path_esa.hpp"
#include "centrolign/score_function.hpp"

namespace centrolign {


// a set of walks of the same sequence in two graphs
struct match_set_t {
    match_set_t() noexcept = default;
    match_set_t(const match_set_t& other) noexcept = default;
    match_set_t(match_set_t&& other) noexcept = default;
    ~match_set_t() = default;
    match_set_t& operator=(const match_set_t& other) noexcept = default;
    match_set_t& operator=(match_set_t&& other) noexcept = default;
    
    std::vector<std::vector<uint64_t>> walks1;
    std::vector<std::vector<uint64_t>> walks2;
    size_t count1 = 0;
    size_t count2 = 0;
};

/*
 * Object that finds matches between two graphs
 */
class MatchFinder {
public:
    
    MatchFinder(const ScoreFunction& score_function) : score_function(&score_function) {}
    MatchFinder() = default;
    ~MatchFinder() = default;
    
    // returns minimal rare matches
    template<class BGraph>
    std::vector<match_set_t> find_matches(const BGraph& graph1, const BGraph& graph2,
                                          const SentinelTableau& tableau1, const SentinelTableau& tableau2) const;
    
    // returns minimal rare matches, with a back translation applied
    template<class BGraph>
    std::vector<match_set_t> find_matches(const BGraph& graph1, const BGraph& graph2,
                                          const SentinelTableau& tableau1, const SentinelTableau& tableau2,
                                          const std::vector<uint64_t>& back_translation1,
                                          const std::vector<uint64_t>& back_translation2) const;
    /*
     * Configurable parameters
     */
    
    // the max count in either of the two graphs
    size_t max_count = 50;
    // throw a GESASizeError if it grows to this times as much as the graph size
    size_t size_limit_factor = 16;
    // query against only embedded paths instead of arbitrary recombinants
    bool path_matches = false;
private:
     
    static const bool debug_match_finder = false;
    
    const ScoreFunction* const score_function = nullptr;
    
    template<class BGraph>
    std::vector<match_set_t> find_matches(const BGraph& graph1, const BGraph& graph2,
                                          const SentinelTableau& tableau1, const SentinelTableau& tableau2,
                                          const std::vector<uint64_t>* back_translation1,
                                          const std::vector<uint64_t>* back_translation2) const;
    
    template<class Index>
    std::vector<match_set_t> query_index(const Index& index) const;
};




/*
 * Template and inline implementations
 */

template<class BGraph>
std::vector<match_set_t> MatchFinder::find_matches(const BGraph& graph1, const BGraph& graph2,
                                                   const SentinelTableau& tableau1, const SentinelTableau& tableau2) const {
    return find_matches(graph1, graph2, tableau1, tableau2, nullptr, nullptr);
}

template<class BGraph>
std::vector<match_set_t> MatchFinder::find_matches(const BGraph& graph1, const BGraph& graph2,
                                                   const SentinelTableau& tableau1, const SentinelTableau& tableau2,
                                                   const std::vector<uint64_t>& back_translation1,
                                                   const std::vector<uint64_t>& back_translation2) const {
    return find_matches(graph1, graph2, tableau1, tableau2, &back_translation1, &back_translation2);
}

template<class BGraph>
std::vector<match_set_t> MatchFinder::find_matches(const BGraph& graph1, const BGraph& graph2,
                                                   const SentinelTableau& tableau1, const SentinelTableau& tableau2,
                                                   const std::vector<uint64_t>* back_translation1,
                                                   const std::vector<uint64_t>* back_translation2) const {
    
    std::vector<const BGraph*> graph_ptrs{&graph1, &graph2};
    std::vector<const std::vector<uint64_t>*> trans_ptrs{back_translation1, back_translation2};
    std::vector<const SentinelTableau*> tableau_ptrs{&tableau1, &tableau2};
    
    // FIXME: this doesn't play well with the graphs being simplified recursively
    size_t size_limit = size_limit_factor * (graph1.node_size() + graph2.node_size());
    
    std::vector<match_set_t> matches;
    
    if (path_matches) {
        if (back_translation1 || back_translation2) {
            throw std::runtime_error("Path restricted match queries have not been implemented with back translations");
        }
        PathESA path_esa(graph_ptrs, tableau_ptrs);
        matches = move(query_index(path_esa));
    }
    else {
        GESA gesa(graph_ptrs, trans_ptrs, size_limit);
        matches = move(query_index(gesa));
    }
    
    return matches;
}

template<class Index>
std::vector<match_set_t> MatchFinder::query_index(const Index& index) const {
    
    logging::log(logging::Debug, "Finding minimal rare matches");
    
    // records of (min count on either graph, total pairs, length, node)
    std::vector<std::tuple<size_t, size_t, size_t, SANode>> matches;
    size_t kept_num_pairs = 0;
    size_t removed_num_pairs = 0;
    for (const auto& match : index.minimal_rare_matches(max_count)) {
        
        const auto& counts = std::get<2>(match);
        
        if (debug_match_finder) {
            auto walked = index.walk_matches(std::get<0>(match), std::get<1>(match));
            std::cerr << "found match node " << std::get<0>(match).begin << ',' << std::get<0>(match).end << " with length " << std::get<1>(match) << ", counts " << counts[0] << " and " << counts[1] << '\n';
        }
        
        // only keep the match if it has positive score
        size_t num_pairs = counts[0] * counts[1];
        if (score_function->anchor_weight(counts[0], counts[1], std::get<1>(match)) > 0.0) {
            matches.emplace_back(std::min(counts[0], counts[1]),
                                 num_pairs, std::get<1>(match), std::get<0>(match));
            kept_num_pairs += num_pairs;
        }
        else {
            removed_num_pairs += num_pairs;
        }
    }
    
    if (logging::level >= logging::Debug) {
        logging::log(logging::Debug, "Completed querying matches, found " + std::to_string(matches.size()) + " unique anchor sequences with max count " + std::to_string(max_count) + ", giving " + std::to_string(kept_num_pairs) + " anchor pairings with positive score and " + std::to_string(removed_num_pairs) + " pairs that were removed.");
    }
    
    logging::log(logging::Debug, "Walking out paths of match sequences");
        
    // walk out the matches into paths
    std::vector<match_set_t> match_sets;
    for (const auto& match : matches) {
        match_sets.emplace_back();
        auto& match_set = match_sets.back();
        for (auto& walked : index.walk_matches(std::get<3>(match), std::get<2>(match))) {
            // add them to the match set
            if (walked.first == 0) {
                match_set.walks1.emplace_back(std::move(walked.second));
            }
            else {
                match_set.walks2.emplace_back(std::move(walked.second));
            }
        }
        // the counts are initially set to the number of walks
        match_set.count1 = match_set.walks1.size();
        match_set.count2 = match_set.walks2.size();
    }
    
    
    return match_sets;
}

}

#endif /* centrolign_match_finder_hpp */
