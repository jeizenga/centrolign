#ifndef centrolign_match_finder_hpp
#define centrolign_match_finder_hpp

#include <deque>
#include <limits>
#include <iostream>
#include <sstream>

#include "centrolign/graph.hpp"
#include "centrolign/modify_graph.hpp"
#include "centrolign/path_esa.hpp"
#include "centrolign/gesa.hpp"
#include "centrolign/score_function.hpp"
#include "centrolign/utility.hpp"
#include "centrolign/simplifier.hpp"

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
    size_t full_length = 0;
};

/*
 *
 */

/*
 * Base class for an object that finds matches between two graphs
 */
class BaseMatchFinder {
public:
    ~BaseMatchFinder() = default;
    
    /*
     * Configurable parameters
     */
    
    // the max count in either of the two graphs
    size_t max_count = 50;
    // use either the Color Set Size index or Range Unique Query index to count match occurrences
    bool use_color_set_size = true;
    
protected:
    
    BaseMatchFinder(const ScoreFunction& score_function) : score_function(&score_function) {}
    BaseMatchFinder() = default;
     
    static const bool debug_match_finder = false;
    
    const ScoreFunction* const score_function = nullptr;
    
    
    template<class Index>
    std::vector<match_set_t> query_index(const Index& index) const;
};

/*
 * Object that finds matches between two graphs using a suffix array of the path sequences
 */
class PathMatchFinder : public BaseMatchFinder {
public:
    PathMatchFinder(const ScoreFunction& score_function) : BaseMatchFinder(score_function) {}
    PathMatchFinder() = default;
    ~PathMatchFinder() = default;
    
    template<class BGraph>
    std::vector<match_set_t> find_matches(const BGraph& graph1, const BGraph& graph2,
                                          const SentinelTableau& tableau1, const SentinelTableau& tableau2) const;
private:
    
    
    static const bool debug_match_finder = false;
    
};

/*
 * Object that finds matches between two graphs using an uncompressed GCSA
 */
class GESAMatchFinder : public BaseMatchFinder {
public:
    GESAMatchFinder(const ScoreFunction& score_function) : BaseMatchFinder(score_function) {}
    GESAMatchFinder() = default;
    ~GESAMatchFinder() = default;
    
    template<class BGraph>
    std::vector<match_set_t> find_matches(const BGraph& graph1, const BGraph& graph2,
                                          const SentinelTableau& tableau1, const SentinelTableau& tableau2) const;
    
    // throw a GESASizeError if it grows to this times as much as the graph size
    size_t size_limit_factor = 16;
    
    // simplifies graph topology in advance of querying matches
    Simplifier simplifier;
    
private:
    
    
    static const bool debug_match_finder = false;
    
    std::vector<match_set_t> index_and_query(ExpandedGraph& expanded1,
                                             ExpandedGraph& expanded2) const;
};

/*
 * Template and inline implementations
 */

template<class BGraph>
std::vector<match_set_t> PathMatchFinder::find_matches(const BGraph& graph1, const BGraph& graph2,
                                                       const SentinelTableau& tableau1, const SentinelTableau& tableau2) const {
    
    std::vector<const BGraph*> graph_ptrs{&graph1, &graph2};
    std::vector<const SentinelTableau*> tableau_ptrs{&tableau1, &tableau2};
    
    PathESA path_esa(graph_ptrs, tableau_ptrs);
    
    return query_index(path_esa);
}

template<class BGraph>
std::vector<match_set_t> GESAMatchFinder::find_matches(const BGraph& graph1, const BGraph& graph2,
                                                       const SentinelTableau& tableau1, const SentinelTableau& tableau2) const {
    
    // do an initial simplification
    ExpandedGraph expanded1 = simplifier.simplify(graph1, tableau1);
    ExpandedGraph expanded2 = simplifier.simplify(graph2, tableau2);
    
    return index_and_query(expanded1, expanded2);
}

template<class Index>
std::vector<match_set_t> BaseMatchFinder::query_index(const Index& index) const {
    
    logging::log(logging::Debug, "Finding minimal rare matches");
    if (logging::level >= logging::Debug) {
        logging::log(logging::Debug, "Match index is occupying " + format_memory_usage(index.memory_size()) + ".");
        logging::log(logging::Debug, "Current memory usage is " + format_memory_usage(current_memory_usage()) + ".");
    }
    
    // records of (min count on either graph, total pairs, length, node)
    std::vector<std::tuple<size_t, size_t, size_t, SANode>> matches;
    size_t kept_num_pairs = 0;
    size_t removed_num_pairs = 0;
    size_t total_length = 0;
    for (const auto& match : index.minimal_rare_matches(max_count, use_color_set_size)) {
        
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
            total_length += (counts[0] + counts[1]) * std::get<1>(match);
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
        match_set.full_length = match_set.walks1.front().size();
    }
    
    if (logging::level >= logging::Debug) {
        logging::log(logging::Debug, "Walked matches currently consuming " + format_memory_usage(total_length * sizeof(uint64_t)) + " of memory.");
        logging::log(logging::Debug, "Current memory usage is " + format_memory_usage(current_memory_usage()) + ".");
    }
    
    return match_sets;
}

}

#endif /* centrolign_match_finder_hpp */
