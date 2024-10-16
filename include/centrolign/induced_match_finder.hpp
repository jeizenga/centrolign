#ifndef centrolign_induced_match_finder_hpp
#define centrolign_induced_match_finder_hpp

#include <vector>
#include <utility>
#include <iostream>
#include <unordered_map>
#include <unordered_set>

#include "centrolign/graph.hpp"
#include "centrolign/modify_graph.hpp"
#include "centrolign/match_finder.hpp"
#include "centrolign/step_index.hpp"

namespace centrolign {

class InducedMatchFinderComponentView; // forward declaration

/*
 * Object that projects self-matches from a full graph into 2-disconnected components
 */
class InducedMatchFinder {
public:
    
    InducedMatchFinder(const BaseGraph& full_graph, const std::vector<match_set_t>& matches,
                       const std::vector<std::pair<uint64_t, uint64_t>>& components, const StepIndex& step_index);
    
    InducedMatchFinder() = default;
    ~InducedMatchFinder() = default;
    
    // return the match finder for one of the components (designated by its index in the input vector)
    InducedMatchFinderComponentView component_view(size_t comp) const;
    
private:
    
    friend class InducedMatchFinderComponentView;
    
    // the location of a matching set on different paths
    struct path_hit_set_t {
        path_hit_set_t() = default;
        path_hit_set_t(const path_hit_set_t& other) noexcept = default;
        path_hit_set_t(path_hit_set_t&& other) noexcept = default;
        ~path_hit_set_t() = default;
        path_hit_set_t& operator=(const path_hit_set_t& other) noexcept = default;
        path_hit_set_t& operator=(path_hit_set_t&& other) noexcept = default;
        
        // by path, the offsets where a match occurs and the index of origin walk in the match set
        std::unordered_map<uint64_t, std::vector<std::pair<size_t, size_t>>> hit_locations;
        // the length of the full match
        size_t length = 0;
        // the count in the entire graph
        size_t deduplicated_count = 0;
    };
    
    const BaseGraph* parent = nullptr;
    
    std::vector<std::vector<path_hit_set_t>> component_path_hits;
};


/*
 * Object that projects self-matches onto a particular 2-disconnected component
 */
class InducedMatchFinderComponentView {
public:
    
    ~InducedMatchFinderComponentView() = default;
    
    // returns minimal rare matches
    template<class BGraph>
    std::vector<match_set_t> find_matches(const BGraph& graph1, const BGraph& graph2,
                                          const SentinelTableau& tableau1, const SentinelTableau& tableau2) const;
    
private:
    
    friend class InducedMatchFinder;
    
    InducedMatchFinderComponentView() = default;
    InducedMatchFinderComponentView(const BaseGraph* parent, const std::vector<InducedMatchFinder::path_hit_set_t>* path_hits) : parent(parent), path_hits(path_hits) {}
    
    // FIXME: copied from Core -- very bad segmentation
    std::tuple<std::string, size_t, size_t> parse_subpath_name(const std::string& subpath_name) const;

    
    const BaseGraph* parent = nullptr;
    const std::vector<InducedMatchFinder::path_hit_set_t>* path_hits = nullptr;
};




/**
 * Template implementations
 */



template<class BGraph>
std::vector<match_set_t> InducedMatchFinderComponentView::find_matches(const BGraph& graph1, const BGraph& graph2,
                                                                       const SentinelTableau& tableau1, const SentinelTableau& tableau2) const {
    
    static const bool debug = false;
    
    std::unordered_set<uint64_t> parent_path_seen;
    size_t parent_path_length1 = 0, parent_path_length2 = 0;
    // translate between the path coordinates (path ID, begin offset, end offset)
    std::vector<std::tuple<uint64_t, size_t, size_t>> path_trans1, path_trans2;
    for (bool do_1 : {true, false}) {
        auto& path_trans = do_1 ? path_trans1 : path_trans2;
        const auto& graph = do_1 ? graph1 : graph2;
        auto& path_length = do_1 ? parent_path_length1 : parent_path_length2;
        path_trans.resize(graph.path_size());
        for (uint64_t path_id = 0; path_id < graph.path_size(); ++path_id) {
            
            std::string path_name;
            size_t path_begin, path_end;
            std::tie(path_name, path_begin, path_end) = parse_subpath_name(graph.path_name(path_id));
            
            uint64_t parent_path_id = parent->path_id(path_name);
            path_trans[path_id] = std::make_tuple(parent_path_id, path_begin, path_end);
            
            if (parent_path_seen.insert(parent_path_id).second) {
                // this is the first time we've seen this full path
                path_length += parent->path(parent_path_id).size();
            }
        }
    }
    
    double approx_count_ratio = double(parent_path_length1) / double(parent_path_length2);
    
    auto assign_count = [](size_t observed1, size_t observed2, size_t target_count, double approx_ratio_12) -> std::pair<size_t, size_t> {
        size_t count2 = round(sqrt(target_count / approx_ratio_12));
        size_t count1 = round(sqrt(target_count * approx_ratio_12));
        if (count1 >= observed1 && count2 < observed2) {
            // count2 would ideally be smaller than the number observed
            count2 = observed2;
            count1 = round(target_count / double(count2));
        }
        else if (count2 >= observed2 && count1 < observed1) {
            // count1 would ideally be smaller than the number observed
            count1 = observed1;
            count2 = round(target_count / double(count1));
        }
        // make sure we count at least as many as we actually saw
        return std::make_pair(std::max(count1, observed1), std::max(count2, observed2));
    };
    
    
    // we'll fill out match sets as we go
    std::vector<match_set_t> matches;
    
    for (const auto& path_hit_set : *path_hits) {
        
        std::unordered_set<size_t> origin_walks_used;
        
        // records of (match begin, match end, on graph 1, path id, path offset)
        std::vector<std::tuple<size_t, size_t, bool, uint64_t, size_t>> intervals;
        size_t observed1 = 0, observed2 = 0;
        for (bool do_1 : {true, false}) {
            const auto& graph = do_1 ? graph1 : graph2;
            const auto& path_trans = do_1 ? path_trans1 : path_trans2;
            auto& observed = do_1 ? observed1 : observed2;
            
            // TODO: i don't love trying every path on every match set, but is there an alternative?
            // if we iterate over the path matches that actually exist, we could have many more paths in the full
            // graph than in a smaller subproblem...
            // the ideal might be to switch off between the iteration orders depending on path_size?
            
            
            
            // the combinations of (node ID, match index) that we have already seen to deduplicate by starting location
            std::unordered_set<std::pair<uint64_t, size_t>> initial_nodes;
            
            for (uint64_t path_id = 0; path_id < graph.path_size(); ++path_id) {
                
                uint64_t parent_path_id;
                size_t path_begin, path_end;
                std::tie(parent_path_id, path_begin, path_end) = path_trans[path_id];
                
                auto it = path_hit_set.hit_locations.find(parent_path_id);
                if (it == path_hit_set.hit_locations.end()) {
                    continue;
                }
                
                // the range of locations that correspond to the subpath through this component
                // note: have to avoid underflow when subtracting off the length
                auto loc_begin = std::lower_bound(it->second.begin(), it->second.end(),
                                                  std::pair<size_t, size_t>(path_begin >= path_hit_set.length ? path_begin - path_hit_set.length : 0, 0));
                auto loc_end = std::upper_bound(it->second.begin(), it->second.end(), std::pair<size_t, size_t>(path_end, 0));
                
                for (auto loc_it = loc_begin; loc_it != loc_end; ++loc_it) {
                                        
                    size_t match_begin = loc_it->first;
                    size_t match_end = match_begin + path_hit_set.length;
                    
                    origin_walks_used.insert(loc_it->second);
                    
                    // the portion overlapping the component
                    size_t begin = match_begin < path_begin ? (path_begin - match_begin) : 0;
                    size_t end = match_end > path_end ? (match_end - path_end) : path_hit_set.length;
                    
                    // the offset on the subproblem's path
                    size_t path_offset = match_begin < path_begin ? 0 : match_begin - path_begin;
                    
                    uint64_t node_id = graph.path(path_id)[path_offset];
                    
                    if (initial_nodes.emplace(node_id, begin).second) {
                        // we haven't seen a match starting here before
                        intervals.emplace_back(begin, end, do_1, path_id, path_offset);
                    }
                    
                    ++observed;
                }
            }
        }
        
        // estimate the global count on the two sequences (with a hack)
        size_t total_count = observed1 * observed2 + path_hit_set.deduplicated_count - origin_walks_used.size();
        size_t count1, count2;
        std::tie(count1, count2) = assign_count(observed1, observed2, total_count, approx_count_ratio);
        if (debug) {
            std::cerr << "assigning counts " << count1 << ", " << count2 << " based on observed " << observed1 << ", " << observed2 << " total " << total_count << ", interval count " << intervals.size() << ", dedupe " << path_hit_set.deduplicated_count << ", walks used " << origin_walks_used.size() << ", ratio " << approx_count_ratio << '\n';
        }
        
        // sort by starting index
        std::sort(intervals.begin(), intervals.end());
        
        if (debug) {
            std::cerr << "convering intervals into matches:\n";
            for (size_t i = 0; i < intervals.size(); ++i) {
                std::cerr << i << ": " << std::get<0>(intervals[i]) << ", " << std::get<1>(intervals[i]) << ", " << std::get<2>(intervals[i]) << ", " << std::get<3>(intervals[i]) << ", " << std::get<4>(intervals[i]) << "\n";
            }
        }
        
        // a heap ordered by the ending index
        std::vector<size_t> active_intervals;
        auto cmp = [&](size_t i, const size_t j) {
            return std::get<1>(intervals[i]) > std::get<1>(intervals[j]);
        };
        
        size_t last = 0; // the previous event
        size_t i = 0; // the index of the next interval among the starts
        // the number of active intervals from each graph
        size_t num_active1 = 0;
        size_t num_active2 = 0;
        while (i < intervals.size() || !active_intervals.empty()) {
            
            // is the next event the start or end of interval(s)?
            bool next_is_start;
            size_t next;
            if (active_intervals.empty() ||
                (!active_intervals.empty() && i < intervals.size() && std::get<0>(intervals[i]) < std::get<1>(intervals[active_intervals.front()]))) {
                next_is_start = true;
                next = std::get<0>(intervals[i]);
            }
            else {
                next_is_start = false;
                next = std::get<1>(intervals[active_intervals.front()]);
            }
            
            if (debug) {
                std::cerr << "interval start " << i << " of " << intervals.size() << ", " << active_intervals.size() << " active intervals with nearest end from " << int(active_intervals.empty() ? -1 : active_intervals.front()) << "\n";
                if (i < intervals.size()) {
                    std::cerr << "next interval start " << std::get<0>(intervals[i]) << ", " << std::get<1>(intervals[i]) << ", " << std::get<2>(intervals[i]) << ", " << std::get<3>(intervals[i]) << ", " << std::get<4>(intervals[i]) << "\n";
                }
                if (!active_intervals.empty()) {
                    std::cerr << "next active interval end " << std::get<0>(intervals[active_intervals.front()]) << ", " << std::get<1>(intervals[active_intervals.front()]) << ", " << std::get<2>(intervals[active_intervals.front()]) << ", " << std::get<3>(intervals[active_intervals.front()]) << ", " << std::get<4>(intervals[active_intervals.front()]) << "\n";
                }
                std::cerr << "next is start? " << next_is_start << ", next at " << next << '\n';
            }
            
            if (num_active1 != 0 && num_active2 != 0 && next != last) {
                // we can emit a match for the currently active intervals before moving to the next event
                if (debug) {
                    std::cerr << "emitting a match for " << num_active1 << " and " << num_active2 << " active intervals from the 2 graphs\n";
                }
                
                matches.emplace_back();
                auto& match_set = matches.back();
                for (size_t idx : active_intervals) {
                    const auto& interval = intervals[idx];
                    auto& walks = std::get<2>(interval) ? match_set.walks1 : match_set.walks2;
                    const auto& graph = std::get<2>(interval) ? graph1 : graph2;
                    const auto& path = graph.path(std::get<3>(interval));
                    
                    walks.emplace_back();
                    auto& walk = walks.back();
                    size_t begin = std::get<4>(interval) + (last - std::get<0>(interval));
                    size_t end = begin + (next - last);
                    for (size_t k = begin; k < end; ++k) {
                        walk.push_back(path[k]);
                    }
                }
                match_set.full_length = path_hit_set.length;
                match_set.count1 = count1;
                match_set.count2 = count2;
            }
            
            last = next;
            
            if (next_is_start) {
                // the next event is a new active interval
                
                // find the range that of intervals that start here
                size_t j = i + 1;
                while (j < intervals.size() && std::get<0>(intervals[j]) == std::get<0>(intervals[i])) {
                    ++j;
                }
                for (size_t k = i; k < j; ++k) {
                    active_intervals.emplace_back(k);
                    if (std::get<2>(intervals[k])) {
                        ++num_active1;
                    }
                    else {
                        ++num_active2;
                    }
                    std::push_heap(active_intervals.begin(), active_intervals.end(), cmp);
                }
                
                i = j;
            }
            else {
                // the next event is the end of an active interval
                
                // move all intervals that end here to the back of the heap vector
                auto heap_end = active_intervals.end();
                std::pop_heap(active_intervals.begin(), heap_end--, cmp);
                
                while (heap_end != active_intervals.begin() &&
                       std::get<1>(intervals[active_intervals.front()]) == std::get<1>(intervals[active_intervals.back()])) {
                    std::pop_heap(active_intervals.begin(), heap_end--, cmp);
                }
                
                for (auto it = heap_end; it != active_intervals.end(); ++it) {
                    if (std::get<2>(intervals[*it])) {
                        --num_active1;
                    }
                    else {
                        --num_active2;
                    }
                }
                
                active_intervals.resize(heap_end - active_intervals.begin());
            }
        }
    }
    
    if (debug) {
        std::cerr << "final matches:\n";
        for (size_t i = 0; i < matches.size(); ++i) {
            std::cerr << "match set " << i << ", full length " << matches[i].full_length << ", count1 " << matches[i].count1 << ", count2 " << matches[i].count2 << '\n';
            std::cerr << "walks on 1:\n";
            for (auto w : matches[i].walks1) {
                std::cerr << '\t';
                for (auto n : w) {
                    std::cerr << n << ' ';
                }
                std::cerr << '\n';
            }
            std::cerr << "walks on 2:\n";
            for (auto w : matches[i].walks2) {
                std::cerr << '\t';
                for (auto n : w) {
                    std::cerr << n << ' ';
                }
                std::cerr << '\n';
            }
        }
    }
    
    return matches;
}

}

#endif /* centrolign_induced_match_finder_hpp */
