#ifndef centrolign_path_esa_hpp
#define centrolign_path_esa_hpp

#include <vector>
#include <cstdint>
#include <algorithm>

#include "centrolign/modify_graph.hpp"
#include "centrolign/utility.hpp"
#include "centrolign/esa.hpp"

namespace centrolign {

/*
 * Enhanced suffix array over paths in a graph
 */
class PathESA : public ESA {
public:
    
    template<class BGraph>
    PathESA(const BGraph& graph, const SentinelTableau& tableau);
    
    template<class BGraph>
    PathESA(const std::vector<const BGraph*>& graphs,
            const std::vector<const SentinelTableau*>& tableaus);
    
    std::vector<std::tuple<SANode, size_t, std::vector<uint64_t>>> minimal_rare_matches(size_t max_count) const;
    
    std::vector<std::pair<size_t, std::vector<uint64_t>>> walk_matches(const SANode& node, size_t length) const;
    
    size_t memory_size() const;
    
protected:
    
    template<class BGraph>
    PathESA(const BGraph* const* const graphs,
            const SentinelTableau* const* const tableaus, size_t size);
    
    // the lexicographic position of the ordinally next suffix
    inline size_t next(size_t i) const;
    
    // construct a suffix array with the SAIS algorithm, assumes a sentinel
    // with minimal value is included at the end of the string
    template<class StringLike>
    static std::vector<size_t> construct_suffix_array(const StringLike& str);
    
    template<class StringLike>
    static std::vector<size_t> construct_lcp_array(const StringLike& str,
                                                   const std::vector<size_t>& suffix_array,
                                                   const std::vector<size_t>& inverse_suffix_array);
    
    // the indexes distinct to this ESA variant
    std::vector<size_t> suffix_array;
    std::vector<size_t> inverse_suffix_array;
    std::string joined_seq;
    
    void print(std::ostream& out) const;
};

/*
 * Template implementations
 */

template<class BGraph>
PathESA::PathESA(const BGraph& graph, const SentinelTableau& tableau) :
    PathESA(&(&graph), std::vector<const SentinelTableau*>(1, &tableau).data(), 1)
{
    // just dispatch
}

template<class BGraph>
PathESA::PathESA(const std::vector<const BGraph*>& graphs,
                 const std::vector<const SentinelTableau*>& tableaus) :
    PathESA(graphs.data(), tableaus.data(), graphs.size())
{
    assert(graphs.size() == tableaus.size());
}

template<class BGraph>
PathESA::PathESA(const BGraph* const* const graphs,
                 const SentinelTableau* const* const tableaus, size_t size) {
    
    const bool debug = false;
    if (debug) {
        logging::level = logging::Debug;
    }
    
    // TODO: could i get rid of the tableaus and just use i + 1 as separator
    // on the corresponding component?
    std::vector<uint64_t> joined_ids;
    std::vector<size_t> index_ranges;
    for (size_t i = 0; i < size; ++i) {
        const auto& graph = *graphs[i];
        const auto& tableau = *tableaus[i];
        for (uint64_t path_id = 0; path_id < graph.path_size(); ++path_id) {
            
            joined_seq.push_back(tableau.src_sentinel + 1);
            joined_ids.push_back(tableau.src_id);
            
            for (auto node_id : graph.path(path_id)) {
                joined_seq.push_back(graph.label(node_id) + 1);
                joined_ids.push_back(node_id);
            }
            
            joined_seq.push_back(tableau.snk_sentinel + 1);
            joined_ids.push_back(tableau.snk_id);
        }
        index_ranges.push_back(joined_seq.size());
    }
    // sentinel that the SA-IS algorithm needs
    joined_seq.push_back(0);
    // arbitrarily assign it to an input graph
    joined_ids.push_back(graphs[size - 1]->node_size());
    index_ranges.back()++;
    
    logging::log(logging::Debug, "Constructing suffix array");
    
    suffix_array = std::move(construct_suffix_array(joined_seq));
    
    logging::log(logging::Debug, "Constructing LCP array");
    
    inverse_suffix_array = std::move(invert(suffix_array));
    lcp_array = std::move(construct_lcp_array(joined_seq, suffix_array, inverse_suffix_array));
    
    if (debug) {
        std::cerr << "constructing component rank tables\n";
    }
    
    // init the nearest rank vectors
    for (size_t c = 0; c < size; ++c) {
        nearest_comp_rank.emplace_back();
        nearest_comp_rank.back().resize(joined_seq.size() + 1);
    }

    component_ranked_ids.resize(size);
    
    leaf_to_comp.resize(suffix_array.size());
    for (size_t i = 0; i < suffix_array.size(); ++i) {
        for (size_t c = 0; c < nearest_comp_rank.size(); ++c) {
            nearest_comp_rank[c][i] = component_ranked_ids[c].size();
        }
        uint64_t node_id = joined_ids[suffix_array[i]];
        size_t comp = std::upper_bound(index_ranges.begin(), index_ranges.end(), suffix_array[i]) - index_ranges.begin();
        component_ranked_ids[comp].push_back(node_id);
        leaf_to_comp[i] = comp;
    }
    // add the final past-the-last entries
    for (size_t c = 0; c < size; ++c) {
        nearest_comp_rank[c][joined_seq.size()] = component_ranked_ids[c].size();
    }
    
    logging::log(logging::Debug, "Constructing child array");
    
    // build structure for traversing downward and assigning annotations
    construct_child_array();
    
    logging::log(logging::Debug, "Constructing suffix links");
    
    // get the suffix links
    auto advance = [&](size_t i) -> size_t {
        return next(i);
    };
    construct_suffix_links(advance);
    
    if (debug) {
        std::cerr << "finished construction\n";
        print(std::cerr);
    }
}

template<class StringLike>
std::vector<size_t> PathESA::construct_lcp_array(const StringLike& str,
                                                 const std::vector<size_t>& suffix_array,
                                                 const std::vector<size_t>& inverse_suffix_array) {
    
    // Kasai's algorithm
    
    std::vector<size_t> lcp_array(suffix_array.size());
    
    // lead 0 to match ESA paper
    lcp_array[0] = 0;
    
    size_t length_matched = 0;
    // note: we skip the last iteration because the end sentinel always goes to the first
    // position, so there is no previous suffix
    for (size_t i = 0, n = str.size() - 1; i < n; ++i) {
        auto pos = inverse_suffix_array[i];
        auto j = suffix_array[pos - 1];
        while (str[i + length_matched] == str[j + length_matched]) {
            ++length_matched;
        }
        lcp_array[pos] = length_matched;
        if (length_matched != 0) {
            // one fewer position will match when we advance in the string
            --length_matched;
        }
    }
    
    return lcp_array;
}


template<class StringLike>
std::vector<size_t> PathESA::construct_suffix_array(const StringLike& str) {
    
    const bool debug = false;
    
    // steps 2 and 3 in both the induction and the full sorting steps (algs 3.3 and 3.4)
    auto finish_sort = [](const StringLike& str,
                          std::vector<size_t>& suffix_array, const std::vector<bool>& s_type,
                          std::vector<size_t>& bucket_ends, std::vector<size_t>& bucket_begins) {
        
        // left-to-right scan putting L type suffixes in their position
        for (size_t i = 0; i < suffix_array.size(); ++i) {
            if (suffix_array[i] != 0 && !s_type[suffix_array[i] - 1]) {
                suffix_array[bucket_begins[str[suffix_array[i] - 1]]++] = suffix_array[i] - 1;
            }
        }
        
        if (debug) {
            std::cerr << "after left to right scan:\n";
            for (size_t i = 0; i < suffix_array.size(); ++i) {
                std::cerr << '\t' << i << ": " << suffix_array[i] << '\n';
            }
        }
        
        // right-to-left scan putting S type suffixes in their position
        for (int64_t i = suffix_array.size() - 1; i >= 0; --i) {
            if (suffix_array[i] != 0 && s_type[suffix_array[i] - 1]) {
                suffix_array[--bucket_ends[str[suffix_array[i] - 1]]] = suffix_array[i] - 1;
            }
        }
        
        if (debug) {
            std::cerr << "after right to left scan:\n";
            for (size_t i = 0; i < suffix_array.size(); ++i) {
                std::cerr << '\t' << i << ": " << suffix_array[i] << '\n';
            }
        }
    };
    
    // we assume the sentinel has already been added
    std::vector<size_t> suffix_array(str.size(), 0);
    
    // figure out where each of character buckets end
    std::vector<size_t> bucket_ends;
    for (auto bucket : str) {
        while (bucket_ends.size() <= bucket) {
            bucket_ends.push_back(0);
        }
        ++bucket_ends[bucket];
    }
    // consolidate and also check for the base case
    bool fully_bucket_sorted = true;
    for (size_t i = 1; i < bucket_ends.size(); ++i) {
        fully_bucket_sorted &= (bucket_ends[i] <= 1);
        bucket_ends[i] += bucket_ends[i - 1];
    }
    // note: we assume without checking that the sentinel (in position 0) has count 1
    
    if (fully_bucket_sorted) {
        if (debug) {
            std::cerr << "input is fully bucket-sorted, ending recursion\n";
        }
        
        // end of recursion, just carry out the bucket sort
        for (size_t i = 0; i < str.size(); ++i) {
            suffix_array[bucket_ends[str[i]] - 1] = i;
        }
    }
    else {
        // we need to construct the LMS suffix array and induce a sort
        
        // identify the S and L type positions
        std::vector<size_t> lms_positions;
        std::vector<bool> s_type(str.size());
        s_type.back() = true;
        for (int64_t i = str.size() - 2; i >= 0; --i) {
            if (str[i] == str[i + 1]) {
                s_type[i] = s_type[i + 1];
            }
            else {
                s_type[i] = (str[i] < str[i + 1]);
                if (!s_type[i] && s_type[i + 1]) {
                    lms_positions.push_back(i + 1);
                }
            }
        }
        
        // put this into ascending order
        std::reverse(lms_positions.begin(), lms_positions.end());
        
        if (debug) {
            std::cerr << "types:\n";
            for (auto s : s_type) {
                std::cerr << ' ' << (s ? 'S' : 'L');
            }
            std::cerr << '\n';
        }
        
        // we will be adding to buckets from both sides
        std::vector<size_t> bucket_begins(bucket_ends.size(), 0);
        for (size_t i = 1; i < bucket_begins.size(); ++i) {
            bucket_begins[i] = bucket_ends[i - 1];
        }
            
        // get a sort that correctly orders the LMS substrings (but possibly not other strings)
        {
            {
                // add the LMS positions to the ends of the buckets in arbitrary order
                auto bucket_ends_copy = bucket_ends;
                for (auto pos : lms_positions) {
                    suffix_array[--bucket_ends_copy[str[pos]]] = pos;
                }
                
                if (debug) {
                    std::cerr << "state after adding LMS to buckets:\n";
                    for (size_t i = 0; i < suffix_array.size(); ++i) {
                        std::cerr << '\t' << i << ": " << suffix_array[i] << '\n';
                    }
                }
            }
            // induce the sort (using copies because we will need these again later)
            auto bucket_ends_copy = bucket_ends;
            auto bucket_begins_copy = bucket_begins;
            finish_sort(str, suffix_array, s_type, bucket_ends_copy, bucket_begins_copy);
        }
        
        // convert LMS strings into a new alphabet
        
        std::vector<uint32_t> lms_names(str.size(), -1);
        // start with the sentinel (which will be at position 0 in the SA)
        lms_names.back() = 0;
        uint32_t curr_name = 0;
        size_t prev_lms_pos = str.size() - 1;
        // continue with the rest of the SA
        for (size_t i = 1; i < suffix_array.size(); ++i) {
            auto pos = suffix_array[i];
            if (pos != 0 && s_type[pos] && !s_type[pos - 1]) {
                // LMS position
                bool full_match = false;
                // check for match in sequence and in S/L type
                auto p1 = prev_lms_pos, p2 = pos;
                while (!full_match && str[p1] == str[p2] && s_type[p1] == s_type[p2]) {
                    ++p1;
                    ++p2;
                    // did we reach another LMS position? (at which point we stop)
                    full_match = (s_type[p1] && !s_type[p1 - 1] && s_type[p2] && !s_type[p2 - 1]);
                }
                
                // TODO: it's possible to check full bucket sorting here (if we ever re-use a name)
                if (!full_match) {
                    // we need another name
                    ++curr_name;
                }
                lms_names[pos] = curr_name;
                
                prev_lms_pos = pos;
            }
        }
        // move the names into the front of the vector
        size_t removed = 0;
        for (size_t i = 0; i < lms_names.size(); ++i) {
            if (lms_names[i] == -1) {
                ++removed;
            }
            else {
                lms_names[i - removed] = lms_names[i];
            }
        }
        lms_names.resize(lms_names.size() - removed);
        
        if (debug) {
            std::cerr << "derived name array for recursive call:\n";
            for (auto n : lms_names) {
                std::cerr << ' ' << n;
            }
            std::cerr << '\n';
        }
        
        // recursively determine the suffix array for these names
        std::vector<size_t> lms_suffix_array = construct_suffix_array(lms_names);
        
        // reset the suffix array ahead of the final sort
        for (size_t i = 0; i < suffix_array.size(); ++i) {
            suffix_array[i] = 0;
        }
        
        // put the LMS positions at the end of their buckets
        {
            auto bucket_ends_copy = bucket_ends;
            for (int64_t i = lms_suffix_array.size() - 1; i >= 0; --i) {
                size_t pos = lms_positions[lms_suffix_array[i]];
                suffix_array[--bucket_ends_copy[str[pos]]] = pos;
            }
            if (debug) {
                std::cerr << "state after sorted LMS to buckets:\n";
                for (size_t i = 0; i < suffix_array.size(); ++i) {
                    std::cerr << '\t' << i << ": " << suffix_array[i] << '\n';
                }
            }
        }
        // for the last sort, we can destructively use up the bucket ends/begins
        finish_sort(str, suffix_array, s_type, bucket_ends, bucket_begins);
    }
    
    return suffix_array;
}

inline size_t PathESA::next(size_t i) const {
    auto pos = suffix_array[i];
    if (pos + 1 < suffix_array.size()) {
        return inverse_suffix_array[pos + 1];
    }
    else {
        return -1;
    }
}


}

#endif /* centrolign_path_esa_hpp */
