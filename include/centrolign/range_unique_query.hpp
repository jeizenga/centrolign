#ifndef centrolign_range_unique_query_hpp
#define centrolign_range_unique_query_hpp

#include <vector>
#include <type_traits>
#include <algorithm>
#include <cctype>
#include <cassert>
#include <cstdlib>
#include <stdexcept>
#include <array>

#include "centrolign/utility.hpp"

namespace centrolign {

// based on
// https://stackoverflow.com/questions/39787455/is-it-possible-to-query-number-of-distinct-integers-in-a-range-in-olg-n

/*
 * Range Unique Query. O(n log n) space and preprocessing time, O(log n) query time
 */
template<size_t N = 2>
class RUQ {
public:
    
    template<typename T>
    RUQ(const std::vector<T>& arr);
    RUQ() = default;
    RUQ(const RUQ& other) noexcept = default;
    RUQ(RUQ&& other) noexcept = default;
    RUQ& operator=(const RUQ& other) noexcept = default;
    RUQ& operator=(RUQ&& other) noexcept = default;
    ~RUQ() = default;
    
    // returns the number of unique values in this half open interval of the array
    uint64_t range_unique(size_t begin, size_t end) const;
    
    inline size_t memory_size() const;
    
private:
    
    static const bool debug = false;
    
    struct LinkedTreeRecord {
        LinkedTreeRecord() = default;
        ~LinkedTreeRecord() = default;
        size_t value;
        // fractional cascading links to the lowest equal or greater value
        std::array<size_t, N> links;
        
        bool operator<(const LinkedTreeRecord& other) const {
            return value < other.value;
        }
    };
    
    uint64_t range_count_less(size_t begin, size_t end, size_t depth,
                              size_t lb_idx, size_t window_begin, size_t window_end) const;
    
    std::vector<std::vector<LinkedTreeRecord>> merge_tree;
    
    // the next power of N higher than the array size (convenient to not have to recompute)
    size_t virtual_size = 0;
};





/*
 * Template implementations
 */

template<size_t N>
template <class T>
RUQ<N>::RUQ(const std::vector<T>& arr) {
    
    if (debug) {
        std::cerr << "beginning construction algorithm\n";
    }
    
    // TODO: this is only for the initial map, which i could also use an
    // unordered_map for...
    static_assert(std::is_integral<T>::value && std::is_unsigned<T>::value,
                  "RUQ can only be built for unsigned integer types");
    static_assert(N >= 2, "RUQ N-ary tree can only be built for N >= 2");
    
    if (arr.empty()) {
        return;
    }
    
    // round up the next power of N for the depth
    size_t log_base_N = 0;
    virtual_size = 1;
    while (virtual_size < arr.size()) {
        ++log_base_N;
        virtual_size *= N;
    }
    merge_tree.resize(log_base_N + 1);
    for (auto& level : merge_tree) {
        level.resize(arr.size());
    }
    
    if (debug) {
        std::cerr << "using " << merge_tree.size() << " levels to sort virtual array of length " << virtual_size << " for real array of size " << arr.size() << '\n';
    }
    
    // init the bottom layer
    {
        auto& bottom_level = merge_tree.back();
        
        // find the index of the next occurrence of each value
        std::vector<size_t> occurence_map;
        for (int64_t i = arr.size() - 1; i >= 0; --i) {
            while (occurence_map.size() <= arr[i]) {
                occurence_map.push_back(arr.size());
            }
            bottom_level[i].value = occurence_map[arr[i]];
            occurence_map[arr[i]] = i;
        }
    }
    
    if (debug) {
        std::cerr << "finished initializing bottom layer\n";
    }
    
    size_t window = N;
    for (int64_t depth = merge_tree.size() - 2; depth >= 0; --depth) {
        if (debug) {
            std::cerr << "building layer at depth " << depth << "\n";
        }
        
        // build the next level up in the tree
        auto& level = merge_tree[depth];
        auto& next_level = merge_tree[depth + 1];
        
        for (size_t w = 0; w < level.size(); w += window) {
            
            // the end of a subwindow
            std::array<size_t, N> end;
            // the index from the window that is currently loaded in the heap
            std::array<size_t, N> curr;
            size_t stride = window / N;
            curr[0] = w;
            end[0] = std::min(w + stride, level.size());
            for (size_t i = 1; i < N; ++i) {
                curr[i] = std::min(curr[i - 1] + stride, level.size());
                end[i] = std::min(end[i - 1] + stride, level.size());
            }
            // the lagging indexes that we use for the link identification
            std::array<size_t, N> curr_link = curr;
            
            // initialize the heap
            std::array<std::pair<size_t, size_t>, N> heap;
            auto heap_end = heap.end();
            for (size_t i = 0; i < N; ++i) {
                if (curr[i] < level.size()) {
                    heap[i] = std::make_pair(next_level[curr[i]].value, i);
                }
                else {
                    --heap_end;
                }
            }
            std::make_heap(heap.begin(), heap_end, std::greater<std::pair<size_t, size_t>>());
            
            for (size_t i = w, n = end.back(); i < n; ++i) {
                // choose the next smallest value from the heap
                auto& rec = level[i];
                rec.value = heap.front().first;
                size_t source = heap.front().second;
                std::pop_heap(heap.begin(), heap_end, std::greater<std::pair<size_t, size_t>>());
                if (++curr[source] < end[source]) {
                    // there are still values left in this n-ary partition
                    (heap_end - 1)->first = next_level[curr[source]].value;
                    std::push_heap(heap.begin(), heap_end, std::greater<std::pair<size_t, size_t>>());
                }
                else {
                    // the n-ary partition is exhausted, abandon the final heap value
                    --heap_end;
                }
                
                // advance the links and save them in the record
                // TODO: could I maintain these in a heap somehow to reduce the number of j values I
                // need to check?
                for (size_t j = 0; j < N; ++j) {
                    while (curr_link[j] < end[j] && next_level[curr_link[j]].value < rec.value) {
                        ++curr_link[j];
                    }
                }
                rec.links = curr_link;
            }
        }
        window *= N;
    }
    
    if (debug) {
        std::cerr << "finished construction algorithm\n";
        for (auto& level : merge_tree) {
            std::cerr << std::string(80, '-') << '\n';
            for (size_t l = 0; l < level.size(); ++l) {
                auto& rec = level[l];
                if (l) {
                    std::cerr << '\t';
                }
                std::cerr << rec.value << " (";
                for (size_t i = 0; i < N; ++i) {
                    if (i) {
                        std::cerr << ' ';
                    }
                    std::cerr << i << ':' << rec.links[i];
                }
                std::cerr << ')';
            }
            std::cerr << '\n';
        }
    }
}

template<size_t N>
uint64_t RUQ<N>::range_unique(size_t begin, size_t end) const {
    
    if (debug) {
        std::cerr << "querying with " << begin << ", " << end << "\n";
    }
    
    // handle this as special case so that we can assume the intersection is non-empty
    if (begin >= end) {
        return 0;
    }
    
    // binary search to find the right bound
    LinkedTreeRecord search_dummy;
    search_dummy.value = end;
    auto it = std::lower_bound(merge_tree.front().begin(), merge_tree.front().end(), search_dummy);
    
    // recursive call to the merge tree
    uint64_t count_less = range_count_less(begin, end, 0, it - merge_tree.front().begin(), 0, virtual_size);
    
    // there is one unique element for each value greater or equal to 'end'
    return end - begin - count_less;
}

template<size_t N>
uint64_t RUQ<N>::range_count_less(size_t begin, size_t end, size_t depth,
                                  size_t lb_idx, size_t window_begin, size_t window_end) const {
    
    if (debug) {
        std::cerr << "recursive query " << window_begin << ":" << window_end << ", depth " << depth << ", lb idx " << lb_idx << "\n";
    }
    
    const auto& level = merge_tree[depth];
    // this is where the window actually ends, but we let the nominal window hang over the
    // edge so that we can keep it lined up with the powers of N
    size_t actual_window_end = std::min(window_end, level.size());
    
    // base case
    if (window_begin >= begin && actual_window_end <= end) {
        // the window we're in right now is fully contained
        if (debug) {
            std::cerr << "reached base case, returning upward " << lb_idx - window_begin << '\n';
        }
        return lb_idx - window_begin;
    }
    
    // find the indexes of the subwindows that contain the interval begin and end
    size_t subwindow_width = (window_end - window_begin) / N;
    // note: these cases are sufficient because we ensure non-zero overlap at each recursive call
    size_t subwindow_begin = begin > window_begin ? (begin - window_begin) / subwindow_width : 0;
    size_t subwindow_end = end < window_end ? (end - window_begin - 1) / subwindow_width + 1 : N;
    
    // recurse into subwindows that overlap with the query interval
    uint64_t count = 0;
    for (size_t i = subwindow_begin, rec_begin = window_begin + i * subwindow_width; i < subwindow_end; ++i) {
        size_t rec_end = rec_begin + subwindow_width;
        size_t next_lb_idx = lb_idx < actual_window_end ? level[lb_idx].links[i] : rec_end;
        count += range_count_less(begin, end, depth + 1, next_lb_idx, rec_begin, rec_end);
        rec_begin = rec_end;
    }
    return count;
}

template<size_t N>
inline size_t RUQ<N>::memory_size() const {
    size_t mem_size = (sizeof(merge_tree)
                       + merge_tree.capacity() * sizeof(decltype(merge_tree)::value_type));
    for (const auto& level : merge_tree) {
        mem_size += level.capacity() * sizeof(decltype(merge_tree)::value_type::value_type);
    }
    return mem_size;
}

}

#endif /* centrolign_range_unique_query_hpp */
