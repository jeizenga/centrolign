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
#include <limits>

#include "centrolign/utility.hpp"

namespace centrolign {

// based on
// https://stackoverflow.com/questions/39787455/is-it-possible-to-query-number-of-distinct-integers-in-a-range-in-olg-n

/*
 * Range Unique Query. O(n log n) space and preprocessing time, O(log n) query time.
 * Optionally uses an N-ary tree with sub-sampling to tune memory consumption against
 * query time.
 */
template<typename UIntSize = size_t, size_t NAry = 2, size_t Sampling = 1>
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
    uint64_t range_unique(UIntSize begin, UIntSize end) const;
    
    inline size_t memory_size() const;
    
private:
    
    static const bool debug = false;
    
    inline UIntSize num_window_sampled_links(UIntSize window_size) const;
    
    uint64_t range_count_less(UIntSize begin, UIntSize end, UIntSize depth,
                              UIntSize lb_idx, UIntSize window_begin, UIntSize window_end) const;
    
    // the top level of the merge sort tree
    std::vector<UIntSize> sorted_next_occurrences;
    // the (possibly subsampled) fractional cascading links
    std::vector<std::vector<std::array<UIntSize, NAry>>> cascading_links;
    // the lower layers that are used (if subsampling) to reconstruct the unsampled cascading links
    std::vector<std::vector<UIntSize>> reconstruction_layers;
    
    // the next power of N higher than the array size (convenient to not have to recompute)
    size_t virtual_size = 0;
};





/*
 * Template implementations
 */

template<typename UIntSize, size_t NAry, size_t Sampling>
template <class T>
RUQ<UIntSize, NAry, Sampling>::RUQ(const std::vector<T>& arr) {
    
    if (debug) {
        std::cerr << "beginning construction algorithm\n";
    }
    
    // TODO: this is only for the initial map, which i could also use an
    // unordered_map for...
    static_assert(std::is_integral<T>::value && std::is_unsigned<T>::value,
                  "RUQ can only be built for unsigned integer types");
    static_assert(NAry >= 2, "RUQ N-ary tree can only be built for N >= 2");
    static_assert(Sampling != 0, "RUQ Sampling must be non-zero");
    
    if (arr.empty()) {
        return;
    }
    
    // round up the next power of N for the depth
    size_t log_base_N = 0;
    virtual_size = 1;
    while (virtual_size < arr.size()) {
        ++log_base_N;
        virtual_size *= NAry;
    }
    
    if (size_t(std::numeric_limits<UIntSize>::max()) < virtual_size) {
        throw std::runtime_error("RUQ integer size is too small to index input.");
    }
    
    // the bottom layer on the merge tree doesn't need fractional cascading links
    cascading_links.resize(log_base_N);
    for (size_t i = 0, window_size = virtual_size; i < cascading_links.size(); ++i) {
        if (Sampling == 1) {
            // the easy case, no subsampling
            cascading_links[i].resize(arr.size());
        }
        else {
            // figure out how many subsampled links there will be
            size_t num_full_windows = arr.size() / window_size;
            size_t samples = (num_full_windows * num_window_sampled_links(window_size)
                              + num_window_sampled_links(arr.size() - num_full_windows * window_size));
            cascading_links[i].resize(samples);
            window_size /= NAry;
        }
    }
    
    if (debug) {
        std::cerr << "using " << log_base_N << " levels to sort virtual array of length " << virtual_size << " for real array of size " << arr.size() << '\n';
    }
    
    // init the bottom layer
    std::vector<UIntSize> next_level(arr.size());
    {
        // find the index of the next occurrence of each value
        std::vector<UIntSize> occurence_map;
        for (int64_t i = arr.size() - 1; i >= 0; --i) {
            while (occurence_map.size() <= arr[i]) {
                occurence_map.push_back(arr.size());
            }
            next_level[i] = occurence_map[arr[i]];
            occurence_map[arr[i]] = i;
        }
    }
    
    if (debug) {
        std::cerr << "finished initializing bottom layer\n";
    }
    
    UIntSize window = NAry;
    for (int64_t depth = log_base_N - 1; depth >= 0; --depth) {
        if (debug) {
            std::cerr << "building layer at depth " << depth << "\n";
        }
        
        // build the next level up in the tree
        std::vector<UIntSize> level(arr.size());
        auto& link_layer = cascading_links[depth];
        
        for (UIntSize w = 0, link_idx = 0; w < level.size(); w += window) {
            
            // the end of a subwindow
            std::array<UIntSize, NAry> end;
            // the index from the window that is currently loaded in the heap
            std::array<UIntSize, NAry> curr;
            UIntSize stride = window / NAry;
            curr[0] = w;
            end[0] = std::min<UIntSize>(w + stride, level.size());
            for (size_t i = 1; i < NAry; ++i) {
                curr[i] = std::min<UIntSize>(curr[i - 1] + stride, level.size());
                end[i] = std::min<UIntSize>(end[i - 1] + stride, level.size());
            }
            // the lagging indexes that we use for the link identification
            std::array<UIntSize, NAry> curr_link = curr;
            
            // initialize the heap with records of (value, which subwindow)
            std::array<std::pair<UIntSize, UIntSize>, NAry> heap;
            auto heap_end = heap.end();
            for (UIntSize i = 0; i < NAry; ++i) {
                if (curr[i] < level.size()) {
                    heap[i] = std::make_pair(next_level[curr[i]], i);
                }
                else {
                    --heap_end;
                }
            }
            std::make_heap(heap.begin(), heap_end, std::greater<std::pair<UIntSize, UIntSize>>());
            
            for (UIntSize i = w, n = end.back(); i < n; ++i) {
                // choose the next smallest value from the heap
                level[i] = heap.front().first;
                UIntSize source = heap.front().second;
                std::pop_heap(heap.begin(), heap_end, std::greater<std::pair<UIntSize, UIntSize>>());
                if (++curr[source] < end[source]) {
                    // there are still values left in this n-ary partition
                    (heap_end - 1)->first = next_level[curr[source]];
                    std::push_heap(heap.begin(), heap_end, std::greater<std::pair<UIntSize, UIntSize>>());
                }
                else {
                    // the n-ary partition is exhausted, abandon the final heap value
                    --heap_end;
                }
                
                // advance the links and save them in the record
                // TODO: could I maintain these in a heap somehow to reduce the number of j values I
                // need to check?
                for (size_t j = 0; j < NAry; ++j) {
                    while (curr_link[j] < end[j] && next_level[curr_link[j]] < level[i]) {
                        ++curr_link[j];
                    }
                }
                if (Sampling == 1) {
                    // we're not subsampling at all, so the indexing is easy
                    link_layer[i] = curr_link;
                }
                else {
                    // only save every n-th set of links
                    if (((i - w) % Sampling) == 0) {
                        link_layer[link_idx++] = curr_link;
                    }
                }
            }
        }
        window *= NAry;
        if (Sampling != 1) {
            reconstruction_layers.push_back(std::move(next_level));
        }
        next_level = std::move(level);
    }
    
    // the array should be fully sorted now
    sorted_next_occurrences = std::move(next_level);
    
    if (Sampling > 1) {
        // reorder the reconstruction layers from highest to lowest
        std::reverse(reconstruction_layers.begin(), reconstruction_layers.end());
    }
    
    if (debug) {
        std::cerr << "finished construction algorithm\n";
        std::cerr << "sorted occurrences:\n";
        for (size_t i = 0; i < sorted_next_occurrences.size(); ++i) {
            if (i) {
                std::cerr << '\t';
            }
            std::cerr << sorted_next_occurrences[i];
        }
        std::cerr << '\n';
        std::cerr << "links:\n";
        for (const auto& layer : cascading_links) {
            for (size_t i = 0; i < layer.size(); ++i) {
                if (i) {
                    std::cerr << '\t';
                }
                std::cerr << '(';
                for (size_t j = 0; j < NAry; ++j) {
                    if (j) {
                        std::cerr << ',';
                    }
                    std::cerr << layer[i][j];
                }
                std::cerr << ')';
            }
            std::cerr << '\n';
        }
        std::cerr << "reconstruction layers (if any):\n";
        for (const auto& layer : reconstruction_layers) {
            for (size_t i = 0; i < layer.size(); ++i) {
                if (i) {
                    std::cerr << '\t';
                }
                std::cerr << layer[i];
            }
            std::cerr << '\n';
        }
    }
}

template<typename UIntSize, size_t NAry, size_t Sampling>
uint64_t RUQ<UIntSize, NAry, Sampling>::range_unique(UIntSize begin, UIntSize end) const {
    
    if (debug) {
        std::cerr << "querying with " << begin << ", " << end << "\n";
    }
    
    // handle this as special case so that we can assume the intersection is non-empty
    if (begin >= end) {
        return 0;
    }
    
    // binary search to find the right bound
    auto it = std::lower_bound(sorted_next_occurrences.begin(), sorted_next_occurrences.end(), end);
    
    // recursive call to the merge tree
    uint64_t count_less = range_count_less(begin, end, 0, it - sorted_next_occurrences.begin(), 0, virtual_size);
    
    // there is one unique element for each value greater or equal to 'end'
    return end - begin - count_less;
}

template<typename UIntSize, size_t NAry, size_t Sampling>
uint64_t RUQ<UIntSize, NAry, Sampling>::range_count_less(UIntSize begin, UIntSize end, UIntSize depth,
                                                         UIntSize lb_idx, UIntSize window_begin, UIntSize window_end) const {
    
    if (debug) {
        std::cerr << "recursive query " << window_begin << ":" << window_end << ", depth " << depth << ", lb idx " << lb_idx << "\n";
    }
    
    // this is where the window actually ends, but we let the nominal window hang over the
    // edge so that we can keep it lined up with the powers of N
    UIntSize actual_window_end = std::min<UIntSize>(window_end, sorted_next_occurrences.size());
    
    // base case
    if (window_begin >= begin && actual_window_end <= end) {
        // the window we're in right now is fully contained
        if (debug) {
            std::cerr << "reached base case, returning upward " << lb_idx - window_begin << '\n';
        }
        return lb_idx - window_begin;
    }
    
    const auto& link_level = cascading_links[depth];
    
    // find the indexes of the subwindows that contain the interval begin and end
    UIntSize window_width = window_end - window_begin;
    UIntSize subwindow_width = window_width / NAry;
    // note: these cases are sufficient because we ensure non-zero overlap at each recursive call
    UIntSize subwindow_begin = begin > window_begin ? (begin - window_begin) / subwindow_width : 0;
    UIntSize subwindow_end = end < window_end ? (end - window_begin - 1) / subwindow_width + 1 : NAry;
    
    // find the index of the nearest sampled link (if necessary)
    UIntSize prev_link_sampled = 0;
    if (Sampling != 1) {
        UIntSize window_link_width = num_window_sampled_links(window_width);
        UIntSize window_link_begin = window_link_width * (window_begin / window_width);
        UIntSize idx_in_window = lb_idx - window_begin;
        prev_link_sampled = window_link_begin + idx_in_window / Sampling;
        if (debug) {
            std::cerr << "window links begin at index " << window_link_begin << ", previous sampled link is " << prev_link_sampled << '\n';
        }
    }
    
    // recurse into subwindows that overlap with the query interval
    uint64_t count = 0;
    for (UIntSize i = subwindow_begin, rec_begin = window_begin + i * subwindow_width; i < subwindow_end; ++i) {
        UIntSize rec_end = rec_begin + subwindow_width;
        UIntSize next_lb_idx;
        if (Sampling == 1) {
            // no subsampling, so we can just grab the value directly
            next_lb_idx = lb_idx < actual_window_end ? link_level[lb_idx][i] : rec_end;
        }
        else {
            // we may need to recompute the exact link
            if (lb_idx == actual_window_end) {
                next_lb_idx = rec_end;
            }
            else {
                const auto& curr_layer = depth ? reconstruction_layers[depth - 1]  : sorted_next_occurrences;
                const auto& next_layer = reconstruction_layers[depth];
                
                // extend from the previous the previous sampled value
                next_lb_idx = link_level[prev_link_sampled][i];
                if (debug) {
                    std::cerr << "in subwindow " << i << " of query " << window_begin << ":" << window_end << ", get initial link to " << next_lb_idx << '\n';
                }
                while (next_lb_idx < rec_end && next_layer[next_lb_idx] < curr_layer[lb_idx]) {
                    ++next_lb_idx;
                    if (debug) {
                        std::cerr << "advance link to " << next_lb_idx << '\n';
                    }
                }
            }
        }
        count += range_count_less(begin, end, depth + 1, next_lb_idx, rec_begin, rec_end);
        rec_begin = rec_end;
    }
    return count;
}

template<typename UIntSize, size_t NAry, size_t Sampling>
inline UIntSize RUQ<UIntSize, NAry, Sampling>::num_window_sampled_links(UIntSize window_size) const {
    return (window_size + Sampling - 1) / Sampling;
}

template<typename UIntSize, size_t NAry, size_t Sampling>
inline size_t RUQ<UIntSize, NAry, Sampling>::memory_size() const {
    size_t mem_size = (sizeof(sorted_next_occurrences)
                       + sorted_next_occurrences.capacity() * sizeof(typename decltype(sorted_next_occurrences)::value_type)
                       + sizeof(cascading_links)
                       + cascading_links.capacity() * sizeof(typename decltype(cascading_links)::value_type)
                       + sizeof(reconstruction_layers)
                       + reconstruction_layers.capacity() * sizeof(typename decltype(reconstruction_layers)::value_type));
    for (const auto& level : cascading_links) {
        mem_size += level.capacity() * sizeof(typename decltype(cascading_links)::value_type::value_type);
    }
    for (const auto& level : reconstruction_layers) {
        mem_size += reconstruction_layers.capacity() * sizeof(typename decltype(reconstruction_layers)::value_type::value_type);
    }
    return mem_size;
}

}

#endif /* centrolign_range_unique_query_hpp */
