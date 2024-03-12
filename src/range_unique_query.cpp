#include "centrolign/range_unique_query.hpp"

#include <cctype>
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <stdexcept>

namespace centrolign {

using namespace std;

uint64_t RUQ::range_unique(size_t begin, size_t end) const {
    
    if (debug) {
        std::cerr << "querying with " << begin << ", " << end << "\n";
    }
    
    // handle this as special cases so that we can assume the intersection is non-empty
    if (begin >= end) {
        return 0;
    }
    
    // binary search to find the right bound
    LinkedTreeRecord search_dummy;
    search_dummy.value = end;
    auto it = std::lower_bound(merge_tree.front().begin(), merge_tree.front().end(), search_dummy);
    
    // recursive call to the merge tree
    uint64_t count_less = range_count_less(begin, end, 0, it - merge_tree.front().begin(),
                                           0, 1 << (merge_tree.size() - 1));
    
    // there is one unique element for each value greater or equal to 'end'
    return end - begin - count_less;
}


uint64_t RUQ::range_count_less(size_t begin, size_t end, size_t depth,
                               size_t lb_idx, size_t window_begin, size_t window_end) const {
    
    if (debug) {
        std::cerr << "recursive query " << window_begin << ":" << window_end << ", depth " << depth << ", lb idx " << lb_idx << "\n";
    }
    
    const auto& level = merge_tree[depth];
    // this is where the window actually ends, but we let the nominal window hang over the
    // edge so that we can keep it lined up with the powers of 2
    size_t actual_window_end = std::min(window_end, level.size());
    
    // base case
    if (window_begin >= begin && actual_window_end <= end) {
        // the window we're in right now is fully contained
        if (debug) {
            std::cerr << "reached base case, returning upward " << lb_idx - window_begin << '\n';
        }
        return lb_idx - window_begin;
    }
    
    size_t window_mid = (window_begin + window_end) / 2;
    
    uint64_t count = 0;
    if (window_mid > begin && window_begin < end) {
        // there is an intersection on the left side
        size_t next_lb_idx = lb_idx != actual_window_end ? level[lb_idx].left_link : window_mid;
        count += range_count_less(begin, end, depth + 1, next_lb_idx, window_begin, window_mid);
    }
    if (actual_window_end > begin && window_mid < end) {
        // there is an intersection of the right side
        size_t next_lb_idx = lb_idx != actual_window_end ? level[lb_idx].right_link : window_end;
        count += range_count_less(begin, end, depth + 1, next_lb_idx, window_mid, window_end);
    }
    return count;
}

size_t RUQ::memory_size() const {
    size_t mem_size = (sizeof(merge_tree) + merge_tree.capacity() * sizeof(std::vector<LinkedTreeRecord>));
    for (const auto& level : merge_tree) {
        mem_size += level.capacity() * sizeof(LinkedTreeRecord);
    }
    return mem_size;
}


}
