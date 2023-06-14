#ifndef centrolign_range_unique_query_hpp
#define centrolign_range_unique_query_hpp

#include <vector>
#include <type_traits>
#include <algorithm>

#include "centrolign/utility.hpp"

namespace centrolign {

// based on
// https://stackoverflow.com/questions/39787455/is-it-possible-to-query-number-of-distinct-integers-in-a-range-in-olg-n

/*
 * Range Unique Query. O(n log n) space and preprocessing time, O(log n) query time
 */
class RUQ {
public:
    
    template<typename T>
    RUQ(const std::vector<T>& arr);
    RUQ() = default;
    RUQ(RUQ&& other) noexcept = default;
    RUQ& operator=(RUQ&& other) noexcept = default;
    ~RUQ() = default;
    
    // returns the number of unique values in this half open interval of the array
    uint64_t range_unique(size_t begin, size_t end) const;
    
private:
    
    static const bool debug = false;
    
    struct LinkedTreeRecord {
        LinkedTreeRecord() = default;
        ~LinkedTreeRecord() = default;
        size_t value;
        // fractional cascading links to the lowest equal or greater value
        size_t left_link;
        size_t right_link;
        
        bool operator<(const LinkedTreeRecord& other) const {
            return value < other.value;
        }
    };
    
    uint64_t range_count_less(size_t begin, size_t end, size_t depth,
                              size_t lb_idx, size_t window_begin, size_t window_end) const;
    
    std::vector<std::vector<LinkedTreeRecord>> merge_tree;
};





/*
 * Template implementations
 */

template <class T>
RUQ::RUQ(const std::vector<T>& arr) {
    
    if (debug) {
        std::cerr << "beginning construction algorithm\n";
    }
    
    // TODO: this is only for the initial map, which i could also use an
    // unordered_map for...
    static_assert(std::is_integral<T>::value && std::is_unsigned<T>::value,
                  "RUQ can only be built for unsigned integer types");
    
    if (arr.empty()) {
        return;
    }
    
    // round up the next power of 2 for the depth
    merge_tree.resize(hi_bit(2 * arr.size() - 1) + 1,
                      std::vector<LinkedTreeRecord>(arr.size()));
    
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
    
    size_t window = 2;
    for (int64_t depth = merge_tree.size() - 2; depth >= 0; --depth) {
        // build the next level up in the tree
        auto& level = merge_tree[depth];
        auto& next_level = merge_tree[depth + 1];
        
        for (size_t w = 0; w < level.size(); w += window) {
            size_t n = std::min(w + window, level.size());
            size_t mid = w + window / 2;
            size_t l = w, ll = w;
            size_t r = mid, rr = mid;
            for (size_t i = w; i < n; ++i) {
                auto& rec = level[i];
                if (l < mid && r < n) {
                    // need to compare
                    if (next_level[l].value < next_level[r].value) {
                        rec.value = next_level[l++].value;
                    }
                    else {
                        rec.value = next_level[r++].value;
                    }
                }
                else if (l < mid) {
                    // right side exhausted
                    rec.value = next_level[l++].value;
                }
                else {
                    // left side exhausted
                    rec.value = next_level[r++].value;
                }
                
                // TODO: some of the if else branches make one or the other loop redundant
                // but it reads better with both of them here
                
                // get the next value links
                while (ll < mid && next_level[ll].value < rec.value) {
                    ++ll;
                }
                while (rr < n && next_level[rr].value < rec.value) {
                    ++rr;
                }
                rec.left_link = ll;
                rec.right_link = rr;
            }
        }
        window *= 2;
    }
    
    if (debug) {
        std::cerr << "finished construction algorithm\n";
        for (size_t i = 0; i < merge_tree.front().size(); ++i) {
            std::cerr << '\t' << i;
        }
        std::cerr << '\n';
        for (auto& level : merge_tree) {
            std::cerr << std::string(80, '-') << '\n';
            for (auto& rec: level) {
                std::cerr << '\t' << rec.value;
            }
            std::cerr << '\n';
            for (auto& rec: level) {
                std::cerr << "\t<" << rec.left_link;
            }
            std::cerr << '\n';
            for (auto& rec: level) {
                std::cerr << "\t>" << rec.right_link;
            }
            std::cerr << '\n';
        }
    }
}

}

#endif /* centrolign_range_unique_query_hpp */
