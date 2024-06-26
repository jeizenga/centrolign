#ifndef centrolign_integer_sort_hpp
#define centrolign_integer_sort_hpp

#include <vector>
#include <cstdint>
#include <iostream>

namespace centrolign {

// returns stably-sorted vector of indexes based on a rank function using
// linear time algorithm (assuming max rank is O(n))
template<class RankGetter>
std::vector<size_t> integer_sort(const std::vector<size_t>& indexes,
                                 const RankGetter& getter) {

    // count up occurrences of each rank
    std::vector<size_t> counts;
    for (size_t i : indexes) {
        size_t next_rank = getter(i) + 1;
        while (counts.size() <= next_rank) {
            counts.push_back(0);
        }
        ++counts[next_rank];
    }
    
    for (size_t i = 1; i < counts.size(); ++i) {
        counts[i] += counts[i - 1];
    }
    
    std::vector<size_t> sorted(indexes.size());
    for (size_t i : indexes) {
        size_t rank = getter(i);
        sorted[counts[rank]++] = i;
    }
    
    return sorted;
}

}

#endif /* centrolign_integer_sort_hpp */
