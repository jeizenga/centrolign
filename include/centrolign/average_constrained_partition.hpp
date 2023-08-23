#ifndef centrolign_average_constrained_partition_hpp
#define centrolign_average_constrained_partition_hpp

#include "centrolign/max_search_tree.hpp"

#include <vector>
#include <utility>
#include <limits>

namespace centrolign {


// takes a vector of (score, weight) pairs as input
// returns a set of half-open intervals that maximize the sum of scores, subject
// to the constraint that each interval has a weighted average value above a minimum
template<class T>
std::vector<std::pair<size_t, size_t>> average_constrained_partition(const std::vector<std::pair<T, T>>& data,
                                                                     const T& min_average) {
    
    
    static const T mininf = std::numeric_limits<T>::lowest();
    
    // to compute summed score
    std::vector<T> prefix_sum(data.size());
    // to check minimum average
    std::vector<T> fractional_prefix_sum(data.size());
    if (!data.empty()) {
        prefix_sum.front() = data.front().first;
        fractional_prefix_sum.front() += data.front().first - data.front().second * min_average
    }
    for (size_t i = 1; i < data.size(); ++i) {
        prefix_sum[i] = prefix_sum[i - 1] + data[i].first;
        fractional_prefix_sum[i] = fractional_prefix_sum[i - 1] + data[i].first - data[i].second * min_average;
    }
    
    // initialize the sparse range max query
    std::vector<std::pair<std::pair<T, size_t>, T>> tree_data;
    tree_data.reserve(data.size());
    for (size_t i = 0; i < data.size(); ++i) {
        tree_data.emplace_back(std::make_pair(fractional_prefix_sum[i], i), mininf);
    }
    MaxSearchTree<std::pair<T, size_t>, T> search_tree(tree_data);
    
    // records of (score if excluded, score if included)
    std::vector<std::pair<T, T>> dp(data.size(), std::make_pair(mininf, mininf));
    // back pointers for the score with the current item included
    // note: the traceback step depends on -1 being the exact initialization here
    std::vector<int64_t> backpointer(dp.size(), -1);
    
    size_t opt_idx = -1;
    
    // boundary condition
    if (!dp.empty()) {
        // no score if don't include first
        dp.front().first = 0;
        
        // score of lone first value, if feasible
        if (fractional_prefix_sum.front() >= 0) {
            dp.front().second = data.front().first;
            opt_idx = 0;
        }
        
        // init query structure
        auto it = search_tree.find(std::make_pair(fractional_prefix_sum.front(), (size_t) 0));
        search_tree.update(it, 0);
    }
    
    // main dynamic programming
    for (size_t i = 1; i < data.size(); ++i) {
        
        // find max if i is not included
        dp[i].first = std::max(dp[i - 1].first, dp[i - 1].second);
        
        // find max if i is included in the partition
        auto max_it = search_tree.range_max(std::pair<T, size_t>(mininf, 0),
                                            std::pair<T, size_t>(fractional_prefix_sum[i + 1], -1));
        if (max_it != search_tree.end() && max_it->second != mininf) {
            // we've completed DP for a valid interval start
            dp[i].second = prefix_sum[i] + max_it->second;
            backpointer[i] = max_it->first.second;
            if (opt_idx == -1 || dp[i].second > dp[opt_idx].second) {
                opt_idx = i;
            }
        }
        
        // enter the query values for intervals starting on next cell
        auto it = search_tree.find(std::make_pair(fractional_prefix_sum[i], i));
        search_tree.update(it, dp[i].first - prefix_sum[i]);
    }
    
    // do traceback
    
    std::vector<std::pair<size_t, size_t>> partition;
    if (opt_idx != -1) {
        // we found at least one valid interval
        bool in_interval = true;
        int64_t tb_idx = opt_idx;
        
        while (tb_idx >= 0) {
            if (in_interval){
                int64_t prev_idx = backpointer[tb_idx];
                // note: this relies on the initialization being -1 to handle the case that
                // the interval starts on 0
                partition.emplace_back(prev_idx + 1, tb_idx + 1); // +1 for past-the-last
                tb_idx = prev_idx;
                in_interval = false;
            }
            else {
                if (tb_idx != 0) {
                    // did we get this value from a finished interval or a gap?
                    in_interval = (dp[tb_idx].first == dp[tb_idx - 1].second);
                }
                --tb_idx;
            }
        }
    }
    
    std::reverse(partition.begin(), partition.end());
    
    return partition;
}

}

#endif /* centrolign_average_constrained_partition_hpp */
