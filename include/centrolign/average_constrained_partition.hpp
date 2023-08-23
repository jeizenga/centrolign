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
    
    
    static const bool debug = false;
    
    static const T mininf = std::numeric_limits<T>::lowest();
    
    if (debug) {
        std::cerr << "starting avg constrained partition, setting up prefix sum vectors\n";
    }
    
    // to compute summed score
    std::vector<T> prefix_sum(data.size());
    // to check minimum average
    std::vector<T> fractional_prefix_sum(data.size());
    if (!data.empty()) {
        prefix_sum.front() = data.front().first;
        fractional_prefix_sum.front() = data.front().first - data.front().second * min_average;
    }
    for (size_t i = 1; i < data.size(); ++i) {
        prefix_sum[i] = prefix_sum[i - 1] + data[i].first;
        fractional_prefix_sum[i] = fractional_prefix_sum[i - 1] + data[i].first - data[i].second * min_average;
    }
    
    if (debug) {
        std::cerr << "prefix sum:\n";
        for (auto v : prefix_sum) {
            std::cerr << ' ' << v;
        }
        std::cerr << '\n';
        std::cerr << "fractional prefix sum:\n";
        for (auto v : fractional_prefix_sum) {
            std::cerr << ' ' << v;
        }
        std::cerr << '\n';
    }
    
    if (debug) {
        std::cerr << "initializing sparse RMQ\n";
    }
    
    
    // records of (score if excluded, score if included)
    std::vector<std::pair<T, T>> dp(data.size() + 1, std::make_pair(mininf, mininf));
    // back pointers for the score with the current item included
    std::vector<int64_t> backpointer(dp.size(), -1);
    
    
    // initialize the sparse range max query
    std::vector<std::pair<std::pair<T, size_t>, T>> tree_data;
    tree_data.reserve(data.size() + 1);
    for (size_t i = 0; i < data.size(); ++i) {
        tree_data.emplace_back(std::make_pair(fractional_prefix_sum[i], i + 1), mininf);
    }
    // the boundary condition
    tree_data.emplace_back(std::pair<T, size_t>(0, 0), 0);
    dp.front().first = 0;
    
    MaxSearchTree<std::pair<T, size_t>, T> search_tree(tree_data);
    
    size_t opt_idx = -1;
    
    if (debug) {
        std::cerr << "doing DP\n";
    }
    
    // main dynamic programming
    for (size_t i = 1; i < dp.size(); ++i) {
        if (debug) {
            std::cerr << "DP iteration " << i << '\n';
        }
        
        // find max if i is not included
        dp[i].first = std::max(dp[i - 1].first, dp[i - 1].second);
        
        // find max if i is included in the partition
        auto max_it = search_tree.range_max(std::pair<T, size_t>(mininf, 0),
                                            std::pair<T, size_t>(fractional_prefix_sum[i - 1], -1));
        if (max_it != search_tree.end() && max_it->second != mininf) {
            // we've completed DP for a valid interval start
            dp[i].second = prefix_sum[i - 1] + max_it->second;
            backpointer[i] = max_it->first.second;
            if (opt_idx == -1 || dp[i].second > dp[opt_idx].second) {
                opt_idx = i;
            }
            if (debug) {
                std::cerr << "found opt partner index at " << max_it->first.second << " with combined DP score of " << dp[i].second << '\n';
            }
        }
        
        // enter the query values for intervals starting on next cell
        auto it = search_tree.find(std::make_pair(fractional_prefix_sum[i - 1], i));
        search_tree.update(it, dp[i].first - prefix_sum[i - 1]);
        if (debug) {
            std::cerr << "record interval start value of " << (dp[i].first - prefix_sum[i - 1]) << " at fractional prefix " << fractional_prefix_sum[i - 1] << '\n';
        }
    }
    
    // do traceback
    
    if (debug) {
        std::cerr << "doing traceback from opt " << opt_idx << "\n";
    }
    
    std::vector<std::pair<size_t, size_t>> partition;
    if (opt_idx != -1) {
        // we found at least one valid interval
        bool in_interval = true;
        int64_t tb_idx = opt_idx;
        
        while (tb_idx > 0) {
            if (in_interval){
                // follow the backpointer and emit an interval
                int64_t prev_idx = backpointer[tb_idx];
                partition.emplace_back(prev_idx, tb_idx);
                tb_idx = prev_idx;
                in_interval = false;
            }
            else {
                // did we get this value from a finished interval or a further gap?
                in_interval = (dp[tb_idx].first == dp[tb_idx - 1].second);
                --tb_idx;
            }
        }
    }
    
    // put in the natural, front-to-back order
    std::reverse(partition.begin(), partition.end());
    
    return partition;
}

}

#endif /* centrolign_average_constrained_partition_hpp */
