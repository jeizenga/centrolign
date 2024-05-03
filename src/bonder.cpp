#include "centrolign/bonder.hpp"

#include <cassert>
#include <limits>

#include "centrolign/max_search_tree.hpp"

namespace centrolign {

// passed in as records of (length, opt segment score, secondary segment score)
std::vector<std::pair<size_t, size_t>>
Bonder::longest_partition(const std::vector<std::tuple<double, double, double>>& shared_subanchors,
                          const std::vector<std::tuple<double, double, double>>& intervening_segments) const {
    
    // because the two segments on the end would never be selected
    assert(intervening_segments.size() == shared_subanchors.size() - 1);
    
    static const double mininf = std::numeric_limits<double>::lowest();
    static const std::pair<double, size_t> inf = std::pair<double, size_t>(std::numeric_limits<double>::max(),
                                                                           std::numeric_limits<size_t>::max());
    
    double opt_pref_sum = 0.0;
    double sec_pref_sum = 0.0;
    std::vector<double> length_prefix_sum(shared_subanchors.size() + 1, 0.0);
    std::vector<double> fractional_difference(shared_subanchors.size() + 1, 0.0);
    
    std::vector<std::pair<std::pair<double, size_t>, double>> tree_data;
    tree_data.reserve(shared_subanchors.size() + 1);
    tree_data.emplace_back(std::pair<double, size_t>(0.0, 0), 0.0);
    
    for (size_t i = 0; i < shared_subanchors.size(); ++i) {
        
        if (i != 0) {
            length_prefix_sum[i + 1] = std::get<0>(intervening_segments[i - 1]);
            opt_pref_sum += std::get<1>(intervening_segments[i - 1]);
            sec_pref_sum += std::get<2>(intervening_segments[i - 1]);
        }
        
        length_prefix_sum[i + 1] += length_prefix_sum[i] + std::get<0>(shared_subanchors[i]);
        opt_pref_sum += std::get<1>(shared_subanchors[i]);
        sec_pref_sum += std::get<2>(shared_subanchors[i]);
        
        fractional_difference[i + 1] = sec_pref_sum - min_opt_proportion * opt_pref_sum;
        
        tree_data.emplace_back(std::pair<double, size_t>(fractional_difference[i + 1], i + 1), mininf);
        
    }
    
    MaxSearchTree<std::pair<double, size_t>, double> tree(tree_data);
    
    // records of (excluded, included)
    std::vector<std::pair<double, double>> dp(length_prefix_sum.size(), std::make_pair(mininf, mininf));
    dp.front().first = 0.0;
    dp.front().second = 0.0;
    std::vector<size_t> backpointer(dp.size(), -1);
    
    for (size_t i = 1; i < dp.size(); ++i) {
        
        dp[i].first = std::max(dp[i - 1].first, dp[i - 1].second);
        
        auto it = tree.range_max(std::pair<double, size_t>(fractional_difference[i], 0), inf);
        if (it != tree.end() && it->second != mininf) {
            backpointer[i] = it->first.second;
            dp[i].second = length_prefix_sum[i] + it->second;
        }
        
        tree.update(tree.find(std::pair<double, size_t>(fractional_difference[i], i)),
                    dp[i].first - length_prefix_sum[i]);
    }
        
    size_t tb_idx = 0;
    for (size_t i = 1; i < dp.size(); ++i) {
        if (dp[i].second > dp[tb_idx].second) {
            tb_idx = i;
        }
    }
    
    return traceback(dp, backpointer, tb_idx);
}


}
