#include "centrolign/bonder.hpp"

#include <cassert>
#include <limits>
#include <iomanip>

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
    
    static const bool debug = false;
    
    double opt_pref_sum = 0.0;
    double sec_pref_sum = 0.0;
    std::vector<double> length_prefix_sum(shared_subanchors.size() + 1, 0.0);
    std::vector<double> excl_length_prefix_sum(shared_subanchors.size() + 1, 0.0);
    std::vector<double> fractional_difference(shared_subanchors.size() + 1, 0.0);
    std::vector<double> excl_fractional_difference(shared_subanchors.size() + 1, 0.0);
    
    std::vector<std::pair<std::pair<double, size_t>, double>> tree_data;
    tree_data.reserve(shared_subanchors.size() + 1);
    tree_data.emplace_back(std::pair<double, size_t>(0.0, 0), 0.0);
    
    for (size_t i = 0; i < shared_subanchors.size(); ++i) {
        
        if (i != 0) {
            
            excl_length_prefix_sum[i] = length_prefix_sum[i] + std::get<0>(intervening_segments[i - 1]);
            
            opt_pref_sum += std::get<1>(intervening_segments[i - 1]);
            sec_pref_sum += std::get<2>(intervening_segments[i - 1]);
            excl_fractional_difference[i] = min_opt_proportion * opt_pref_sum - sec_pref_sum;
            if (debug) {
                std::cerr << "(bw) i " << i << ", opt sum " << opt_pref_sum << ", sec sum " << sec_pref_sum << ", excl frac diff " << excl_fractional_difference[i] << '\n';
            }
            
            tree_data.emplace_back(std::pair<double, size_t>(excl_fractional_difference[i], i), mininf);
        }
        
        
        length_prefix_sum[i + 1] += excl_length_prefix_sum[i] + std::get<0>(shared_subanchors[i]);
        opt_pref_sum += std::get<1>(shared_subanchors[i]);
        sec_pref_sum += std::get<2>(shared_subanchors[i]);
        
        fractional_difference[i + 1] = min_opt_proportion * opt_pref_sum - sec_pref_sum;
        
        if (debug) {
            std::cerr << "i " << i << ", opt sum " << opt_pref_sum << ", sec sum " << sec_pref_sum << ", frac diff " << fractional_difference[i + 1] << '\n';
        }
        
//        tree_data.emplace_back(std::pair<double, size_t>(fractional_difference[i + 1], i + 1), mininf);
        
    }
    
    if (debug) {
        std::cerr << std::setprecision(3);
        std::cerr << "i\tLP\tELP\tFD\tEFD\tSh\tIn\n";
        for (size_t i = 0; i < fractional_difference.size(); ++i) {
            std::cerr << i << '\t' << length_prefix_sum[i] << '\t' << excl_length_prefix_sum[i] << '\t' << fractional_difference[i] << '\t' << excl_fractional_difference[i];
            if (i != 0) {
                std::cerr << '\t' << std::get<0>(shared_subanchors[i - 1]) << ',' << std::get<1>(shared_subanchors[i - 1]) << ',' << std::get<2>(shared_subanchors[i - 1]);
            }
            if (i != 0 && i + 1 != fractional_difference.size()) {
                std::cerr << '\t' << std::get<0>(intervening_segments[i - 1]) << ',' << std::get<1>(intervening_segments[i - 1]) << ',' << std::get<2>(intervening_segments[i - 1]);
            }
            std::cerr << '\n';
        }
    }
    
    MaxSearchTree<std::pair<double, size_t>, double> tree(tree_data);
    
    // records of (excluded, included)
    std::vector<std::pair<double, double>> dp(length_prefix_sum.size(), std::make_pair(mininf, mininf));
    dp.front().first = 0.0;
    dp.front().second = 0.0;
    std::vector<size_t> backpointer(dp.size(), -1);
    
    for (size_t i = 1; i < dp.size(); ++i) {
        
        dp[i].first = std::max(dp[i - 1].first, dp[i - 1].second);
        
        if (debug) {
            std::cerr << "at " << i << ", querying with " << fractional_difference[i] << '\n';
            std::cerr << "RMQ state\n";
            for (auto kv : tree) {
                std::cerr << '\t' << kv.first.first << ": " << kv.second << '\t' << "(" << kv.first.second << ")\n";
            }
        }
        
        auto it = tree.range_max(std::pair<double, size_t>(fractional_difference[i], 0), inf);
        
        if (it != tree.end() && it->second != mininf) {
            if (debug) {
                std::cerr << "got predecessor " << it->first.second << " with key " << it->first.first << " and val " << it->second << '\n';
            }
            
            backpointer[i] = it->first.second;
            dp[i].second = length_prefix_sum[i] + it->second - min_length;
        }
        else if (debug) {
            std::cerr << "did not find any results\n";
        }
        
        if (i + 1 != dp.size()) {
            tree.update(tree.find(std::pair<double, size_t>(excl_fractional_difference[i], i)),
                        dp[i].first - excl_length_prefix_sum[i]);
        }
    }
    
    if (debug) {
        std::cerr << "final DP:\n";
        for (size_t i = 0; i < dp.size(); ++i) {
            std::cerr << '\t' << i << '\t' << dp[i].first << '\t' << dp[i].second << '\n';
        }
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
