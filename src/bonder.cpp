#include "centrolign/bonder.hpp"

#include <cassert>
#include <limits>
#include <iomanip>

#include "centrolign/max_search_tree.hpp"

namespace centrolign {

int Bonder::output_num = 0;

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
                std::cerr << "(bw) i " << i << ", opt sum " << opt_pref_sum << ", sec sum " << sec_pref_sum << ", excl frac diff " << excl_fractional_difference[i] << ", excl len pref " << excl_length_prefix_sum[i] << '\n';
            }
            
            tree_data.emplace_back(std::pair<double, size_t>(excl_fractional_difference[i], i), mininf);
        }
        
        
        length_prefix_sum[i + 1] = excl_length_prefix_sum[i] + std::get<0>(shared_subanchors[i]);
        opt_pref_sum += std::get<1>(shared_subanchors[i]);
        sec_pref_sum += std::get<2>(shared_subanchors[i]);
        
        fractional_difference[i + 1] = min_opt_proportion * opt_pref_sum - sec_pref_sum;
        
        if (debug) {
            std::cerr << "i " << i << ", opt sum " << opt_pref_sum << ", sec sum " << sec_pref_sum << ", frac diff " << fractional_difference[i + 1] << ", len pref " << length_prefix_sum[i + 1] << '\n';
        }
    }
    
    if (debug) {
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
    
    size_t tb_idx = 0;
    for (size_t i = 1; i < dp.size(); ++i) {
        
        dp[i].first = std::max(dp[i - 1].first, dp[i - 1].second);
        
        if (debug) {
            std::cerr << "at " << i << ", querying with " << fractional_difference[i] << '\n';
//            std::cerr << "RMQ state\n";
//            for (auto kv : tree) {
//                std::cerr << '\t' << kv.first.first << ": " << kv.second << '\t' << "(" << kv.first.second << ")\n";
//            }
        }
        
        auto it = tree.range_max(std::pair<double, size_t>(fractional_difference[i], 0), inf);
        
        if (it != tree.end() && it->second != mininf) {
            if (debug) {
                std::cerr << "got predecessor " << it->first.second << " with key " << it->first.first << " and val " << it->second << '\n';
            }
            
            backpointer[i] = it->first.second;
            dp[i].second = length_prefix_sum[i] + it->second - min_length;
            if (dp[i].second > dp[tb_idx].second) {
                tb_idx = i;
                if (debug) {
                    std::cerr << "this is the new traceback opt with score " << dp[i].second << "\n";
                }
            }
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
    
    return traceback(dp, backpointer, tb_idx);
}

std::vector<std::pair<size_t, size_t>>
Bonder::longest_windowed_partition(const std::vector<std::tuple<double, double, double>>& shared_subanchors,
                                   const std::vector<std::tuple<double, double, double>>& intervening_segments) const {
    
    
    static const bool debug = false;
    
    assert(intervening_segments.size() == shared_subanchors.size() - 1);
    
    static const double mininf = std::numeric_limits<double>::lowest();
    static const std::pair<double, size_t> inf = std::pair<double, size_t>(std::numeric_limits<double>::max(),
                                                                           std::numeric_limits<size_t>::max());
    
    size_t joined_size = (shared_subanchors.size() + intervening_segments.size());
    std::vector<bool> meets_constraint_left_adj(joined_size);
    std::vector<bool> meets_constraint_right_adj(joined_size);
    std::vector<int64_t> leftward_partner(joined_size);
    std::vector<int64_t> rightward_partner(joined_size);
    
    for (bool forward : {true, false}) {
        // the score and weight of the full interval at the current iteration
        double curr_window_length = 0.0;
        double window_opt_score = 0.0;
        double window_sec_score = 0.0;
        // the past-the-last position of the interval of this iteration
        int64_t end = forward ? 0 : joined_size - 1;
        int64_t incr = forward ? 1 : -1;
        auto& meets_constraint = forward ? meets_constraint_left_adj : meets_constraint_right_adj;
        auto& partner = forward ? rightward_partner : leftward_partner;
        for (int64_t i = end; i < partner.size() && i >= 0; i += incr) {
            while (end < partner.size() && end >= 0 && curr_window_length < window_length) {
                bool between = (end % 2 == 1);
                const auto& segment = between ? intervening_segments[end / 2] : shared_subanchors[end / 2];
                curr_window_length += std::get<0>(segment);
                window_opt_score += std::get<1>(segment);
                window_sec_score += std::get<2>(segment);
                end += incr;
            }
            partner[i] = end;
            
            if ((end < 0 || end >= partner.size()) && curr_window_length < window_length) {
                // we don't have a full window anymore, so the previous full interval applies to the
                // rest of the vector
                meets_constraint[i] = meets_constraint[i - incr];
            }
            else {
                double final_length, final_opt_score, final_sec_score;
                if (end % 2 == 0) {
                    std::tie(final_length, final_opt_score, final_sec_score) = intervening_segments[(end - incr) / 2];
                }
                else {
                    std::tie(final_length, final_opt_score, final_sec_score) = shared_subanchors[(end - incr) / 2];
                }
                
                double opt_window_score = (window_opt_score - final_opt_score +
                                           (window_length - (curr_window_length - final_length)) / final_length * final_opt_score);
                double sec_window_score = (window_sec_score - final_sec_score +
                                           (window_length - (curr_window_length - final_length)) / final_length * final_sec_score);
                
                meets_constraint[i] = (sec_window_score > min_opt_proportion * opt_window_score);
                if (i % 2 == 1 && partner[i] == partner[i] + incr) {
                    // the entire window is in between two shared segments, which we treat as a violation of the constraint
                    meets_constraint[i] = false;
                }
            }
            
            const auto& segment = (i % 2 == 1) ? intervening_segments[i / 2] : shared_subanchors[i / 2];
            curr_window_length -= std::get<0>(segment);
            window_opt_score -= std::get<1>(segment);
            window_sec_score -= std::get<2>(segment);
        }
    }
    
    double opt_pref_sum = 0.0;
    double sec_pref_sum = 0.0;
    std::vector<double> length_prefix_sum(shared_subanchors.size() + 1, 0.0);
    std::vector<double> excl_length_prefix_sum(shared_subanchors.size() + 1, 0.0);
    std::vector<double> fractional_difference(shared_subanchors.size() + 1, 0.0);
    std::vector<double> excl_fractional_difference(shared_subanchors.size() + 1, 0.0);
    
    std::vector<int> left_adj_constraint_prefix_sum(shared_subanchors.size() + 1, 0);
    std::vector<int> right_adj_constraint_prefix_sum(shared_subanchors.size() + 1, 0);
    
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
                std::cerr << "(bw) i " << i << ", opt sum " << opt_pref_sum << ", sec sum " << sec_pref_sum << ", excl frac diff " << excl_fractional_difference[i] << ", excl len pref " << excl_length_prefix_sum[i] << '\n';
            }
            
            tree_data.emplace_back(std::pair<double, size_t>(excl_fractional_difference[i], i), mininf);
        }
        
        
        length_prefix_sum[i + 1] = excl_length_prefix_sum[i] + std::get<0>(shared_subanchors[i]);
        opt_pref_sum += std::get<1>(shared_subanchors[i]);
        sec_pref_sum += std::get<2>(shared_subanchors[i]);
        
        fractional_difference[i + 1] = min_opt_proportion * opt_pref_sum - sec_pref_sum;
        
        bool meets_left = meets_constraint_left_adj[2 * i];
        bool meets_right = meets_constraint_right_adj[2 * i];
        if (i != 0) {
            meets_left = meets_left || meets_constraint_left_adj[2 * i - 1];
            meets_right = meets_right || meets_constraint_right_adj[2 * i - 1];
        }
        left_adj_constraint_prefix_sum[i + 1] = left_adj_constraint_prefix_sum[i] + int(meets_left);
        right_adj_constraint_prefix_sum[i + 1] = right_adj_constraint_prefix_sum[i] + int(meets_right);
    }
    
    MaxSearchTree<std::pair<double, size_t>, double> tree(tree_data);
    
    // records of (excluded, included)
    std::vector<std::pair<double, double>> dp(length_prefix_sum.size(), std::make_pair(mininf, mininf));
    dp.front().first = 0.0;
    dp.front().second = 0.0;
    std::vector<size_t> backpointer(dp.size(), -1);
    
    size_t tb_idx = 0;
    
    // for maintaining the shorter-than-window sparse query structure
    size_t window_begin = 0;
    double curr_window_length = 0.0;
    
    // these start as null until an index falls out of the active window
    size_t outside_window_argmax = -1;
    size_t argmax_partner = -1;
    // index of first right-adjusted window to fully the right of the points outside the active window
    size_t k = 0;
    // index of the first left-adjusted window that extends beyond the active window
    size_t l = 0;
    
    // figure out where we need to stop iterating l to get to the last full window, but no further
    size_t final_l = shared_subanchors.size();
    {
        double tail_length = 0.0;
        while (final_l != 0 && tail_length + std::get<0>(shared_subanchors[final_l - 1]) < window_length) {
            tail_length += std::get<0>(shared_subanchors[final_l - 1]);
            --final_l;
        }
    }
    if (debug) {
        std::cerr << "l iteration will be stopped at " << final_l << '\n';
    }
    
    for (size_t i = 1; i < dp.size(); ++i) {
        
        
        if (debug) {
            std::cerr << "DP iteration " << i << '\n';
        }
        
        // find the first position whose left-adjusted window finishes after i's right side
        // note 1: the partner is past-the-last, and i is 1 greater than the data index
        // note 2: we stop after the first interval that reaches the end so that we don't get confused by small partial
        // windows that don't belong in the longer-than-window calculations
        while (l < final_l && rightward_partner[l] <= i) {
            ++l;
        }
        
        // note: l is 1 past the position we want to check, but we also need to switch from data to DP index
        if (outside_window_argmax != -1 &&
            (left_adj_constraint_prefix_sum[outside_window_argmax] != left_adj_constraint_prefix_sum[l]
             || right_adj_constraint_prefix_sum[argmax_partner] != right_adj_constraint_prefix_sum[i])) {
            // there is a window in between that falls below the minimum difference, we have to start over
            // for the values outside the window
            outside_window_argmax = -1;
        }
        
        // move the current window to the right
        curr_window_length += std::get<0>(shared_subanchors[i - 1]);
        
        while (window_begin < shared_subanchors.size() && curr_window_length > window_length) {
            // we can trim the leftmost value from the interval and maintain window length
            curr_window_length -= std::get<0>(shared_subanchors[window_begin]);
            
            // remove it from the sparse query data structure as well (by returning it to a null value)
            auto it = tree.find(std::make_pair(excl_fractional_difference[window_begin], window_begin));
            tree.update(it, mininf);
            
            // find the nearest position whose right-adjusted window comes after here
            while (k < shared_subanchors.size() && leftward_partner[k] + 1 < it->first.second) {
                ++k;
            }
            
            if ((left_adj_constraint_prefix_sum[it->first.second] == left_adj_constraint_prefix_sum[l]
                 && right_adj_constraint_prefix_sum[k] == right_adj_constraint_prefix_sum[i]) &&
                (outside_window_argmax == -1 ||
                 dp[it->first.second].first - length_prefix_sum[it->first.second] > dp[outside_window_argmax].first - length_prefix_sum[outside_window_argmax])) {
                
                outside_window_argmax = it->first.second;
                argmax_partner = k;
            }
            
            ++window_begin;
        }
        
        // find max if i is not included
        dp[i].first = std::max(dp[i - 1].first, dp[i - 1].second);
        
        // case 1: the segment is shorter than the window length
        
        auto it = tree.range_max(std::pair<double, size_t>(fractional_difference[i], 0), inf);
        
        if (it != tree.end() && it->second != mininf) {
            // there's a valid interval within the window
            backpointer[i] = it->first.second;
            dp[i].second = length_prefix_sum[i] + it->second - min_length;
            if (dp[i].second > dp[tb_idx].second) {
                tb_idx = i;
                if (debug) {
                    std::cerr << "this is the new traceback opt with score " << dp[i].second << "\n";
                }
            }
        }
        
        // case 2: the segment is longer than the window length
        
        if (outside_window_argmax != -1) {
            double outside_window_length = dp[outside_window_argmax].first + length_prefix_sum[i] - length_prefix_sum[outside_window_argmax] - min_length;
            if (debug) {
                std::cerr << "outside window opt score " << outside_window_length << " occurs at " << outside_window_argmax << '\n';
            }
            if (outside_window_length > dp[i].second) {
                dp[i].second = outside_window_length;
                backpointer[i] = outside_window_argmax;
                if (debug) {
                    std::cerr << "choose outside-window score as the dp value\n";
                }
            }
        }
        
        if (dp[i].second > dp[tb_idx].second) {
            tb_idx = i;
            if (debug) {
                std::cerr << "found new opt at " << tb_idx << " with score " << dp[i].second << '\n';
            }
        }
        
        if (i + 1 != dp.size()) {
            // enter into search tree
            tree.update(tree.find(std::pair<double, size_t>(excl_fractional_difference[i], i)),
                        dp[i].first - excl_length_prefix_sum[i]);
        }
        
    }
    
    return traceback(dp, backpointer, tb_idx);
}

}
