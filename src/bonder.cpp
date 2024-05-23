#include "centrolign/bonder.hpp"

#include <cassert>
#include <limits>
#include <iomanip>
#include <cmath>

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
            
            if ((end < 0 || end >= partner.size()) && curr_window_length < window_length) {
                // we don't have a full window anymore, so the previous full interval applies to the
                // rest of the vector
                partner[i] = end;
                if (i - incr >= 0 && i - incr < meets_constraint.size()) {
                    meets_constraint[i] = meets_constraint[i - incr];
                }
                else {
                    // the entire vector fits in one window
                    meets_constraint[i] = (window_sec_score > min_opt_proportion * window_opt_score);
                }
            }
            else {
                partner[i] = end - incr;
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
            }
            
            if (break_intervening_windows && (i % 2 == 1) && partner[i] == i) {
                // the entire window is in between two shared segments, which we treat as a violation of the constraint
                meets_constraint[i] = false;
            }
            
            const auto& segment = (i % 2 == 1) ? intervening_segments[i / 2] : shared_subanchors[i / 2];
            curr_window_length -= std::get<0>(segment);
            window_opt_score -= std::get<1>(segment);
            window_sec_score -= std::get<2>(segment);
        }
    }
    
    if (debug) {
        std::cerr << std::setprecision(3);
        std::cerr << "window walking results\n";
        std::cerr << "i" << '\t' << "Cl" << '\t' << "Cr" << '\t' << "Pr" << '\t' << "Pl" << '\t' << "L" << '\t' << "OS" << '\t' << "SS" << '\n';
        for (size_t i = 0; i < meets_constraint_left_adj.size(); ++i) {
            const auto& r = (i % 2 == 0) ? shared_subanchors[i / 2] : intervening_segments[i / 2];
            std::cerr << i << '\t' << meets_constraint_left_adj[i] << '\t' << meets_constraint_right_adj[i] << '\t' << rightward_partner[i] << '\t' << leftward_partner[i] << '\t' << std::get<0>(r) << '\t' << std::get<1>(r) << '\t' << std::get<2>(r) << '\n';
        }
    }
    
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
            
            tree_data.emplace_back(std::pair<double, size_t>(excl_fractional_difference[i], i), mininf);
        }
        
        
        length_prefix_sum[i + 1] = excl_length_prefix_sum[i] + std::get<0>(shared_subanchors[i]);
        opt_pref_sum += std::get<1>(shared_subanchors[i]);
        sec_pref_sum += std::get<2>(shared_subanchors[i]);
        
        fractional_difference[i + 1] = min_opt_proportion * opt_pref_sum - sec_pref_sum;
    }
    
    // we a spacer for a silent intervening segment before the first, which simplifies indexing
    std::vector<int> left_adj_constraint_prefix_sum(joined_size + 1, 0);
    std::vector<int> right_adj_constraint_prefix_sum(joined_size + 1, 0);
    
    for (size_t i = 1; i < left_adj_constraint_prefix_sum.size(); ++i) {
        left_adj_constraint_prefix_sum[i] = left_adj_constraint_prefix_sum[i - 1] + int(!meets_constraint_left_adj[i - 1]);
        right_adj_constraint_prefix_sum[i] = right_adj_constraint_prefix_sum[i - 1] + int(!meets_constraint_right_adj[i - 1]);
    }
    
    if (debug) {
        std::cerr << "prefix data structures\n";
        std::cerr << "i" << '\t' << "LP" << '\t' << "LPx" << '\t' << "FD" << '\t' << "FDx" << '\n';
        for (size_t i = 0; i < shared_subanchors.size() + 1; ++i) {
            std::cerr << i << '\t' << length_prefix_sum[i] << '\t' << excl_length_prefix_sum[i] << '\t' << fractional_difference[i] << '\t' << excl_fractional_difference[i] << '\n';//'\t' << left_adj_constraint_prefix_sum[i] << '\t' << excl_left_adj_constraint_prefix_sum[i] << '\t' << right_adj_constraint_prefix_sum[i] << '\t' << excl_right_adj_constraint_prefix_sum[i] << '\n';
        }
        std::cerr << "constraint prefix structures\n";
        std::cerr << "i" << '\t' << "LC" << '\t' << "RC" << '\n';
        for (size_t i = 0; i < left_adj_constraint_prefix_sum.size(); ++i) {
            std::cerr << i << '\t' << left_adj_constraint_prefix_sum[i] << '\t' << right_adj_constraint_prefix_sum[i] << '\n';
        }
    }
    
    MaxSearchTree<std::pair<double, size_t>, double> tree(tree_data);
    
    // records of (excluded, included)
    std::vector<std::pair<double, double>> dp(length_prefix_sum.size(), std::make_pair(mininf, mininf));
    dp.front().first = 0.0;
    dp.front().second = 0.0;
    std::vector<size_t> backpointer(dp.size(), -1);
    
    size_t tb_idx = 0;
    
    // for maintaining the shorter-than-window sparse query structure
    int64_t window_begin = 0;
    double curr_window_length = 0.0;
    
    // these start as null until an index falls out of the active window
    size_t outside_window_argmax = -1;
    size_t argmax_partner = -1;
    // index of first right-adjusted window to fully the right of the points outside the active window
    size_t k = 0;
    // index of the first left-adjusted window that extends beyond the active window
    int64_t l = 0;
    
    // figure out where we need to stop iterating l to get to the last full window, but no further
    size_t final_l = joined_size;
    {
        double tail_length = 0.0;
        while (final_l != 0 && tail_length < window_length) {
            if (final_l % 2 == 1) {
                tail_length += std::get<0>(shared_subanchors[final_l / 2]);
            }
            else {
                tail_length += std::get<0>(intervening_segments[final_l / 2 - 1]);
            }
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
        while (l < final_l && rightward_partner[l + 1] <= 2 * (i - 1)){// && rightward_partner[l + 1] <= 2 * i - 1) {
            if (debug) {
                std::cerr << "advance l from " << l << " to " << l + 1 << " with rightward partner " << rightward_partner[l + 1] << ", which is covered by " << 2 * (i - 1) << '\n';
            }
            ++l;
        }
        
        if (debug && outside_window_argmax != -1) {
            std::cerr << "checking if outside window argmax, current value " << outside_window_argmax << " is interrupted by windows left adj? " << (left_adj_constraint_prefix_sum[2 * outside_window_argmax] != left_adj_constraint_prefix_sum[l + 1]) << " from 2*argmax " << (2 * outside_window_argmax) << " and l incr " << (l + 1) << ", right adj? " << (right_adj_constraint_prefix_sum[argmax_partner] != right_adj_constraint_prefix_sum[2 * i - 1]) << " from argmax partner " << argmax_partner << " and i shifted " << (2 * i - 1) << '\n';
        }
        
        // note: l is 1 past the position we want to check, but we also need to switch from data to DP index
        if (outside_window_argmax != -1 &&
            (left_adj_constraint_prefix_sum[2 * outside_window_argmax] != left_adj_constraint_prefix_sum[l + 1]
             || right_adj_constraint_prefix_sum[argmax_partner] != right_adj_constraint_prefix_sum[2 * i - 1])) {
            // there is a window in between that falls below the minimum difference, we have to start over
            // for the values outside the window
            if (debug) {
                std::cerr << "setting argmax to null\n";
            }
            outside_window_argmax = -1;
        }
        
        // move the current window to the right
        if (i > 1) {
            curr_window_length += std::get<0>(intervening_segments[i - 2]);
        }
        curr_window_length += std::get<0>(shared_subanchors[i - 1]);
        
        while (window_begin < shared_subanchors.size() && curr_window_length > window_length) {
            // we can trim the leftmost value from the interval and maintain window length
            curr_window_length -= std::get<0>(shared_subanchors[window_begin]);
            if (window_begin != intervening_segments.size()) {
                curr_window_length -= std::get<0>(intervening_segments[window_begin]);
            }
            
            if (debug) {
                std::cerr << "moving current window begin " << window_begin << " out of the shorter-than-window query structure\n";
            }
            
            // remove it from the sparse query data structure as well (by returning it to a null value)
            auto it = tree.find(std::make_pair(excl_fractional_difference[window_begin], window_begin));
            tree.update(it, mininf);
            
            // find the nearest position whose right-adjusted window comes after here
            while (k < leftward_partner.size() && leftward_partner[k] < 2 * window_begin) {
                if (debug) {
                    std::cerr << "advance k with leftward partner " << leftward_partner[k] << " (this is left of " << (2 * window_begin) << ") from " << k << " to " << k + 1 << '\n';
                }
                ++k;
            }
            
            if (debug) {
                std::cerr << "checking window begin " << (2 * window_begin) << " vs l (incr) " << (l + 1) << ": " << (left_adj_constraint_prefix_sum[2 * window_begin] == left_adj_constraint_prefix_sum[l + 1]) << " and k " << k << " vs i " << (2 * i - 1) << ": " << (right_adj_constraint_prefix_sum[k] == right_adj_constraint_prefix_sum[2 * i - 1]) << " and query value " << dp[window_begin].first - excl_length_prefix_sum[window_begin] << '\n';
            }
            
            if ((left_adj_constraint_prefix_sum[2 * window_begin] == left_adj_constraint_prefix_sum[l + 1]
                 && right_adj_constraint_prefix_sum[k] == right_adj_constraint_prefix_sum[2 * i - 1]) &&
                (outside_window_argmax == -1 ||
                 dp[window_begin].first - excl_length_prefix_sum[window_begin] > dp[outside_window_argmax].first - excl_length_prefix_sum[outside_window_argmax])) {
                if (debug) {
                    std::cerr << "index " << window_begin << " meets constraints and beats current outside window argmax of " << outside_window_argmax << ", setting it as the new argmax\n";
                }
                outside_window_argmax = window_begin;
                argmax_partner = k;
            }
            
            ++window_begin;
        }
        
        if (debug) {
            std::cerr << "tree state:\n";
            for (const auto& r : tree) {
                if (r.second != mininf) {
                    std::cerr << '\t' << r.first.second << '\t' << r.first.first << '\t' << r.second << '\n';
                }
            }
        }
        
        // find max if i is not included
        dp[i].first = std::max(dp[i - 1].first, dp[i - 1].second);
        if (debug) {
            std::cerr << "excluded DP val is " << dp[i].first << '\n';
        }
        
        // case 1: the segment is shorter than the window length
        
        auto it = tree.range_max(std::pair<double, size_t>(fractional_difference[i], 0), inf);
        
        if (it != tree.end() && it->second != mininf) {
            // there's a valid interval within the window
            backpointer[i] = it->first.second;
            dp[i].second = length_prefix_sum[i] + it->second - min_length;
            if (debug) {
                std::cerr << "inside window opt is at " << it->first.second << " with score " << dp[i].second << '\n';
            }
        }
        
        // case 2: the segment is longer than the window length
        
        if (outside_window_argmax != -1) {
            double outside_window_length = dp[outside_window_argmax].first + length_prefix_sum[i] - excl_length_prefix_sum[outside_window_argmax] - min_length;
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
    
    if (debug) {
        std::cerr << "final DP state\n";
        std::cerr << "i" << '\t' << "X" << '\t' << "I" << '\t' << "BP" << '\n';
        for (size_t i = 0; i < dp.size(); ++i) {
            std::cerr << i << '\t' << dp[i].first << '\t' << dp[i].second << '\t' << backpointer[i] << '\n';
        }
    }
    
    return traceback(dp, backpointer, tb_idx);
}

std::vector<std::pair<size_t, size_t>>
Bonder::longest_deviation_constrained_partition(const std::vector<std::tuple<double, double, double>>& shared_subanchors,
                                                const std::vector<std::tuple<double, double, double>>& intervening_segments,
                                                const std::vector<std::pair<int64_t, int64_t>>& deviation) const {
    
    assert(intervening_segments.size() == shared_subanchors.size() - 1);
    
    static const double mininf = std::numeric_limits<double>::lowest();
    
    // records of (excluded, included)
    std::vector<std::pair<double, double>> dp(shared_subanchors.size() + 1, std::make_pair(mininf, mininf));
    dp.front().first = 0.0;
    dp.front().second = 0.0;
    std::vector<size_t> backpointer(dp.size(), -1);
    
    size_t tb_idx = 0;
    
    for (size_t i = 1; i < dp.size(); ++i) {
        
        dp[i].first = std::max(dp[i - 1].first, dp[i - 1].second);
        
        double running_length = 0.0;
        double running_opt_score = 0.0;
        double running_sec_score = 0.0;
        int64_t running_opt_dev = 0;
        int64_t running_sec_dev = 0;
        for (size_t j = i - 1; j < dp.size(); --j) {
            const auto& shared = shared_subanchors[j];
            running_length += std::get<0>(shared);
            running_opt_score += std::get<1>(shared);
            running_sec_score += std::get<2>(shared);
            if (j + 1 != i) {
                const auto& between = intervening_segments[j];
                running_length += std::get<0>(between);
                running_opt_score += std::get<1>(between);
                running_sec_score += std::get<2>(between);
                const auto& dev = deviation[j];
                running_opt_dev += dev.first;
                running_sec_dev += dev.second;
            }
            
            if (running_sec_score >= min_opt_proportion * running_opt_score &&
                abs(running_opt_dev - running_sec_dev) <= sqrt(running_length) * deviation_drift_factor) {
                
                double score = dp[j].first + running_length - min_length;
                
                if (score > dp[i].second) {
                    dp[i].second = score;
                    backpointer[i] = j;
                }
            }
        }
        
        if (dp[i].second > dp[tb_idx].second) {
            tb_idx = i;
        }
    }
    
    return traceback(dp, backpointer, tb_idx);
}


}
