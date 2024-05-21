#ifndef centrolign_partitioner_hpp
#define centrolign_partitioner_hpp

#include <vector>

#include "centrolign/anchorer.hpp"
#include "centrolign/score_function.hpp"
#include "centrolign/utility.hpp"

namespace centrolign {

struct anchor_t;

/*
 * Class to identify good segments within an anchor chain
 */
class Partitioner : public Extractor {
public:
    Partitioner(const ScoreFunction& score_function) : score_function(&score_function) {}
    Partitioner() = default;
    ~Partitioner() = default;
    
    
    template<class BGraph, class XMerge>
    std::vector<std::vector<anchor_t>> partition_anchors(std::vector<anchor_t>& anchor_chain,
                                                         const BGraph& graph1, const BGraph& graph2,
                                                         const SentinelTableau& tableau1, const SentinelTableau& tableau2,
                                                         const XMerge& xmerge1, const XMerge& xmerge2) const;
    
    // different constraint algorithms
    enum ConstraintMethod {Null, Unconstrained, MinAverage, MinWindowAverage};
    
    // which constraint to apply to to segments
    ConstraintMethod constraint_method = MinWindowAverage;
    // apply a minimum score to include a segment in the partition
    bool use_min_segment_score = true;
    
    // the minimium unscaled score for a segment
    double minimum_segment_score = 15000.0;
    // the minimum unscaled basewise average score for a segment
    double minimum_segment_average = 0.1;
    // window length used for min window average
    double window_length = 10000.0;
    // the parameter of a Holder generalized mean used to measure distance between anchors
    double generalized_length_mean = -0.5;
    
protected:
    
    const ScoreFunction* const score_function = nullptr;
    
    template<class T>
    std::vector<std::pair<size_t, size_t>> maximum_weight_partition(const std::vector<T>& data) const;
    
    template<class T>
    std::vector<std::pair<size_t, size_t>> average_constrained_partition(const std::vector<std::pair<T, T>>& data) const;
    
    template<class T>
    std::vector<std::pair<size_t, size_t>> window_average_constrained_partition(const std::vector<std::pair<T, T>>& data) const;
    
    template<class T>
    std::vector<std::pair<size_t, size_t>> traceback(const std::vector<std::pair<T, T>>& dp,
                                                     const std::vector<size_t>& backpointer,
                                                     size_t tb_idx) const;
    
};



/*
 * Template implementations
 */


template<class BGraph, class XMerge>
std::vector<std::vector<anchor_t>> Partitioner::partition_anchors(std::vector<anchor_t>& anchor_chain,
                                                                  const BGraph& graph1, const BGraph& graph2,
                                                                  const SentinelTableau& tableau1, const SentinelTableau& tableau2,
                                                                  const XMerge& xmerge1, const XMerge& xmerge2) const {
    
    std::vector<std::pair<size_t, size_t>> partition;
    // count how many matches we used from each set
    std::vector<size_t> num_anchors_from_set;
    for (const auto& anchor : anchor_chain) {
        while (num_anchors_from_set.size() <= anchor.match_set) {
            num_anchors_from_set.push_back(0);
        }
        ++num_anchors_from_set[anchor.match_set];
    }
    // we reduce the count penalty for match sets that were used multiple times in this chain (this helps
    // with recent, highly-identical duplications)
    auto anchor_score = [&](const anchor_t& anchor) -> double {
        return score_function->anchor_weight(anchor.count1 - num_anchors_from_set[anchor.match_set] + 1,
                                             anchor.count2 - num_anchors_from_set[anchor.match_set] + 1,
                                             anchor.walk1.size());
    };
    
    if (constraint_method == Null) {
        // put the whole chain into the partition in one block
        partition.emplace_back(0, anchor_chain.size());
    }
    else if (constraint_method == Unconstrained) {
        // there is no active average constraint, use the simple algorithm
        
        std::vector<double> partition_data;
        partition_data.reserve(anchor_chain.size());
        for (const auto& anchor : anchor_chain) {
            partition_data.push_back(anchor_score(anchor));
        }
        
        partition = std::move(maximum_weight_partition(partition_data));
    }
    else {
        
        auto graphs_between = extract_graphs_between(anchor_chain, graph1, graph2, tableau1, tableau2,
                                                     xmerge1, xmerge2);
        
        // construct a vector of anchors and linking regions to do the partition on
        std::vector<std::pair<double, double>> partition_data(anchor_chain.size() + graphs_between.size());
        for (size_t i = 0; i < partition_data.size(); ++i) {
            if (i % 2 == 0) {
                // this is a segment corresponding to a gap
                const auto& graph_pair = graphs_between[i / 2];
                std::vector<double> graph_sizes;
                
                // compute the minimum distance across either of the stitch graphs
                // TODO: what if instead a minimum-biased mean, like harmonic?
                for (auto subgraph_ptr : {&graph_pair.first, &graph_pair.second}) {
                    const auto& subgraph = *subgraph_ptr;
                    if (subgraph.subgraph.node_size() == 0) {
                        graph_sizes.push_back(0.00001);
                    }
                    else {
                        graph_sizes.push_back(source_sink_minmax(subgraph).first + 1);
                    }
                }
                
                partition_data[i].first = 0.0;
                // a measure of central tendency that's min-biased, but can rise to 2^4 * min
                partition_data[i].second = generalized_mean(graph_sizes.begin(), graph_sizes.end(),
                                                            generalized_length_mean);
            }
            else {
                const auto& anchor = anchor_chain[i / 2];
                partition_data[i].first = anchor_score(anchor);
                partition_data[i].second = anchor.walk1.size();
            }
        }
        
        // execute the partition algorithm
        if (constraint_method == MinAverage) {
            partition = std::move(average_constrained_partition(partition_data));
        }
        else if (constraint_method == MinWindowAverage) {
            partition = std::move(window_average_constrained_partition(partition_data));
        }
        else {
            throw std::runtime_error("Unrecognized partition constraint algorithm " + std::to_string((int) constraint_method));
        }
        
        
        // convert into intervals of only anchors
        for (auto& interval : partition) {
            interval.first /= 2;
            interval.second = (interval.second + 1) / 2; // round up to preserve past-the-last-ness
        }
    }
    
    std::vector<std::vector<anchor_t>> partition_segments;
    partition_segments.reserve(partition.size());
    
    for (const auto& interval : partition) {
        partition_segments.emplace_back();
        auto& segment = partition_segments.back();
        for (size_t i = interval.first; i < interval.second; ++i) {
            segment.emplace_back(std::move(anchor_chain[i]));
        }
    }
    
    static const bool instrument = false;
    if (instrument) {
        logging::log(logging::Debug, "Adjusted partitioning params: min score = " + std::to_string(score_function->score_scale * minimum_segment_score) + ", min average = " + std::to_string(score_function->score_scale * minimum_segment_average) + ", window length = " + std::to_string(window_length));
        for (size_t i = 0; i < partition.size(); ++i) {
            auto p = partition[i];
            auto& segment = partition_segments[i];
            double weight = 0.0;
            for (const auto& anchor : segment) {
                weight += anchor.score;
            }
            auto& first = segment.front();
            auto& last = segment.back();
            std::cerr << '|' << '\t' << p.first << '\t' << p.second << '\t' << weight << '\t' << first.walk1.front() << '\t' << last.walk1.back() << '\t' << first.walk2.front() << '\t' << last.walk2.back() << '\t' << (last.walk1.back() - first.walk1.front()) << '\t' << (last.walk2.back() - first.walk2.front());
            if (i != 0) {
                auto& prev_segment = partition_segments[i - 1];
                std::cerr << '\t' << (first.walk1.front() - prev_segment.back().walk1.back()) << '\t' << (first.walk2.front() - prev_segment.back().walk2.back());
            }
            else {
                std::cerr << '\t' << first.walk1.front() << '\t' << first.walk2.front();
            }
            std::cerr << '\n';
        }
    }
    
    return partition_segments;
}

template<class T>
std::vector<std::pair<size_t, size_t>> Partitioner::traceback(const std::vector<std::pair<T, T>>& dp,
                                                              const std::vector<size_t>& backpointer, size_t tb_idx) const {
    
    std::vector<std::pair<size_t, size_t>> partition;
    
    bool in_interval = true;
    while (tb_idx > 0) {
        if (in_interval) {
            size_t prev = backpointer[tb_idx];
            partition.emplace_back(prev, tb_idx);
            tb_idx = prev;
            in_interval = false;
        }
        else {
            in_interval = (dp[tb_idx].first == dp[tb_idx - 1].second);
            --tb_idx;
        }
    }
    
    // put in forward order
    std::reverse(partition.begin(), partition.end());
    
    return partition;
}

template<class T>
std::vector<std::pair<size_t, size_t>> Partitioner::maximum_weight_partition(const std::vector<T>& data) const {
    
    static const T mininf = std::numeric_limits<T>::lowest();
    
    // adjust the parameters by the scale
    T min_score = minimum_segment_score * score_function->score_scale;
    
    std::vector<T> prefix_sum(data.size() + 1, 0);
    for (size_t i = 0; i < data.size(); ++i) {
        prefix_sum[i + 1] = prefix_sum[i] + data[i];
    }
    
    // records of (score excluded, score included)
    std::vector<std::pair<T, T>> dp(data.size() + 1, std::make_pair(mininf, mininf));
    std::vector<size_t> backpointer(dp.size(), -1);
    
    // boundary conditions
    dp[0].first = 0;
    dp[0].second = 0;
    size_t prefix_argmax = 0;
    size_t tb_idx = 0;
    
    // main DP loop
    for (size_t i = 1; i < dp.size(); ++i) {
        // score if excluded
        dp[i].first = std::max(dp[i - 1].first, dp[i - 1].second);
        
        // score if included
        dp[i].second = dp[prefix_argmax].first + prefix_sum[i] - prefix_sum[prefix_argmax] - min_score;
        backpointer[i] = prefix_argmax;
        
        // remember if this is the best end to a partial partition found so far
        if (dp[i].first - prefix_sum[i] > dp[prefix_argmax].first - prefix_sum[prefix_argmax]) {
            prefix_argmax = i;
        }
        if (dp[i].second > dp[tb_idx].second) {
            tb_idx = i;
        }
    }
    
    return traceback(dp, backpointer, tb_idx);
}

template<class T>
std::vector<std::pair<size_t, size_t>> Partitioner::average_constrained_partition(const std::vector<std::pair<T, T>>& data) const {
    
    static const bool debug = false;
    
    static const T mininf = std::numeric_limits<T>::lowest();
    
    // adjust the parameters by the scale
    T min_score = minimum_segment_score * score_function->score_scale;
    T min_average = minimum_segment_average * score_function->score_scale;
    
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
    
    // records of (score if excluded, score if included)
    std::vector<std::pair<T, T>> dp(data.size() + 1, std::make_pair(mininf, mininf));
    // back pointers for the score with the current item included
    std::vector<size_t> backpointer(dp.size(), -1);
    
    
    // initialize the sparse range max query
    std::vector<std::pair<std::pair<T, size_t>, T>> tree_data;
    tree_data.reserve(data.size() + 1);
    for (size_t i = 0; i < data.size(); ++i) {
        tree_data.emplace_back(std::make_pair(fractional_prefix_sum[i], i + 1), mininf);
    }
    // the boundary condition
    tree_data.emplace_back(std::pair<T, size_t>(0, 0), 0);
    dp.front().first = 0;
    dp.front().second = 0;
    
    MaxSearchTree<std::pair<T, size_t>, T> search_tree(tree_data);
    
    size_t opt_idx = 0;
    
    // main dynamic programming
    for (size_t i = 1; i < dp.size(); ++i) {
        // find max if i is not included
        dp[i].first = std::max(dp[i - 1].first, dp[i - 1].second);
        
        // find max if i is included in the partition
        auto max_it = search_tree.range_max(std::pair<T, size_t>(mininf, 0),
                                            std::pair<T, size_t>(fractional_prefix_sum[i - 1], -1));
        if (max_it != search_tree.end() && max_it->second != mininf) {
            // we've completed DP for a valid interval start
            dp[i].second = prefix_sum[i - 1] + max_it->second - min_score;
            backpointer[i] = max_it->first.second;
            if (dp[i].second > dp[opt_idx].second) {
                opt_idx = i;
            }
        }
        
        // enter the query values for intervals starting on next cell
        auto it = search_tree.find(std::make_pair(fractional_prefix_sum[i - 1], i));
        search_tree.update(it, dp[i].first - prefix_sum[i - 1]);
    }
    
    // do traceback
    
    return traceback(dp, backpointer, opt_idx);
}

template<class T>
std::vector<std::pair<size_t, size_t>> Partitioner::window_average_constrained_partition(const std::vector<std::pair<T, T>>& data) const {
    
    static const bool debug = false;
    
    static const T mininf = std::numeric_limits<T>::lowest();
    
    // adjust the parameters by the scale
    T min_score = minimum_segment_score * score_function->score_scale;
    T min_average = minimum_segment_average * score_function->score_scale;
    
    std::vector<bool> meets_constraint_left_adj(data.size());
    std::vector<bool> meets_constraint_right_adj(data.size());
    std::vector<int64_t> leftward_partner(data.size());
    std::vector<int64_t> rightward_partner(data.size());
    
    for (bool forward : {true, false}) {
        // the score and weight of the full interval at the current iteration
        T window_score = 0.0;
        T window_weight = 0.0;
        // the past-the-last position of the interval of this iteration
        int64_t end = forward ? 0 : data.size() - 1;
        int64_t incr = forward ? 1 : -1;
        auto& meets_constraint = forward ? meets_constraint_left_adj : meets_constraint_right_adj;
        auto& partner = forward ? rightward_partner : leftward_partner;
        for (int64_t i = end; i < data.size() && i >= 0; i += incr) {
            while (end < data.size() && end >= 0 && window_weight < window_length) {
                window_score += data[end].first;
                window_weight += data[end].second;
                end += incr;
            }
            partner[i] = end;
            
            if ((end < 0 || end >= data.size()) && window_weight < window_length) {
                // we don't have a full window anymore, so the previous full interval applies to the
                // rest of the vector
                //std::cerr << "i " << i << ", incr " << incr << ", len " << meets_constraint.size() << '\n';
                if (i - incr >= 0 && i - incr < data.size()) {
                    meets_constraint[i] = meets_constraint[i - incr];
                }
                else {
                    // the entire vector is in one window
                    meets_constraint[i] = (window_score >= min_average * window_weight);
                }
            }
            else {
                T final_score, final_weight;
                std::tie(final_score, final_weight) = data[end - incr];
                
                // this algebraically works out to weighting the final interval proportional to how much it overlaps
                // it's a bit tortured to avoid divisions, in case we ever apply this to an integer type
                meets_constraint[i] = (final_weight * window_score + (window_length - window_weight) * final_score >= final_weight * min_average * window_length);
            }
            
            window_score -= data[i].first;
            window_weight -= data[i].second;
        }
    }
    
    if (debug) {
        std::cerr << "data\n";
        for (size_t i = 0; i < data.size(); ++i) {
            std::cerr << i << '\t' << data[i].first << '\t' << data[i].second << '\n';
        }
        std::cerr << "left-adjusted windows:\n";
        for (size_t i = 0; i < data.size(); ++i) {
            std::cerr << i << '\t' << rightward_partner[i] << '\t' << meets_constraint_left_adj[i] << '\n';
        }
        std::cerr << "right-adjusted windows:\n";
        for (size_t i = 0; i < data.size(); ++i) {
            std::cerr << leftward_partner[i] << '\t' << i << '\t' << meets_constraint_right_adj[i] << '\n';
        }
    }
    
    // to compute summed score
    std::vector<T> prefix_sum(data.size() + 1);
    // to compute average constraint satisfaction
    std::vector<T> fractional_prefix_sum(data.size() + 1);
    // to determine where constraint failures occur
    std::vector<int> left_adj_constraint_prefix_sum(data.size() + 1);
    std::vector<int> right_adj_constraint_prefix_sum(data.size() + 1);
    for (size_t i = 0; i < data.size(); ++i) {
        prefix_sum[i + 1] = prefix_sum[i] + data[i].first;
        fractional_prefix_sum[i + 1] = fractional_prefix_sum[i] + data[i].first - data[i].second * min_average;
        left_adj_constraint_prefix_sum[i + 1] = left_adj_constraint_prefix_sum[i] + (int) !meets_constraint_left_adj[i];
        right_adj_constraint_prefix_sum[i + 1] = right_adj_constraint_prefix_sum[i] + (int) !meets_constraint_right_adj[i];
    }
    
    if (debug) {
        std::cerr << "constraint prefix sums:\n";
        for (size_t i = 0; i < left_adj_constraint_prefix_sum.size(); ++i) {
            std::cerr << i << '\t' << left_adj_constraint_prefix_sum[i] << '\t' << right_adj_constraint_prefix_sum[i] << '\n';
        }
    }
    
    
    // initialize the sparse range max query
    std::vector<std::pair<std::pair<T, size_t>, T>> tree_data;
    tree_data.reserve(fractional_prefix_sum.size() + 1);
    for (size_t i = 0; i < fractional_prefix_sum.size(); ++i) {
        tree_data.emplace_back(std::make_pair(fractional_prefix_sum[i], i), mininf);
    }
    tree_data.front().second = 0;
    // we want this for the first queries when the window is still populating, but we need to get rid of it too...
    MaxSearchTree<std::pair<T, size_t>, T> search_tree(tree_data);
    
    // records of (score if excluded, score if included)
    std::vector<std::pair<T, T>> dp(data.size() + 1, std::make_pair(mininf, mininf));
    // back pointers for the score with the current item included
    std::vector<size_t> backpointer(dp.size(), -1);
    // boundary condition
    dp.front().first = 0;
    dp.front().second = 0;
    
    size_t tb_idx = 0;
    
    // for maintaining the shorter-than-window sparse query structure
    size_t window_begin = 0;
    T window_weight = 0.0;

    // these start as null until an index falls out of the active window
    size_t outside_window_argmax = -1;
    size_t argmax_partner = -1;
    // index of first right-adjusted window to fully the right of the points outside the active window
    size_t k = 0;
    // index of the first left-adjusted window that extends beyond the active window
    size_t l = 0;
    
    // figure out where we need to stop iterating l to get to the last full window, but no further
    size_t final_l = data.size();
    {
        T tail_weight = 0.0;
        while (final_l != 0 && tail_weight + data[final_l - 1].second < window_length) {
            tail_weight += data[final_l - 1].second;
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
            if (debug) {
                std::cerr << "l value " << l << " has rightward partner " << rightward_partner[l] << ", advancing\n";
            }
            ++l;
        }
        
        // note: l is 1 past the position we want to check, but we also need to switch from data to DP index
        if (outside_window_argmax != -1 &&
            (left_adj_constraint_prefix_sum[outside_window_argmax] != left_adj_constraint_prefix_sum[l]
             || right_adj_constraint_prefix_sum[argmax_partner] != right_adj_constraint_prefix_sum[i])) {
            // there is a window in between that falls below the minimum average, we have to start over
            // for the values outside the window
            if (debug) {
                std::cerr << "encountered failing window between argmax " << outside_window_argmax << " (k " << argmax_partner << ") and " << i << "(l " << l << ")" << '\n';
            }
            outside_window_argmax = -1;
        }
        else if (debug && outside_window_argmax != -1) {
            std::cerr << "argmax " << outside_window_argmax << " with partner " << argmax_partner << " is valid in this iteration\n";
        }
        
        // move the current window to the right
        window_weight += data[i - 1].second;
        
        if (debug) {
            std::cerr << "add data point " << (i - 1) << " to window beginning at " << window_begin << " for total weight " << window_weight << " of window length " << window_length << '\n';
        }
        
        while (window_begin < data.size() && window_weight > window_length) {
            // we can trim the leftmost value from the interval and maintain window length
            window_weight -= data[window_begin].second;
            
            // remove it from the sparse query data structure as well (by returning it to a null key)
            auto it = search_tree.find(std::make_pair(fractional_prefix_sum[window_begin], window_begin));
            search_tree.update(it, mininf);
            
            if (debug) {
                std::cerr << "move index " << it->first.second << ", key (" << it->first.first << ',' << it->first.second << ") out of active window for total weight " << window_weight << '\n';
            }
            
            // find the nearest position whose right-adjusted window comes after here
            while (k < data.size() && leftward_partner[k] + 1 < it->first.second) {
                if (debug) {
                    std::cerr << "k value " << k << " has leftward partner " << leftward_partner[k] << ", advancing\n";
                }
                ++k;
            }
            
            if (debug) {
                std::cerr << "check constraint for indexes j " << it->first.second << ", k " << k << ", l " << l  << ", i " << i << ", sum values " << left_adj_constraint_prefix_sum[it->first.second] << ' ' << left_adj_constraint_prefix_sum[l] << ' ' << right_adj_constraint_prefix_sum[k] << ' ' << right_adj_constraint_prefix_sum[i] << ", dp value " << (dp[it->first.second].first - prefix_sum[it->first.second]) << " vs current argmax " << (outside_window_argmax == -1 ? std::string(".") : std::to_string(dp[outside_window_argmax].first - prefix_sum[outside_window_argmax])) << '\n';
            }
            
            if ((left_adj_constraint_prefix_sum[it->first.second] == left_adj_constraint_prefix_sum[l]
                 && right_adj_constraint_prefix_sum[k] == right_adj_constraint_prefix_sum[i]) &&
                (outside_window_argmax == -1 ||
                 dp[it->first.second].first - prefix_sum[it->first.second] > dp[outside_window_argmax].first - prefix_sum[outside_window_argmax])) {
                
                outside_window_argmax = it->first.second;
                argmax_partner = k;
                if (debug) {
                    std::cerr << "found new argmax at " << outside_window_argmax << '\n';
                }
            }
            
            ++window_begin;
        }
        
        // find max if i is not included
        dp[i].first = std::max(dp[i - 1].first, dp[i - 1].second);
        if (debug) {
            std::cerr << "excluded DP value set to " << dp[i].first << " from prev excluded " << dp[i - 1].first << ", included " << dp[i - 1].second << '\n';
        }
        
        // case 1: the segment is shorter than the window length
        
        if (debug) {
            std::cerr << "query key: " << fractional_prefix_sum[i] << '\n';
            std::cerr << "sparse query contents:\n";
            for (auto r : search_tree) {
                if (r.second != mininf) {
                    std::cerr << " (" << r.first.second << ',' << r.first.first << "):" << r.second;
                }
            }
            std::cerr << '\n';
        }
        
        auto max_it = search_tree.range_max(std::pair<T, size_t>(mininf, 0),
                                            std::pair<T, size_t>(fractional_prefix_sum[i], -1));
        if (max_it != search_tree.end() && max_it->second != mininf) {
            // there's a valid interval within the window
            dp[i].second = prefix_sum[i] + max_it->second - min_score;
            backpointer[i] = max_it->first.second;
            if (debug) {
                std::cerr << "update dp value to " << dp[i].second << " using inside-window query at " << max_it->first.second << '\n';
            }
        }
        
        // case 2: the segment is longer than the window length
        
        if (outside_window_argmax != -1) {
            T outside_window_score = dp[outside_window_argmax].first + prefix_sum[i] - prefix_sum[outside_window_argmax] - min_score;
            if (debug) {
                std::cerr << "outside window opt score " << outside_window_score << " occurs at " << outside_window_argmax << '\n';
            }
            if (outside_window_score > dp[i].second) {
                dp[i].second = outside_window_score;
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
        
        // enter into search tree
        auto it = search_tree.find(std::make_pair(fractional_prefix_sum[i], i));
        search_tree.update(it, dp[i].first - prefix_sum[i]);
    }
    
    if (debug) {
        std::cerr << "final DP state:\n";
        for (size_t i = 0; i < dp.size(); ++i) {
            std::cerr << i << '\t';
            if (dp[i].first == mininf) {
                std::cerr << '.';
            }
            else {
                std::cerr << dp[i].first;
            }
            std::cerr << '\t';
            if (dp[i].second == mininf) {
                std::cerr << '.';
            }
            else {
                std::cerr << dp[i].second;
            }
            std::cerr << '\t';
            if (backpointer[i] == -1) {
                std::cerr << '.';
            }
            else {
                std::cerr << backpointer[i];
            }
            std::cerr << '\n';
        }
    }
    
    return traceback(dp, backpointer, tb_idx);
}

}

#endif /* centrolign_partitioner_hpp */
