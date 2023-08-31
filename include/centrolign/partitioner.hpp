#ifndef centrolign_partitioner_hpp
#define centrolign_partitioner_hpp

#include <vector>

#include "centrolign/anchorer.hpp"

namespace centrolign {

struct anchor_t;

/*
 * Class to identify good segments within an anchor chain
 */
class Partitioner : public Extractor {
public:
    Partitioner() = default;
    ~Partitioner() = default;
    
    
    template<class BGraph, class XMerge>
    std::vector<std::vector<anchor_t>> partition_anchors(std::vector<anchor_t>& anchor_chain,
                                                         const BGraph& graph1, const BGraph& graph2,
                                                         const SentinelTableau& tableau1, const SentinelTableau& tableau2,
                                                         const XMerge& xmerge1, const XMerge& xmerge2) const;
    
    // different constraint algorithms
    enum ConstraintMethod {Unconstrained, MinAverage, MinWindowAverage};
    
    // which constraint to apply to to segments
    ConstraintMethod constraint_method = MinAverage;
    // apply a minimum score to include a segment in the partition
    bool use_min_segment_score = true;
    
    // records the intrinsic scale of the scoring function on these sequences
    double score_scale = 0.280039; // ~ chr12 value
    
    // the minimium unscaled score for a segment
    double minimum_segment_score = 15000.0;
    // the minimum unscaled basewise average score for a segment
    double minimum_segment_average = 0.5;
    // window length used for min window average
    double window_length = 25000.0;
protected:
    
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
    
    if (constraint_method == Unconstrained) {
        // there is no active average constraint, use the simple algorithm
        
        std::vector<double> partition_data;
        partition_data.reserve(anchor_chain.size());
        for (const auto& anchor : anchor_chain) {
            partition_data.push_back(anchor.score);
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
                size_t graph_size = std::numeric_limits<size_t>::max();
                
                // compute the minimum distance across either of the stitch graphs
                // TODO: what if instead a minimum-biased mean, like harmonic?
                for (auto subgraph_ptr : {&graph_pair.first, &graph_pair.second}) {
                    const auto& subgraph = *subgraph_ptr;
                    if (subgraph.subgraph.node_size() == 0) {
                        // we still count empty graphs as having size 1 to avoid divide by 0 issues
                        graph_size = 1;
                    }
                    else {
                        graph_size = std::min<size_t>(graph_size, source_sink_minmax(subgraph).first);
                    }
                }
                
                partition_data[i].first = 0.0;
                partition_data[i].second = graph_size;
            }
            else {
                const auto& anchor = anchor_chain[i / 2];
                partition_data[i].first = anchor.score;
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
    
    static const bool instrument = true;
    if (instrument) {
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
                std::cerr << '\t' << 0 << '\t' << 0;
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
    T min_score = minimum_segment_score * score_scale;
    
    std::vector<T> prefix_sum(data.size() + 1, 0);
    for (size_t i = 0; i < data.size(); ++i) {
        prefix_sum[i + 1] = prefix_sum[i] + data[i];
    }
    
    // records of (score excluded, score included)
    std::vector<std::pair<T, T>> dp(data.size() + 1, std::make_pair(mininf, mininf));
    std::vector<size_t> backpointer(dp.size(), -1);
    
    // boundary conditions
    dp[0].first = 0;
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
    T min_score = minimum_segment_score * score_scale;
    T min_average = minimum_segment_average * score_scale;
    
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
    
    static const T mininf = std::numeric_limits<T>::lowest();
    
    // adjust the parameters by the scale
    T min_score = minimum_segment_score * score_scale;
    T min_average = minimum_segment_average * score_scale;
    
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
            while (end < data.size() && window_weight < window_length) {
                window_score += data[end].first;
                window_weight += data[end].second;
                end += incr;
            }
            partner[i] = end;
            
            T final_score, final_weight;
            std::tie(final_score, final_weight) = data[end - incr];
            
            // this algebraically works out to weighting the final interval proportional to how much it overlaps
            // it's a bit tortured to avoid divisions, in case we ever apply this to an integer type
            meets_constraint[i] = (final_weight * window_score + (window_length - window_weight) * final_score >= final_weight * min_average);
            
            window_score -= data[i].first;
            window_weight -= data[i].second;
        }
    }
    
    // to compute summed score
    std::vector<T> prefix_sum(data.size() + 1);
    // to compute average constraint satisfaction
    std::vector<T> fractional_prefix_sum(data.size() + 1);
    std::vector<int> left_adj_constraint_prefix_sum(data.size() + 1);
    std::vector<int> right_adj_constraint_prefix_sum(data.size() + 1);
    for (size_t i = 0; i < data.size(); ++i) {
        prefix_sum[i + 1] = prefix_sum[i] + data[i].first;
        fractional_prefix_sum[i + 1] = fractional_prefix_sum[i] + data[i].first - data[i].second * min_average;
        left_adj_constraint_prefix_sum[i + 1] = left_adj_constraint_prefix_sum[i] + (int) !meets_constraint_left_adj[i];
        right_adj_constraint_prefix_sum[i + 1] = right_adj_constraint_prefix_sum[i] + (int) !meets_constraint_right_adj[i];
    }
    
    // initialize the sparse range max query
    std::vector<std::pair<std::pair<T, size_t>, T>> tree_data;
    tree_data.reserve(fractional_prefix_sum.size() + 1);
    for (size_t i = 0; i < fractional_prefix_sum.size(); ++i) {
        tree_data.emplace_back(std::make_pair(fractional_prefix_sum[i], i), mininf);
    }
    // we want this for the first queries when the window is still populating, but we need to get rid of it too...
    MaxSearchTree<std::pair<T, size_t>, T> search_tree(tree_data);
    
    // records of (score if excluded, score if included)
    std::vector<std::pair<T, T>> dp(data.size() + 1, std::make_pair(mininf, mininf));
    // back pointers for the score with the current item included
    std::vector<size_t> backpointer(dp.size(), -1);
    // boundary condition
    dp.front().first = 0;
    
    size_t tb_idx = 0;
    
    // for maintaining the shorter-than-window sparse query structure
    size_t window_begin = 0;
    T window_weight = 0.0;

    size_t outside_window_argmax = 0;
    size_t argmax_partner = 0;
    size_t k = 0;
    size_t l = 0;
    
    for (size_t i = 1; i < dp.size(); ++i) {
                
        // find the first position whose left-adjusted window finishes after i's right side
        // note: the partner is past-the-last, and i is 1 greater than the data index
        while (l < data.size() && rightward_partner[l] < i) {
            ++l;
        }
        
        // note: l is 1 past the position we want to check, but we also need to switch from data to DP index
        if (outside_window_argmax != -1 &&
            (left_adj_constraint_prefix_sum[outside_window_argmax] != left_adj_constraint_prefix_sum[l]
             || right_adj_constraint_prefix_sum[argmax_partner] != right_adj_constraint_prefix_sum[i])) {
            // there is a window in between that falls below the minimum average, we have to start over
            // for the values outside the window
            outside_window_argmax = -1;
        }
        
        // move the current window to the right
        window_weight += data[i - 1].second;
        T front_weight = window_begin == 0 ? 0.0 : data[window_begin].second;
        while (window_weight >= window_length + front_weight) {
            // we can trim the leftmost value from the interval and maintain window length
            window_weight -= front_weight;
            
            // remove it from the sparse query data structure as well (by returning it to a null key)
            auto it = search_tree.find(std::make_pair(fractional_prefix_sum[window_begin], window_begin));
            search_tree.update(it, mininf);
            
            // find the nearest position whose right-adjusted window comes after here
            // note: we need to offset the "past-the-first" leftward partner to get the first index actually included
            // in k's window, and we need to offset the tree key by 1 to go from DP indexes to data index
            while (k < data.size() && leftward_partner[k] + 1 < (int64_t) it->first.second - 1) {
                ++k;
            }
            
            if (outside_window_argmax == -1 ||
                dp[it->first.second].first - prefix_sum[it->first.second] > dp[outside_window_argmax].first - prefix_sum[outside_window_argmax]) {
                
                outside_window_argmax = it->first.second;
                argmax_partner = k;
            }
            
            ++window_begin;
            front_weight = data[window_begin - 1].second;
        }
        
        // find max if i is not included
        dp[i].first = std::max(dp[i - 1].first, dp[i - 1].second);
        
        // case 1: the segment is shorter than the window length
        
        auto max_it = search_tree.range_max(std::pair<T, size_t>(mininf, 0),
                                            std::pair<T, size_t>(fractional_prefix_sum[i], -1));
        if (max_it != search_tree.end() && max_it->second != mininf) {
            // there's a valid interval within the window a valid interval start
            dp[i].second = prefix_sum[i] + max_it->second - min_score;
            backpointer[i] = max_it->first.second;
        }
        
        // case 2: the segment is longer than the window length
        
        if (outside_window_argmax != -1) {
            T outside_window_score = dp[outside_window_argmax].first + prefix_sum[i] - prefix_sum[outside_window_argmax] - min_score;
            if (outside_window_score > dp[i].second) {
                dp[i].second = outside_window_score;
                backpointer[i] = outside_window_argmax;
            }
        }
        
        if (dp[i].second > dp[tb_idx].second) {
            tb_idx = i;
        }
        
        // enter into search tree
        auto it = search_tree.find(std::make_pair(fractional_prefix_sum[i], i));
        search_tree.update(it, dp[i].first - prefix_sum[i]);
    }
    
    return traceback(dp, backpointer, tb_idx);
}

}

#endif /* centrolign_partitioner_hpp */
