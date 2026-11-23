#include "centrolign/anchorer.hpp"

#include <algorithm>
#include <sstream>
#include <iomanip>

#include "centrolign/topological_order.hpp"
#include "centrolign/minmax_distance.hpp"

namespace centrolign {

using namespace std;

const bool Anchorer::debug_anchorer = false;

pair<int64_t, int64_t> Extractor::source_sink_minmax(const SubGraphInfo& extraction) {
    
    auto dists = minmax_distance(extraction.subgraph, &extraction.sources);
    pair<int64_t, int64_t> mm(numeric_limits<int64_t>::max(), -1);
    for (auto node_id : extraction.sinks) {
        mm.first = min(mm.first, dists[node_id].first);
        mm.second = max(mm.second, dists[node_id].second);
    }
    return mm;
}

SparseDPProfiler::SparseDPProfiler(const std::vector<uint64_t>& partition_begin, const std::vector<uint64_t>& partition_end) :
    partition_begin(partition_begin), partition_end(partition_end), thread_stats(partition_begin.size()), signal_time(0) {
    class_batches.fill(0);
    class_nanos.fill(0);
}

void SparseDPProfiler::finish_batch(size_t batch_size, uint64_t signal_time, uint64_t own_end, uint64_t wait_end) {

    wait_nanos += wait_end - own_end;

    // the time lost to uneven load is how much longer the slowest thread took than the average
    uint64_t max_busy = 0;
    double mean_busy = 0.0;
    for (const auto& stats : thread_stats) {
        max_busy = std::max(max_busy, stats.last_busy_nanos);
        mean_busy += stats.last_busy_nanos;
    }
    mean_busy /= thread_stats.size();
    imbalance_nanos += max_busy - mean_busy;

    ++num_batches;
    num_queries += batch_size;

    // size classes: 1, 2-7, 8-63, 64+
    size_t size_class = batch_size <= 1 ? 0 : (batch_size < 8 ? 1 : (batch_size < 64 ? 2 : 3));
    ++class_batches[size_class];
    class_nanos[size_class] += wait_end - signal_time;
}

void SparseDPProfiler::log(const std::string& name, size_t num_nodes, uint64_t begin, uint64_t end) const {

    if (thread_stats.empty()) {
        // there were no chains to partition
        return;
    }

    auto secs = [](double nanos) {
        std::stringstream strm;
        strm << std::fixed << std::setprecision(2) << nanos / 1e9;
        return strm.str();
    };

    uint64_t min_chains = -1, max_chains = 0;
    for (size_t t = 0; t < partition_begin.size(); ++t) {
        min_chains = std::min(min_chains, partition_end[t] - partition_begin[t]);
        max_chains = std::max(max_chains, partition_end[t] - partition_begin[t]);
    }

    uint64_t total_update_nanos = 0, total_query_nanos = 0, total_barrier_nanos = 0, total_wake_nanos = 0;
    uint64_t total_updates = 0, total_query_pairs = 0;
    uint64_t min_busy = -1, max_busy = 0;
    double mean_busy = 0.0;
    for (const auto& stats : thread_stats) {
        uint64_t busy = stats.update_nanos + stats.query_nanos;
        min_busy = std::min(min_busy, busy);
        max_busy = std::max(max_busy, busy);
        mean_busy += busy;
        total_update_nanos += stats.update_nanos;
        total_query_nanos += stats.query_nanos;
        total_barrier_nanos += stats.barrier_nanos;
        total_wake_nanos += stats.wake_nanos;
        total_updates += stats.num_updates;
        total_query_pairs += stats.num_query_pairs;
    }
    mean_busy /= thread_stats.size();

    uint64_t total_nanos = end - begin;
    double mean_wake_micros = 0.0;
    if (thread_stats.size() > 1 && num_batches != 0) {
        mean_wake_micros = total_wake_nanos / (1000.0 * num_batches * (thread_stats.size() - 1));
    }
    double query_pair_nanos = total_query_pairs == 0 ? 0.0 : double(total_query_nanos) / total_query_pairs;
    int64_t huge_page_mem = huge_page_memory_usage();

    logging::log(logging::Debug, name + " DP profile: " + std::to_string(thread_stats.size()) + " threads over " + std::to_string(partition_end.back()) + " chains in graph2 (" + std::to_string(min_chains) + "-" + std::to_string(max_chains) + " per thread for updates), " + std::to_string(num_nodes) + " nodes, " + std::to_string(num_batches) + " synchronized batches, " + std::to_string(num_queries) + " queries, " + std::to_string(total_updates) + " tree updates, " + secs(total_nanos) + " s total" + (huge_page_mem >= 0 ? ", " + format_memory_usage(huge_page_mem) + " of the process in huge pages." : "."));
    logging::log(logging::Debug, name + " DP main thread: " + secs(build_nanos) + " s building batches, " + secs(thread_stats.front().update_nanos + thread_stats.front().query_nanos) + " s on its own updates and queries, " + secs(wait_nanos) + " s waiting for workers, " + secs(reduce_nanos) + " s reducing results.");
    logging::log(logging::Debug, name + " DP threads: busy " + secs(min_busy) + " s min, " + secs(mean_busy) + " s mean, " + secs(max_busy) + " s max (mean " + std::to_string(int(round(100.0 * mean_busy / std::max<uint64_t>(total_nanos, 1)))) + "% of total); " + secs(total_update_nanos) + " s updating and " + secs(total_query_nanos) + " s querying summed over threads; " + secs(total_barrier_nanos) + " s waiting between updates and queries summed over threads; " + secs(imbalance_nanos) + " s lost to uneven load; mean wake delay " + std::to_string(int(round(mean_wake_micros))) + " us; " + std::to_string(int(round(query_pair_nanos))) + " ns per query-chain pair.");
    std::string classes;
    std::array<std::string, num_size_classes> class_names{"1", "2-7", "8-63", "64+"};
    for (size_t i = 0; i < num_size_classes; ++i) {
        if (i != 0) {
            classes += ", ";
        }
        classes += class_names[i] + ": " + std::to_string(class_batches[i]) + " batches in " + secs(class_nanos[i]) + " s";
    }
    logging::log(logging::Debug, name + " DP batches by number of queries (time from signal to all threads finished): " + classes + ".");
}

std::vector<size_t> Extractor::get_logging_indexes(size_t size) {
    std::vector<size_t> logging_indexes;
    for (size_t i = 1; i < 10; ++i) {
        logging_indexes.push_back((size * i) / 10);
    }
    auto end = std::unique(logging_indexes.begin(), logging_indexes.end());
    logging_indexes.resize(end - logging_indexes.begin());
    return logging_indexes;
}

uint64_t Anchorer::AnchorGraph::add_node(size_t set, size_t idx1, size_t idx2, double weight,
                                         double initial_weight, double final_weight) {
    nodes.emplace_back(set, idx1, idx2, weight, initial_weight, final_weight);
    return nodes.size() - 1;
}

void Anchorer::AnchorGraph::add_edge(uint64_t from_id, uint64_t to_id, double weight) {
    nodes[from_id].edges.push_back(to_id);
    nodes[from_id].edge_weights.push_back(weight);
    ++nodes[to_id].in_degree;
}

const vector<size_t>& Anchorer::AnchorGraph::next(uint64_t node_id) const {
    return nodes[node_id].edges;
}

size_t Anchorer::AnchorGraph::next_size(uint64_t node_id) const {
    return nodes[node_id].edges.size();
}

size_t Anchorer::AnchorGraph::previous_size(uint64_t node_id) const {
    return nodes[node_id].in_degree;
}

size_t Anchorer::AnchorGraph::node_size() const {
    return nodes.size();
}

tuple<size_t, size_t, size_t> Anchorer::AnchorGraph::label(uint64_t node_id) const {
    const auto& node = nodes[node_id];
    return make_tuple(node.set, node.idx1, node.idx2);
}

vector<uint64_t> Anchorer::AnchorGraph::heaviest_weight_path(double min_score) const {
    
    // set up DP structures
    std::vector<size_t> backpointer(nodes.size(), -1);
    std::vector<double> dp(nodes.size());
    for (size_t i = 0; i < nodes.size(); ++i) {
        dp[i] = nodes[i].initial_weight;
    }
    
    // do dynamic programming
    uint64_t max_id = -1;
    double max_weight = min_score;
    for (auto node_id : topological_order(*this)) {
        
        // add nodes own weight in
        const auto& node = nodes[node_id];
        auto& dp_here = dp[node_id];
        if (dp_here == numeric_limits<double>::lowest()) {
            continue;
        }
        dp_here += node.weight;
        
        if (debug_anchorer) {
            cerr << "DP at " << node_id << " set to " << dp_here << '\n';
        }
        
        if (node.final_weight != numeric_limits<double>::lowest() &&
            dp_here + node.final_weight > max_weight) {
            max_id = node_id;
            max_weight = dp_here + node.final_weight;
        }
        
        // propagate forward
        for (size_t i = 0; i < node.edges.size(); ++i) {
            auto next_id = node.edges[i];
            auto edge_weight = node.edge_weights[i];
            if (dp_here + edge_weight > dp[next_id]) {
                dp[next_id] = dp_here + edge_weight;
                backpointer[next_id] = node_id;
            }
        }
    }
    
    if (debug_anchorer) {
        cerr << "max DP occurs at " << max_id << " with DP value " << max_weight << '\n';
    }
    
    // traceback
    vector<uint64_t> traceback;
    if (max_id != -1) {
        traceback.push_back(max_id);
        while (backpointer[traceback.back()] != -1) {
            if (debug_anchorer) {
                cerr << "trace back to " << backpointer[traceback.back()] << '\n';
            }
            traceback.push_back(backpointer[traceback.back()]);
        }
        reverse(traceback.begin(), traceback.end());
    }
    
    if (debug_anchorer) {
        cerr << "heaviest weight path consists of " << traceback.size() << " anchors, which have total weight " << max_weight << '\n';
    }
    
    return traceback;
}


std::vector<size_t> Anchorer::assign_reanchor_budget(const std::vector<std::pair<SubGraphInfo, SubGraphInfo>>& fill_in_graphs) const {
    
    // the sum of the matrix sizes
    size_t total = 0;
    for (const auto& stitch_pair : fill_in_graphs) {
        total += (stitch_pair.first.subgraph.node_size() + 1) * (stitch_pair.second.subgraph.node_size() + 1);
    }
    
    std::vector<size_t> budget;
    budget.reserve(fill_in_graphs.size());
    for (const auto& stitch_pair : fill_in_graphs) {
        // proportionate to this matrix's contribution to the sum of matrix sizes
        size_t size = (stitch_pair.first.subgraph.node_size() + 1) * (stitch_pair.second.subgraph.node_size() + 1);
        budget.push_back(ceil(double(max_num_match_pairs) * double(size) / double(total)));
    }
    
    return budget;
}


void Anchorer::merge_fill_in_chains(std::vector<anchor_t>& anchors,
                                    const std::vector<std::vector<anchor_t>> fill_in_chains,
                                    const std::vector<std::pair<SubGraphInfo, SubGraphInfo>>& fill_in_graphs,
                                    const std::vector<std::vector<std::pair<size_t, std::pair<std::vector<size_t>, std::vector<size_t>>>>>& match_origin) const {
    
    std::vector<anchor_t> merged;
    
    assert(anchors.size() + 1 == fill_in_chains.size());
    
    for (size_t i = 0; i < fill_in_chains.size(); ++i) {
        int remember_idx;
        if (i != 0) {
            // copy the original anchor
            auto& anchor = anchors[i - 1];
            if (!merged.empty()) {
                // update the gap before the original anchor
                anchor.gap_before = merged.back().gap_after;
                anchor.gap_score_before = merged.back().gap_score_after;
            }
            merged.emplace_back(std::move(anchor));
        }
        
        const auto& back_trans1 = fill_in_graphs[i].first.back_translation;
        const auto& back_trans2 = fill_in_graphs[i].second.back_translation;
        
        const auto& to_merge = fill_in_chains[i];
        for (size_t j = 0; j < to_merge.size(); ++j) {
            const auto& anchor = to_merge[j];
            if (j == 0 && !merged.empty()) {
                // update the gap after the original anchor
                merged.back().gap_score_after = anchor.gap_score_before;
                merged.back().gap_after = anchor.gap_before;
            }
            // copy the anchor in this in-between stitching chain
            merged.emplace_back();
            auto& translated = merged.back();
            translated.walk1.reserve(anchor.walk1.size());
            translated.walk2.reserve(anchor.walk2.size());
            // the node IDs must be translated from subgraph to main graph
            for (auto node_id : anchor.walk1) {
                translated.walk1.push_back(back_trans1[node_id]);
            }
            for (auto node_id : anchor.walk2) {
                translated.walk2.push_back(back_trans2[node_id]);
            }
            // we preserve the original queried count from the match index
            translated.count1 = anchor.count1;
            translated.count2 = anchor.count2;
            const auto& origin_set = match_origin[i][anchor.match_set];
            
            // back translate the identity of this anchor
            translated.match_set = origin_set.first;
            translated.idx1 = origin_set.second.first[anchor.idx1];
            translated.idx2 = origin_set.second.second[anchor.idx2];
            translated.full_length = anchor.full_length;
            translated.score = anchor.score;
            translated.gap_before = anchor.gap_before;
            translated.gap_score_before = anchor.gap_score_before;
            translated.gap_after = anchor.gap_after;
            translated.gap_score_after = anchor.gap_score_after;
        }
    }
    
    anchors = std::move(merged);
}

}
