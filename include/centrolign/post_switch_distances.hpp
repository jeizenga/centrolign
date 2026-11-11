#ifndef centrolign_post_switch_distances_hpp
#define centrolign_post_switch_distances_hpp

#include <vector>
#include <cstdint>
#include <thread>
#include <algorithm>

#include "centrolign/utility.hpp"

namespace centrolign {

/*
 * The D arrays from Chandra, et al. 2021, which represent the distance walked
 * to reach the node after leaving a given path
 */
template<class DistVector>
class PostSwitchDistances {
public:
    
    template<class BGraph, class XMerge>
    PostSwitchDistances(const BGraph& graph, const XMerge& xmerge, uint64_t threads = 1) noexcept;
    PostSwitchDistances(PostSwitchDistances<DistVector>&& other) noexcept = default;
    PostSwitchDistances(const PostSwitchDistances<DistVector>& other) noexcept = default;
    PostSwitchDistances() = default;
    ~PostSwitchDistances() = default;
    
    PostSwitchDistances& operator=(const PostSwitchDistances<DistVector>& other) noexcept = default;
    PostSwitchDistances& operator=(PostSwitchDistances<DistVector>&& other) noexcept = default;
    
    inline size_t distance(uint64_t to_node_id, uint64_t from_path_id) const;
    
    inline size_t memory_size() const;
    
private:
    
    std::vector<DistVector> distances;
    
};


/*
 * Template implementations
 */

template<class DistVector>
template<class BGraph, class XMerge>
PostSwitchDistances<DistVector>::PostSwitchDistances(const BGraph& graph, const XMerge& xmerge, uint64_t threads) noexcept {
    
    // we reserve 0 as a sentinel for no distance and add 1 to all other values
    distances.reserve(xmerge.chain_size());
    for (uint64_t path_id = 0; path_id < xmerge.chain_size(); ++path_id) {
        distances.emplace_back(graph.node_size());
    }
    
    auto order = topological_order(graph);
    
    // the paths are independent, so each thread computes the rows of a contiguous range of them
    auto compute_rows = [&](uint64_t path_begin, uint64_t path_end) {
        // note: this DP is different from the paper because i have the predecessor on the same
        // path being the previous node rather than the node itself
        for (auto node_id : order) {
            for (uint64_t p = path_begin; p < path_end; ++p) {
                auto& path_row = distances[p];
                for (auto prev_id : graph.previous(node_id)) {
                    auto pred = xmerge.predecessor_index(node_id, p);
                    if (xmerge.index_on(prev_id, p) == pred) {
                        // switching paths you here immediately, no distance
                        path_row[node_id] = 1;
                        break;
                    }
                    else if (xmerge.predecessor_index(prev_id, p) == pred) {
                        size_t dist_thru = size_t(path_row[prev_id]) + graph.label_size(prev_id);
                        if (size_t(path_row[node_id]) == 0 || size_t(path_row[node_id]) > dist_thru) {
                            path_row[node_id] = dist_thru;
                        }
                    }
                }
            }
        }
    };
    
    uint64_t num_threads = std::max<uint64_t>(1, std::min<uint64_t>(threads, xmerge.chain_size()));
    std::vector<std::thread> workers;
    for (uint64_t t = 1; t < num_threads; ++t) {
        workers.emplace_back(compute_rows, (xmerge.chain_size() * t) / num_threads, (xmerge.chain_size() * (t + 1)) / num_threads);
    }
    compute_rows(0, xmerge.chain_size() / num_threads);
    for (auto& worker : workers) {
        worker.join();
    }
}

template<class DistVector>
inline size_t PostSwitchDistances<DistVector>::distance(uint64_t to_node_id, uint64_t from_path_id) const {
    auto d = distances[from_path_id][to_node_id];
    return d == 0 ? -1 : d;
}

template<class DistVector>
inline size_t PostSwitchDistances<DistVector>::memory_size() const {
    size_t mem = 0;
    for (const auto& row : distances) {
        mem += get_vec_memory_size(row);
    }
    return mem;
}



}

#endif /* centrolign_post_switch_distances_hpp */
