#include <numeric>
#include <limits>
#include <unordered_set>
#include <sstream>
#include <chrono>
#include <atomic>
#include <unordered_map>

#include "centrolign/esa.hpp"
#include <algorithm>

namespace centrolign {

using namespace std;

const bool ESA::debug_esa = false;
const size_t ESA::chunk_split_depth;

vector<SANode> ESA::children(const SANode& parent) const {
    vector<SANode> to_return;
    if (!parent.is_leaf()) {
        
        // get the first l-index
        size_t next_l_index = first_l_index(parent);
        // follow next l-index pointers until the last one
        to_return.emplace_back(parent.begin, next_l_index - 1);
        while (child_array_is_l_index(next_l_index)) {
            size_t curr = next_l_index;
            next_l_index = child_array[next_l_index];
            to_return.emplace_back(curr, next_l_index - 1);
        }
        to_return.emplace_back(next_l_index, parent.end);
    }
    return to_return;
}


void ESA::construct_child_array() {
    
    // one shorter, because the final position can only point "up",
    // but the "up" value gets entered in the previous index
    // TODO: could i get rid of the special cases by adding a special up value
    // to the final position?
    child_array.resize(lcp_array.size() - 1, -1);
    
    if (debug_esa) {
        cerr << "computing next l indexes\n";
    }
    
    // get the next l-indexes
    vector<size_t> stack;
    stack.push_back(0);
    for (size_t i = 1; i < lcp_array.size(); ++i) {
        // ensure invariant that LCPs are increasing over the stacked indexes, with intervening
        // LCPs being strictly longer
        while (lcp_array[stack.back()] > lcp_array[i]) {
            stack.pop_back();
        }
        if (lcp_array[i] == lcp_array[stack.back()]) {
            // we are at the next l-index
            child_array[stack.back()] = i;
            // clear this index out of the way now that we've found its next l-index
            stack.pop_back();
        }
        stack.push_back(i);
    }
    
    
    if (debug_esa) {
        cerr << "finished computing next l indexes, child array state:\n";
        for (auto i : child_array) {
            cerr << ' ';
            if (i == -1) {
                cerr << '.';
            }
            else {
                cerr << i;
            }
        }
        cerr << '\n';
        cerr << "computing up and down values\n";
    }
    // get the up and down indexes
    stack.clear();
    stack.push_back(0);
    for (size_t i = 1; i < lcp_array.size(); ++i) {
        if (debug_esa) {
            cerr << "iteration " << i << "\nstack state:";
            for (auto x : stack) {
                cerr << ' ' << x;
            }
            cerr << '\n';
        }
        size_t last_idx = -1;
        while (lcp_array[stack.back()] > lcp_array[i]) {
            last_idx = stack.back();
            stack.pop_back();
            if (child_array[stack.back()] == -1 && // not already holding a next l-index
                lcp_array[i] <= lcp_array[stack.back()] && // is latest, later blocked by equality at i
                lcp_array[last_idx] != lcp_array[stack.back()]) { // inequality from increasing stack invariant is strict
                                                                  // set childtab[stack.back()].down
                if (debug_esa) {
                    cerr << "record [" << stack.back() <<"].down = " << last_idx << '\n';
                }
                child_array[stack.back()] = last_idx;
            }
        }
        if (last_idx != -1) {
            // final item popped is earliest, must be strict inequality, all intervening have longer LCP
            // by the stack invariant
            // set childtab[i].up (but stored in position i - 1, which is guaranteed to be empty)
            if (debug_esa) {
                cerr << "record [" << i << "].up = " << last_idx << '\n';
            }
            child_array[i - 1] = last_idx;
            last_idx = -1;
        }
        stack.push_back(i);
    }
    // note: this extra step of stack clearing is necessary if you don't have an end sentinel that
    // sorts to the end like in the paper (because then the final interval might not have fallen
    // out of scope yet)
    if (debug_esa) {
        cerr << "attempting to clear stack above LCP 0\n";
    }
    while (lcp_array[stack.back()] > 0) {
        size_t last_idx = stack.back();
        stack.pop_back();
        if (child_array[stack.back()] == -1 && // not already holding a next l-index
            lcp_array[last_idx] != lcp_array[stack.back()]) { // inequality from increasing stack invariant is strict
                                                              // set childtab[stack.back()].down
            if (debug_esa) {
                cerr << "record [" << stack.back() <<"].down = " << last_idx << '\n';
            }
            child_array[stack.back()] = last_idx;
        }
    }
    
    if (debug_esa) {
        cerr << "finished computing up and down values, child array state:\n";
        for (auto i : child_array) {
            cerr << ' ';
            if (i == -1) {
                cerr << '.';
            }
            else {
                cerr << i;
            }
        }
        cerr << '\n';
    }
}

void ESA::split_chunks(uint64_t threads, std::vector<size_t>& splits, std::vector<size_t>& item_begins) const {
    
    splits.clear();
    splits.push_back(0);
    for (size_t i = 1; i < lcp_array.size(); ++i) {
        if (lcp_array[i] < chunk_split_depth) {
            splits.push_back(i);
        }
    }
    splits.push_back(lcp_array.size());
    const size_t num_chunks = splits.size() - 1;
    
    // group consecutive chunks into work items of a target size
    const size_t target_size = std::max<size_t>(1, lcp_array.size() / (16 * threads));
    item_begins.clear();
    item_begins.push_back(0);
    for (size_t k = 0; k < num_chunks; ++k) {
        if (splits[k + 1] - splits[item_begins.back()] >= target_size) {
            item_begins.push_back(k + 1);
        }
    }
    if (item_begins.back() != num_chunks) {
        item_begins.push_back(num_chunks);
    }
}

namespace {

// TODO: temporary, used to explore the parallelization of the color set size index

/*
 * A hierarchical range min query over an array, with low memory and parallel construction,
 * at the cost of O(B log_B N) query time
 */
class HierarchicalRMQ {
public:
    
    HierarchicalRMQ(const std::vector<size_t>& arr, uint64_t threads) : arr(arr) {
        // reserve so that we don't invalidate the reference to the level below
        levels.reserve(64);
        const std::vector<size_t>* below = &arr;
        while (below->size() > 2 * block_size) {
            const auto& lower = *below;
            levels.emplace_back((lower.size() + block_size - 1) / block_size);
            auto& level = levels.back();
            size_t num_chunks = std::max<size_t>(1, std::min<size_t>(level.size(), 4 * threads));
            parallel_for(num_chunks, threads, [&](uint64_t c) {
                for (size_t i = level.size() * c / num_chunks, n = level.size() * (c + 1) / num_chunks; i < n; ++i) {
                    size_t val = std::numeric_limits<size_t>::max();
                    for (size_t j = i * block_size, m = std::min(lower.size(), j + block_size); j < m; ++j) {
                        val = std::min(val, lower[j]);
                    }
                    level[i] = val;
                }
            });
            below = &level;
        }
    }
    
    // the position of a minimum in [begin, end)
    size_t range_arg_min(size_t begin, size_t end) const {
        return range_min(0, begin, end).second;
    }
    
    size_t memory_size() const {
        size_t size = sizeof(*this) + levels.capacity() * sizeof(std::vector<size_t>);
        for (const auto& level : levels) {
            size += level.capacity() * sizeof(size_t);
        }
        return size;
    }
    
private:
    
    static const size_t block_size = 32;
    
    inline const std::vector<size_t>& level(size_t k) const {
        return k == 0 ? arr : levels[k - 1];
    }
    
    // returns the min value and its index in the range [begin, end) at level k
    std::pair<size_t, size_t> range_min(size_t k, size_t begin, size_t end) const {
        const auto& a = level(k);
        std::pair<size_t, size_t> best(std::numeric_limits<size_t>::max(), -1);
        auto scan = [&](size_t b, size_t e) {
            for (size_t i = b; i < e; ++i) {
                if (a[i] < best.first) {
                    best = std::make_pair(a[i], i);
                }
            }
        };
        size_t block_begin = (begin + block_size - 1) / block_size;
        size_t block_end = end / block_size;
        if (k == levels.size() || end - begin <= 2 * block_size || block_begin >= block_end) {
            scan(begin, end);
        }
        else {
            scan(begin, block_begin * block_size);
            scan(block_end * block_size, end);
            auto above = range_min(k + 1, block_begin, block_end);
            if (above.first < best.first) {
                // descend into the block to find the index at this level
                size_t i = above.second * block_size;
                while (a[i] != above.first) {
                    ++i;
                }
                best = std::make_pair(above.first, i);
            }
        }
        return best;
    }
    
    const std::vector<size_t>& arr;
    std::vector<std::vector<size_t>> levels;
};

const size_t HierarchicalRMQ::block_size;

// an order-independent hash of the counts
size_t counts_hash(const std::vector<std::vector<size_t>>& counts, uint64_t threads) {
    size_t hsh = 0;
    for (size_t c = 0; c < counts.size(); ++c) {
        const auto& comp_counts = counts[c];
        size_t num_blocks = std::max<size_t>(1, std::min<size_t>(comp_counts.size(), 4 * threads));
        std::vector<size_t> block_hashes(num_blocks, 0);
        parallel_for(num_blocks, threads, [&](uint64_t b) {
            for (size_t i = comp_counts.size() * b / num_blocks, n = comp_counts.size() * (b + 1) / num_blocks; i < n; ++i) {
                if (comp_counts[i] != 0) {
                    // splitmix64 finalizer
                    size_t x = (i * counts.size() + c) * 0x9e3779b97f4a7c15ull + comp_counts[i];
                    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull;
                    x = (x ^ (x >> 27)) * 0x94d049bb133111ebull;
                    block_hashes[b] += x ^ (x >> 31);
                }
            }
        });
        for (auto h : block_hashes) {
            hsh += h;
        }
    }
    return hsh;
}

void zero_counts(std::vector<std::vector<size_t>>& counts, uint64_t threads) {
    for (auto& comp_counts : counts) {
        size_t num_blocks = std::max<size_t>(1, std::min<size_t>(comp_counts.size(), 4 * threads));
        parallel_for(num_blocks, threads, [&](uint64_t b) {
            std::fill(comp_counts.begin() + comp_counts.size() * b / num_blocks,
                      comp_counts.begin() + comp_counts.size() * (b + 1) / num_blocks, 0);
        });
    }
}

inline double seconds_since(const std::chrono::steady_clock::time_point& begin) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
}

}

void ESA::count_duplicates_serial(std::vector<std::vector<size_t>>& counts) const {
    
    // TODO: temporary logging to assess a parallel strategy that resolves pairs with shallow LCAs serially
    static const std::vector<size_t> log_depths{2, 4, 6, 8, 12, 16, 24, 32};
    std::vector<size_t> num_shallow_pairs(log_depths.size(), 0);
    size_t num_pairs = 0;
    
    // for each component, the leaf index where we last saw each ID
    std::vector<std::vector<size_t>> previous_occurrence(component_ranked_ids.size());
    
    // records of (lcp, left, first l-index) for the internal nodes that contain the current leaf,
    // which are ordered by both depth and left boundary
    std::vector<std::tuple<size_t, size_t, size_t>> stack;
    if (lcp_array.size() > 1) {
        stack.emplace_back(0, 0, first_l_index(root()));
    }
    for (size_t l = 0; l < lcp_array.size(); ++l) {
        
        if (l != 0) {
            // update the stack to the internal nodes that contain this leaf
            size_t left = l - 1;
            while (std::get<0>(stack.back()) > lcp_array[l]) {
                left = std::get<1>(stack.back());
                stack.pop_back();
            }
            if (std::get<0>(stack.back()) < lcp_array[l]) {
                // this is the first l-index of the new internal node
                stack.emplace_back(lcp_array[l], left, l);
            }
        }
        
        auto c = leaf_to_comp[l];
        auto id = component_ranked_ids[c][nearest_comp_rank[c][l]];
        auto& comp_prev_occurrences = previous_occurrence[c];
        while (comp_prev_occurrences.size() <= id) {
            comp_prev_occurrences.push_back(-1);
        }
        
        if (comp_prev_occurrences[id] != -1) {
            // this is a duplicate occurrence
            
            // the LCA is the deepest internal node containing this leaf that also contains the
            // previous occurrence
            auto l_prev = comp_prev_occurrences[id];
            auto it = std::upper_bound(stack.begin(), stack.end(), l_prev,
                                       [](size_t val, const std::tuple<size_t, size_t, size_t>& rec) {
                return val < std::get<1>(rec);
            });
            --it;
            
            ++num_pairs;
            for (size_t k = 0; k < log_depths.size(); ++k) {
                num_shallow_pairs[k] += (std::get<0>(*it) < log_depths[k]);
            }
            
            // increase the number of duplicates on this component in the LCA's subtree
            counts[c][std::get<2>(*it)]++;
        }
        comp_prev_occurrences[id] = l;
    }
    
    if (logging::level >= logging::Debug) {
        std::stringstream strm;
        strm << "Color set size found " << num_pairs << " repeated pairs. Pairs with LCA shallower than depth:";
        for (size_t k = 0; k < log_depths.size(); ++k) {
            strm << ' ' << log_depths[k] << " (" << num_shallow_pairs[k] << ")";
        }
        logging::log(logging::Debug, strm.str());
    }
}

void ESA::count_duplicates_chunked(std::vector<std::vector<size_t>>& counts, uint64_t threads) const {
    
    // pairs of leaves whose LCA is deep are resolved within the chunk, and pairs whose LCA is shallow
    // or whose previous occurrence is in an earlier chunk are queued and resolved serially
    
    auto begin_time = std::chrono::steady_clock::now();
    
    std::vector<size_t> splits, item_begins;
    split_chunks(threads, splits, item_begins);
    const size_t num_items = item_begins.size() - 1;
    
    // for each work item, the queued leaves as (leaf, previous occurrence or -1 if it is in an earlier
    // chunk), and the last occurrence of each ID in each chunk, both in chunk order
    std::vector<std::vector<std::pair<size_t, size_t>>> item_queues(num_items);
    std::vector<std::vector<size_t>> item_lasts(num_items);
    std::atomic<size_t> num_resolved(0);
    std::atomic<size_t> num_shallow_queued(0);
    
    parallel_for(num_items, threads, [&](uint64_t w) {
        auto& queue = item_queues[w];
        auto& lasts = item_lasts[w];
        size_t resolved = 0, shallow_queued = 0;
        
        // the last occurrence of each (component, ID) in the current chunk
        std::unordered_map<std::pair<size_t, uint64_t>, size_t> last_occurrence;
        // records of (lcp, left, first l-index) above a floor at the chunk's start
        std::vector<std::tuple<size_t, size_t, size_t>> stack;
        for (size_t k = item_begins[w]; k < item_begins[w + 1]; ++k) {
            size_t b = splits[k], e = splits[k + 1];
            last_occurrence.clear();
            stack.clear();
            stack.emplace_back(chunk_split_depth - 1, b, -1);
            for (size_t l = b; l < e; ++l) {
                if (l != b) {
                    // update the stack to the internal nodes that contain this leaf
                    size_t left = l - 1;
                    while (std::get<0>(stack.back()) > lcp_array[l]) {
                        left = std::get<1>(stack.back());
                        stack.pop_back();
                    }
                    if (std::get<0>(stack.back()) < lcp_array[l]) {
                        stack.emplace_back(lcp_array[l], left, l);
                    }
                }
                
                auto c = leaf_to_comp[l];
                auto id = component_ranked_ids[c][nearest_comp_rank[c][l]];
                auto key = std::make_pair(size_t(c), id);
                auto f = last_occurrence.find(key);
                if (f == last_occurrence.end()) {
                    // the previous occurrence, if any, is in an earlier chunk
                    queue.emplace_back(l, -1);
                    last_occurrence[key] = l;
                }
                else {
                    auto l_prev = f->second;
                    auto it = std::upper_bound(stack.begin(), stack.end(), l_prev,
                                               [](size_t val, const std::tuple<size_t, size_t, size_t>& rec) {
                        return val < std::get<1>(rec);
                    });
                    --it;
                    if (it != stack.begin()) {
                        // the LCA is a deep node in this chunk, so its first l-index is too
                        counts[c][std::get<2>(*it)]++;
                        ++resolved;
                    }
                    else {
                        // the LCA is shallow
                        queue.emplace_back(l, l_prev);
                        ++shallow_queued;
                    }
                    f->second = l;
                }
            }
            for (const auto& rec : last_occurrence) {
                lasts.push_back(rec.second);
            }
        }
        num_resolved += resolved;
        num_shallow_queued += shallow_queued;
    });
    
    double parallel_time = seconds_since(begin_time);
    auto serial_begin_time = std::chrono::steady_clock::now();
    
    // resolve the queue with the shallow nodes
    size_t num_cross_pairs = 0;
    std::vector<std::vector<size_t>> previous_occurrence(component_ranked_ids.size());
    // records of (lcp, left, first l-index) for the shallow nodes
    std::vector<std::tuple<size_t, size_t, size_t>> stack;
    if (lcp_array.size() > 1) {
        stack.emplace_back(0, 0, first_l_index(root()));
    }
    for (size_t w = 0; w < num_items; ++w) {
        const auto& queue = item_queues[w];
        const auto& lasts = item_lasts[w];
        size_t q = 0, r = 0;
        for (size_t k = item_begins[w]; k < item_begins[w + 1]; ++k) {
            if (k != 0) {
                // update the stack to the shallow nodes that contain this chunk, where the
                // previous chunk's top node (if any) is the leftmost deep node we leave
                size_t i = splits[k];
                size_t left = splits[k - 1];
                while (std::get<0>(stack.back()) > lcp_array[i]) {
                    left = std::get<1>(stack.back());
                    stack.pop_back();
                }
                if (std::get<0>(stack.back()) < lcp_array[i]) {
                    stack.emplace_back(lcp_array[i], left, i);
                }
            }
            
            for (; q < queue.size() && queue[q].first < splits[k + 1]; ++q) {
                size_t l = queue[q].first;
                auto c = leaf_to_comp[l];
                size_t l_prev = queue[q].second;
                if (l_prev == -1) {
                    auto id = component_ranked_ids[c][nearest_comp_rank[c][l]];
                    auto& comp_prev_occurrences = previous_occurrence[c];
                    if (id < comp_prev_occurrences.size()) {
                        l_prev = comp_prev_occurrences[id];
                    }
                    if (l_prev == -1) {
                        // this is the first occurrence
                        continue;
                    }
                    ++num_cross_pairs;
                }
                auto it = std::upper_bound(stack.begin(), stack.end(), l_prev,
                                           [](size_t val, const std::tuple<size_t, size_t, size_t>& rec) {
                    return val < std::get<1>(rec);
                });
                --it;
                counts[c][std::get<2>(*it)]++;
            }
            for (; r < lasts.size() && lasts[r] < splits[k + 1]; ++r) {
                size_t l = lasts[r];
                auto c = leaf_to_comp[l];
                auto id = component_ranked_ids[c][nearest_comp_rank[c][l]];
                auto& comp_prev_occurrences = previous_occurrence[c];
                while (comp_prev_occurrences.size() <= id) {
                    comp_prev_occurrences.push_back(-1);
                }
                comp_prev_occurrences[id] = l;
            }
        }
    }
    
    double serial_time = seconds_since(serial_begin_time);
    
    if (logging::level >= logging::Debug) {
        size_t queue_size = 0, lasts_size = 0, list_mem = 0;
        for (size_t w = 0; w < num_items; ++w) {
            queue_size += item_queues[w].size();
            lasts_size += item_lasts[w].size();
            list_mem += (item_queues[w].capacity() * sizeof(std::pair<size_t, size_t>)
                         + item_lasts[w].capacity() * sizeof(size_t));
        }
        std::stringstream strm;
        strm << "Chunked CSS duplicate counting with " << (splits.size() - 1) << " chunks in " << num_items << " work items: resolved " << num_resolved.load() << " pairs in parallel in " << parallel_time << " s, then " << num_shallow_queued.load() << " in-chunk pairs with shallow LCAs and " << num_cross_pairs << " cross-chunk pairs serially in " << serial_time << " s. Queue has " << queue_size << " entries and last occurrences have " << lasts_size << " entries, occupying " << format_memory_usage(list_mem) << ".";
        logging::log(logging::Debug, strm.str());
    }
}

void ESA::count_duplicates_rmq(std::vector<std::vector<size_t>>& counts, uint64_t threads) const {
    
    // each pair of leaves is resolved independently with an RMQ over the LCP array, which identifies
    // an l-index of the LCA (but not necessarily the first one)
    
    auto begin_time = std::chrono::steady_clock::now();
    
    HierarchicalRMQ rmq(lcp_array, threads);
    
    double build_time = seconds_since(begin_time);
    auto scan_begin_time = std::chrono::steady_clock::now();
    
    auto add_pair = [&](size_t c, size_t l_prev, size_t l) {
        size_t i = rmq.range_arg_min(l_prev + 1, l + 1);
        __atomic_fetch_add(&counts[c][i], 1, __ATOMIC_RELAXED);
    };
    
    // find pairs within blocks of leaves, and record the first and last occurrences in each block
    const size_t num_blocks = std::max<size_t>(1, std::min<size_t>(lcp_array.size(), 4 * threads));
    std::vector<std::vector<size_t>> block_firsts(num_blocks), block_lasts(num_blocks);
    std::atomic<size_t> num_in_block_pairs(0);
    parallel_for(num_blocks, threads, [&](uint64_t b) {
        auto& firsts = block_firsts[b];
        auto& lasts = block_lasts[b];
        size_t num_pairs = 0;
        std::unordered_map<std::pair<size_t, uint64_t>, size_t> last_occurrence;
        for (size_t l = lcp_array.size() * b / num_blocks, n = lcp_array.size() * (b + 1) / num_blocks; l < n; ++l) {
            auto c = leaf_to_comp[l];
            auto id = component_ranked_ids[c][nearest_comp_rank[c][l]];
            auto key = std::make_pair(size_t(c), id);
            auto f = last_occurrence.find(key);
            if (f == last_occurrence.end()) {
                firsts.push_back(l);
                last_occurrence[key] = l;
            }
            else {
                add_pair(c, f->second, l);
                f->second = l;
                ++num_pairs;
            }
        }
        for (const auto& rec : last_occurrence) {
            lasts.push_back(rec.second);
        }
        num_in_block_pairs += num_pairs;
    });
    
    double scan_time = seconds_since(scan_begin_time);
    auto link_begin_time = std::chrono::steady_clock::now();
    
    // link the first occurrences in each block to the last occurrences in earlier blocks
    std::vector<std::pair<size_t, size_t>> cross_pairs;
    {
        std::vector<std::vector<size_t>> previous_occurrence(component_ranked_ids.size());
        for (size_t b = 0; b < num_blocks; ++b) {
            for (auto l : block_firsts[b]) {
                auto c = leaf_to_comp[l];
                auto id = component_ranked_ids[c][nearest_comp_rank[c][l]];
                const auto& comp_prev_occurrences = previous_occurrence[c];
                if (id < comp_prev_occurrences.size() && comp_prev_occurrences[id] != -1) {
                    cross_pairs.emplace_back(comp_prev_occurrences[id], l);
                }
            }
            for (auto l : block_lasts[b]) {
                auto c = leaf_to_comp[l];
                auto id = component_ranked_ids[c][nearest_comp_rank[c][l]];
                auto& comp_prev_occurrences = previous_occurrence[c];
                while (comp_prev_occurrences.size() <= id) {
                    comp_prev_occurrences.push_back(-1);
                }
                comp_prev_occurrences[id] = l;
            }
        }
    }
    size_t list_mem = cross_pairs.capacity() * sizeof(std::pair<size_t, size_t>);
    for (size_t b = 0; b < num_blocks; ++b) {
        list_mem += (block_firsts[b].capacity() + block_lasts[b].capacity()) * sizeof(size_t);
    }
    
    double link_time = seconds_since(link_begin_time);
    auto cross_begin_time = std::chrono::steady_clock::now();
    
    {
        size_t num_chunks = std::max<size_t>(1, std::min<size_t>(cross_pairs.size(), 4 * threads));
        parallel_for(num_chunks, threads, [&](uint64_t k) {
            for (size_t i = cross_pairs.size() * k / num_chunks, n = cross_pairs.size() * (k + 1) / num_chunks; i < n; ++i) {
                add_pair(leaf_to_comp[cross_pairs[i].second], cross_pairs[i].first, cross_pairs[i].second);
            }
        });
    }
    
    double cross_time = seconds_since(cross_begin_time);
    auto fold_begin_time = std::chrono::steady_clock::now();
    
    // move the counts to the first l-index of each node
    auto fold = [&](size_t first, size_t i) {
        for (auto& comp_counts : counts) {
            comp_counts[first] += comp_counts[i];
            comp_counts[i] = 0;
        }
    };
    std::vector<size_t> splits, item_begins;
    split_chunks(threads, splits, item_begins);
    parallel_for(item_begins.size() - 1, threads, [&](uint64_t w) {
        // records of (lcp, first l-index) above a floor at the chunk's start
        std::vector<std::pair<size_t, size_t>> stack;
        for (size_t k = item_begins[w]; k < item_begins[w + 1]; ++k) {
            stack.clear();
            stack.emplace_back(chunk_split_depth - 1, -1);
            for (size_t i = splits[k] + 1; i < splits[k + 1]; ++i) {
                while (stack.back().first > lcp_array[i]) {
                    stack.pop_back();
                }
                if (stack.back().first < lcp_array[i]) {
                    stack.emplace_back(lcp_array[i], i);
                }
                else {
                    // a later l-index of a deep node
                    fold(stack.back().second, i);
                }
            }
        }
    });
    {
        // the shallow nodes
        std::vector<std::pair<size_t, size_t>> stack;
        if (lcp_array.size() > 1) {
            stack.emplace_back(0, first_l_index(root()));
        }
        for (size_t k = 1; k + 1 < splits.size(); ++k) {
            size_t i = splits[k];
            while (stack.back().first > lcp_array[i]) {
                stack.pop_back();
            }
            if (stack.back().first < lcp_array[i]) {
                stack.emplace_back(lcp_array[i], i);
            }
            else if (stack.back().second != i) {
                fold(stack.back().second, i);
            }
        }
    }
    
    double fold_time = seconds_since(fold_begin_time);
    
    if (logging::level >= logging::Debug) {
        std::stringstream strm;
        strm << "RMQ CSS duplicate counting built RMQ occupying " << format_memory_usage(rmq.memory_size()) << " in " << build_time << " s, resolved " << num_in_block_pairs.load() << " pairs within " << num_blocks << " blocks in parallel in " << scan_time << " s, linked " << cross_pairs.size() << " cross-block pairs serially in " << link_time << " s, resolved them in parallel in " << cross_time << " s, and moved counts to first l-indexes in " << fold_time << " s. Block lists occupy " << format_memory_usage(list_mem) << ".";
        logging::log(logging::Debug, strm.str());
    }
}

std::vector<std::vector<size_t>> ESA::index_color_set_size(uint64_t threads) const {
    
    static const bool debug = false;
    
    // these will be filled out with the number of duplicate IDs in the subtree
    // note: we will only query this with LCAs, which can't be leaves, so we skip
    // the leaf annoataions
    std::vector<std::vector<size_t>> component_repeat_counts(component_ranked_ids.size());
    for (auto& counts : component_repeat_counts) {
        counts.resize(leaf_to_comp.size(), 0);
    }
    
    // TODO: temporary comparison of strategies for counting duplicates
    {
        auto begin_time = std::chrono::steady_clock::now();
        count_duplicates_serial(component_repeat_counts);
        double serial_time = seconds_since(begin_time);
        
        if (threads > 1) {
            size_t serial_hash = counts_hash(component_repeat_counts, threads);
            
            zero_counts(component_repeat_counts, threads);
            begin_time = std::chrono::steady_clock::now();
            count_duplicates_chunked(component_repeat_counts, threads);
            double chunked_time = seconds_since(begin_time);
            bool chunked_match = (counts_hash(component_repeat_counts, threads) == serial_hash);
            
            zero_counts(component_repeat_counts, threads);
            begin_time = std::chrono::steady_clock::now();
            count_duplicates_rmq(component_repeat_counts, threads);
            double rmq_time = seconds_since(begin_time);
            bool rmq_match = (counts_hash(component_repeat_counts, threads) == serial_hash);
            
            if (!chunked_match || !rmq_match) {
                logging::log(logging::Basic, std::string("warning: parallel CSS duplicate counting disagrees with serial counting (chunked ") + (chunked_match ? "agrees" : "disagrees") + ", RMQ " + (rmq_match ? "agrees" : "disagrees") + "), falling back to serial counts");
                zero_counts(component_repeat_counts, threads);
                count_duplicates_serial(component_repeat_counts);
            }
            
            if (logging::level >= logging::Debug) {
                std::stringstream strm;
                strm << "CSS duplicate counting took " << serial_time << " s serially, " << chunked_time << " s with chunking, and " << rmq_time << " s with RMQ using " << threads << " threads.";
                logging::log(logging::Debug, strm.str());
            }
        }
        else if (logging::level >= logging::Debug) {
            logging::log(logging::Debug, "CSS duplicate counting took " + std::to_string(serial_time) + " s serially.");
        }
    }
    
    auto begin_time = std::chrono::steady_clock::now();
    
    std::vector<size_t> splits, item_begins;
    split_chunks(threads, splits, item_begins);
    
    // bottom up traversal to add up subtree totals of duplicates
    chunked_bottom_up_traversal<uint8_t>(splits, item_begins, threads,
                                         [&](const SANode& node, const std::vector<SANode>& children,
                                             const std::vector<uint8_t>& child_values, size_t item) {
        
        // aggregate the children in the parent
        size_t j = st_node_annotation_idx(node).second;
        for (const auto& child : children) {
            size_t j_child = st_node_annotation_idx(child).second;
            for (size_t c = 0; c < component_repeat_counts.size(); ++c) {
                component_repeat_counts[c][j] += component_repeat_counts[c][j_child];
            }
        }
        return uint8_t(0);
    });
    
    double sum_time = seconds_since(begin_time);
    begin_time = std::chrono::steady_clock::now();
    
    // convert from duplicate count to non-duplicate count in place
    chunked_bottom_up_traversal<uint8_t>(splits, item_begins, threads,
                                         [&](const SANode& node, const std::vector<SANode>& children,
                                             const std::vector<uint8_t>& child_values, size_t item) {
        
        size_t j = st_node_annotation_idx(node).second;
        for (size_t c = 0; c < component_repeat_counts.size(); ++c) {
            size_t total = nearest_comp_rank[c][node.end + 1] - nearest_comp_rank[c][node.begin];
            component_repeat_counts[c][j] = total - component_repeat_counts[c][j];
        }
        return uint8_t(0);
    });
    
    if (logging::level >= logging::Debug) {
        std::stringstream strm;
        strm << "CSS subtree sums took " << sum_time << " s and conversion to totals took " << seconds_since(begin_time) << " s.";
        logging::log(logging::Debug, strm.str());
    }
    
    return component_repeat_counts;
}

}
