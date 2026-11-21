#include <numeric>
#include <limits>
#include <unordered_set>

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

std::vector<std::vector<size_t>> ESA::index_color_set_size(uint64_t threads) const {
    
    static const bool debug = false;
    
    // these will be filled out with the number of duplicate IDs in the subtree
    // note: we will only query this with LCAs, which can't be leaves, so we skip
    // the leaf annoataions
    std::vector<std::vector<size_t>> component_repeat_counts(component_ranked_ids.size());
    for (auto& counts : component_repeat_counts) {
        counts.resize(leaf_to_comp.size(), 0);
    }
    
    {
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
                
                // increase the number of duplicates on this component in the LCA's subtree
                component_repeat_counts[c][std::get<2>(*it)]++;
            }
            comp_prev_occurrences[id] = l;
        }
    }
    
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
    
    return component_repeat_counts;
}

}
