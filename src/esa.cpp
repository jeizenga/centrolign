#include <numeric>
#include <limits>
#include <unordered_set>

#include "centrolign/esa.hpp"
#include "centrolign/range_min_query.hpp"

namespace centrolign {

using namespace std;

const bool ESA::debug_esa = false;

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

std::vector<std::vector<size_t>> ESA::index_color_set_size() const {
    
    static const bool debug = false;
    
    // these will be filled out with the number of duplicate IDs in the subtree
    // note: we will only query this with LCAs, which can't be leaves, so we skip
    // the leaf annoataions
    std::vector<std::vector<size_t>> component_repeat_counts(component_ranked_ids.size());
    for (auto& counts : component_repeat_counts) {
        counts.resize(leaf_to_comp.size(), 0);
    }
    
    {
        // build an Eulerian traversal of the tree to support LCA
        std::vector<SANode> eulerian_traversal;
        std::vector<size_t> eulerian_depth;
        // note: we only query from leaves, so we skip the internal nodes here
        std::vector<size_t> position;
        position.resize(leaf_to_comp.size());
        {
            // records of (node, children, next child, depth)
            std::vector<std::tuple<SANode, std::vector<SANode>, size_t, size_t>> stack;
            
            // DFS from the root
            stack.emplace_back(root(), children(root()), 0, 0);
            while (!stack.empty()) {
                auto& top = stack.back();
                if (debug) {
                    std::cerr << "eulerian tour at " << std::get<0>(top).begin << ',' << std::get<0>(top).end << ", idx " << std::get<2>(top) << ", depth " << std::get<3>(top) << ", children";
                    for (auto child : std::get<1>(top)) {
                        std::cerr << ' ' << child.begin << ',' << child.end;
                    }
                    std::cerr << '\n';
                }
                if (std::get<0>(top).is_leaf()) {
                    position[std::get<0>(top).begin] = eulerian_traversal.size();
                }
                eulerian_depth.push_back(std::get<3>(top));
                eulerian_traversal.push_back(std::get<0>(top));
                if (std::get<2>(top) == std::get<1>(top).size()) {
                    stack.pop_back();
                }
                else {
                    auto next = std::get<1>(top)[std::get<2>(top)++];
                    stack.emplace_back(next, children(next), 0, std::get<3>(top) + 1);
                }
            }
        }
        
        if (debug) {
            std::cerr << "final Eulerian traversal of tree:\n";
            for (size_t i = 0; i < eulerian_depth.size(); ++i) {
                std::cerr << i << '\t' << eulerian_traversal[i].begin << '\t' << eulerian_traversal[i].end << '\t' << eulerian_depth[i] << '\n';
            }
        }
        
        // use eulerian traversal to make an LCA index
        RMQ<size_t> lca_rmq(eulerian_depth);
        
        // for each component, the leaf index where we last saw each ID
        std::vector<std::vector<size_t>> previous_occurrence(component_ranked_ids.size());
        
        for (size_t l = 0; l < lcp_array.size(); ++l) {
            
            auto c = leaf_to_comp[l];
            auto id = component_ranked_ids[c][nearest_comp_rank[c][l]];
            auto& comp_prev_occurrences = previous_occurrence[c];
            while (comp_prev_occurrences.size() <= id) {
                comp_prev_occurrences.push_back(-1);
            }
            
            if (comp_prev_occurrences[id] != -1) {
                // this is a duplicate occurrence
                
                // get the LCA
                auto l_prev = comp_prev_occurrences[id];
                auto pos = position[l];
                auto pos_prev = position[l_prev];
                auto lca = eulerian_traversal[lca_rmq.range_arg_min(std::min(pos, pos_prev), std::max(pos, pos_prev) + 1)];
                
                // increase the number of duplicates on this component in the LCA's subtree
                component_repeat_counts[c][st_node_annotation_idx(lca).second]++;
            }
            comp_prev_occurrences[id] = l;
        }
    }
    
    // bottom up traversal to add up subtree totals of duplicates
    {
        // aggregate the children in the parent
        auto add_child_duplicates = [&](const SANode& node, const std::vector<SANode>& children) {
            
            size_t j = st_node_annotation_idx(node).second;
            
            for (const auto& child : children) {
                if (child.is_leaf()) {
                    continue;
                }
                size_t j_child = st_node_annotation_idx(child).second;
                for (size_t c = 0; c < component_repeat_counts.size(); ++c) {
                    component_repeat_counts[c][j] += component_repeat_counts[c][j_child];
                }
            }
        };
        
        // records of (lcp, left, children)
        std::vector<std::tuple<size_t, size_t, std::vector<SANode>>> stack;
        stack.emplace_back(0, 0, std::vector<SANode>());
        for (size_t i = 1; i < lcp_array.size(); ++i) {
            
            SANode last_node(-1, -1);
            
            // figure out which internal nodes we're leaving
            size_t left = i - 1;
            while (std::get<0>(stack.back()) > lcp_array[i]) {
                
                auto& top = stack.back();
                
                // emit an internal node
                last_node = SANode(std::get<1>(top), i - 1);
                add_child_duplicates(last_node, std::get<2>(top));
                
                left = std::get<1>(top);
                stack.pop_back();
                if (std::get<0>(stack.back()) >= lcp_array[i]) {
                    // record this as a child of the parent
                    std::get<2>(stack.back()).push_back(last_node);
                    last_node = SANode(-1, -1);
                }
            }
            if (std::get<0>(stack.back()) < lcp_array[i]) {
                stack.emplace_back(lcp_array[i], left, std::vector<SANode>());
                if (last_node != SANode(-1, -1)) {
                    std::get<2>(stack.back()).push_back(last_node);
                }
            }
        }
        // clear the stack
        while (!stack.empty()) {
            auto& top = stack.back();
            
            // emit an internal node
            SANode node(std::get<1>(top), lcp_array.size() - 1);
            add_child_duplicates(node, std::get<2>(top));
            stack.pop_back();
            
            // record this as a child of the parent
            if (!stack.empty()) {
                std::get<2>(stack.back()).push_back(node);
            }
        }
    }
    
    {
        // convert from duplicate count to non-duplicate count in place
        auto convert_to_total = [&](const SANode& node) {
            
            size_t j = st_node_annotation_idx(node).second;
            
            for (size_t c = 0; c < component_repeat_counts.size(); ++c) {
                size_t total = nearest_comp_rank[c][node.end + 1] - nearest_comp_rank[c][node.begin];
                component_repeat_counts[c][j] = total - component_repeat_counts[c][j];
            }
        };
        
        // records of (lcp, left)
        std::vector<std::tuple<size_t, size_t>> stack;
        stack.emplace_back(0, 0);
        for (size_t i = 1; i < lcp_array.size(); ++i) {
            
            SANode last_node(-1, -1);
            
            // figure out which internal nodes we're leaving
            size_t left = i - 1;
            while (std::get<0>(stack.back()) > lcp_array[i]) {
                
                auto& top = stack.back();
                
                // emit an internal node
                last_node = SANode(std::get<1>(top), i - 1);
                convert_to_total(last_node);
                
                left = std::get<1>(top);
                stack.pop_back();
                if (std::get<0>(stack.back()) >= lcp_array[i]) {
                    // record this as a child of the parent
                    last_node = SANode(-1, -1);
                }
            }
            if (std::get<0>(stack.back()) < lcp_array[i]) {
                stack.emplace_back(lcp_array[i], left);
  
            }
        }
        // clear the stack
        while (!stack.empty()) {
            auto& top = stack.back();
            
            // emit an internal node
            SANode node(std::get<1>(top), lcp_array.size() - 1);
            convert_to_total(node);
            stack.pop_back();
        }
    }
    
    return component_repeat_counts;
}


}
