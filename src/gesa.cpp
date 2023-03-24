#include <numeric>

#include "centrolign/gesa.hpp"

namespace centrolign {

using namespace std;

bool GESA::debug_gesa = false;

size_t GESA::component_size() const {
    return component_sizes.size();
}

GESANode GESA::root() const {
    return GESANode(0, lcp_array.size() - 1);
}

GESANode GESA::link(const GESANode& node) const {
    size_t i, j;
    tie(i, j) = st_node_annotation_idx(node);
    return suffix_links[i][j];
}

vector<GESANode> GESA::children(const GESANode& parent) const {
    vector<GESANode> to_return;
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

void GESA::construct_child_array() {
    
    // one shorter, because the final position can only point "up",
    // but the "up" value gets entered in the previous index
    // TODO: could i get rid of the special cases by adding a special up value
    // to the final position?
    child_array.resize(lcp_array.size() - 1, -1);
    
    if (debug_gesa) {
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
    
    
    if (debug_gesa) {
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
        if (debug_gesa) {
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
//            if (debug_gesa) {
//                cerr << "destack " << last_idx << ", stack state:";
//                for (auto x : stack) {
//                    cerr << ' ' << x;
//                }
//                cerr << '\n';
//            }
            if (child_array[stack.back()] == -1 && // not already holding a next l-index
                lcp_array[i] <= lcp_array[stack.back()] && // is latest, later blocked by equality at i
                lcp_array[last_idx] != lcp_array[stack.back()]) { // inequality from increasing stack invariant is strict
                // set childtab[stack.back()].down
                if (debug_gesa) {
                    cerr << "record [" << stack.back() <<"].down = " << last_idx << '\n';
                }
                child_array[stack.back()] = last_idx;
            }
        }
        if (last_idx != -1) {
            // final item popped is earliest, must be strict inequality, all intervening have longer LCP
            // by the stack invariant
            // set childtab[i].up (but stored in position i - 1, which is guaranteed to be empty)
            if (debug_gesa) {
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
    if (debug_gesa) {
        cerr << "attempting to clear stack above LCP 0\n";
    }
    while (lcp_array[stack.back()] > 0) {
        size_t last_idx = stack.back();
        stack.pop_back();
//        if (debug_gesa) {
//            cerr << "destack " << last_idx << ", stack state:";
//            for (auto x : stack) {
//                cerr << ' ' << x;
//            }
//            cerr << '\n';
//        }
        if (child_array[stack.back()] == -1 && // not already holding a next l-index
            lcp_array[last_idx] != lcp_array[stack.back()]) { // inequality from increasing stack invariant is strict
                                                              // set childtab[stack.back()].down
            if (debug_gesa) {
                cerr << "record [" << stack.back() <<"].down = " << last_idx << '\n';
            }
            child_array[stack.back()] = last_idx;
        }
    }
    
    if (debug_gesa) {
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

void GESA::construct_suffix_links() {
    
    if (debug_gesa) {
        cerr << "computing suffix links\n";
    }
    
    // TODO: decide if I want the leaf suffix links
    
    // records of (lcp, left, right)
    vector<tuple<size_t, size_t, size_t>> stack;
    
    // do a top down traversal to construct lists of l-intervals
    vector<vector<GESANode>> internal_l_interval_lists, leaf_l_interval_lists;
    stack.emplace_back(0, 0, -1);
    size_t first_l = 1;
    if (lcp_array.size() > 1) {
        first_l = lcp_array[1] + 1;
    }
    while (leaf_l_interval_lists.size() <= first_l) {
        internal_l_interval_lists.emplace_back();
        leaf_l_interval_lists.emplace_back();
    }
    if (debug_gesa) {
        cerr << "emit leaf [" << 0 << ", " << 0 << "], l = " << first_l << '\n';
    }
    leaf_l_interval_lists[first_l].emplace_back(0, 0);
    for (size_t i = 1; i < lcp_array.size(); ++i) {
        if (debug_gesa) {
            cerr << "iteration " << i << ", stack state:\n";
            for (auto x : stack) {
                cerr << " (" << get<0>(x) << ',' << get<1>(x) << ',' << (int64_t) get<2>(x) << ')';
            }
            cerr << '\n';
        }
        // get the unique prefix of the leaf
        size_t leaf_l = lcp_array[i];
        if (i + 1 < lcp_array.size()) {
            leaf_l = std::max(leaf_l, lcp_array[i + 1]);
        }
        ++leaf_l;
        while (leaf_l_interval_lists.size() <= leaf_l) {
            internal_l_interval_lists.emplace_back();
            leaf_l_interval_lists.emplace_back();
        }
        if (debug_gesa) {
            cerr << "emit leaf [" << i << ", " << i << "], l = " << leaf_l << '\n';
        }
        leaf_l_interval_lists[leaf_l].emplace_back(i, i);
        
        // figure out which internal nodes we're leaving
        size_t left = i - 1;
        while (get<0>(stack.back()) > lcp_array[i]) {
            auto& top = stack.back();
            get<2>(top) = i - 1;
            if (debug_gesa) {
                cerr << "emit internal [" << get<1>(top) << ", " << get<2>(top) << "], l = " << get<0>(top) << '\n';
            }
            while (internal_l_interval_lists.size() <= get<0>(top)) {
                internal_l_interval_lists.emplace_back();
                leaf_l_interval_lists.emplace_back();
            }
            internal_l_interval_lists[get<0>(top)].emplace_back(get<1>(top), get<2>(top));
            left = get<1>(top);
            stack.pop_back();
        }
        if (lcp_array[i] > get<0>(stack.back())) {
            stack.emplace_back(lcp_array[i], left, -1);
        }
    }
    if (debug_gesa) {
        cerr << "clearing the stack\n";
    }
    // clear the stack
    while (!stack.empty()) {
        auto& top = stack.back();
        get<2>(top) = lcp_array.size() - 1;
        if (debug_gesa) {
            cerr << "emit internal [" << get<1>(top) << ", " << get<2>(top) << "], l = " << get<0>(top) << '\n';
        }
        while (internal_l_interval_lists.size() <= get<0>(top)) {
            internal_l_interval_lists.emplace_back();
            leaf_l_interval_lists.emplace_back();
        }
        internal_l_interval_lists[get<0>(top)].emplace_back(get<1>(top), get<2>(top));
        stack.pop_back();
    }
    
    // merge the leaf and internal node lists
    vector<vector<GESANode>> l_interval_lists(internal_l_interval_lists.size());
    for (size_t l = 0; l < l_interval_lists.size(); ++l) {
        auto& merge_list = l_interval_lists[l];
        auto& list1 = internal_l_interval_lists[l];
        auto& list2 = leaf_l_interval_lists[l];
        size_t i = 0, j = 0;
        while (i < list1.size() && j < list2.size()) {
            if (list1[i] < list2[j]) {
                merge_list.emplace_back(list1[i++]);
            }
            else {
                merge_list.push_back(list2[j++]);
            }
        }
        for (; i < list1.size(); ++i) {
            merge_list.push_back(list1[i]);
        }
        for (; j < list2.size(); ++j) {
            merge_list.push_back(list2[j]);
        }
    }
    
    if (debug_gesa) {
        cerr << "l interval lists:\n";
        for (size_t l = 0; l < l_interval_lists.size(); ++l) {
            cerr << l << ':';
            for (auto node : l_interval_lists[l]) {
                cerr << " [" << node.begin << ',' << node.end << "]";
            }
            cerr << '\n';
        }
    }
    
    suffix_links[0].resize(lcp_array.size());
    suffix_links[1].resize(lcp_array.size());
    for (size_t l = 1; l < l_interval_lists.size(); ++l) {
        if (debug_gesa) {
            cerr << "forming links for nodes at depth " << l << "\n";
        }
        auto& link_interval_list = l_interval_lists[l - 1];
        for (const GESANode& node : l_interval_lists[l]) {
            const auto& node_edges = edges[node.begin];
            size_t i, j;
            tie(i, j) = st_node_annotation_idx(node);
            if (node_edges.empty()) {
                // this is a sink node, so removing a character gives the empty string
                suffix_links[i][j] = root();
                if (debug_gesa) {
                    cerr << "link [" << node.begin << ',' << node.end << "] -> [" << root().begin << ',' << root().end << "], storing at " << i << ", " << j << "\n";
                }
            }
            else {
                // walk one edge forward to find a node that should be in the interval
                size_t next_rank = node_edges.front();
                
                // binary search to find the containing interval
                size_t lo = 0, hi = link_interval_list.size() - 1;
                while (lo != hi) {
                    size_t mid = (lo + hi) / 2;
                    if (next_rank < link_interval_list[mid].begin) {
                        hi = mid - 1;
                    }
                    else if (next_rank > link_interval_list[mid].end) {
                        lo = mid + 1;
                    }
                    else {
                        hi = lo = mid;
                    }
                }
                if (debug_gesa) {
                    cerr << "link [" << node.begin << ',' << node.end << "] -> [" << link_interval_list[lo].begin << ',' << link_interval_list[lo].end << "], storing at " << i << ", " << j << "\n";
                }
                
                // record the suffix link
                suffix_links[i][j] = link_interval_list[lo];
            }
        }
    }
}

void GESA::compute_subtree_counts() {
    
    if (debug_gesa) {
        cerr << "computing subtree counts\n";
    }
    
    // map the original node IDs to their components
    size_t total_size = accumulate(component_sizes.begin(), component_sizes.end(), 0);
    vector<size_t> node_to_comp(total_size, 0);
    uint64_t node_id = 0;
    for (size_t i = 0; i < component_sizes.size(); ++i) {
        auto& component_table = component_subtree_counts[i];
        for (size_t j = 0, n = component_sizes[i]; j < n; ++j) {
            node_to_comp[node_id] = i;
            ++node_id;
        }
    }
    
    // initialize the table
    component_subtree_counts.resize(component_size());
    for (size_t i = 0; i < component_sizes.size(); ++i) {
        component_subtree_counts[i][0].resize(lcp_array.size(), 0);
        component_subtree_counts[i][1].resize(lcp_array.size(), 0);
    }
    
    // set the leaf values to 1 for their component
    for (size_t r = 0; r < ranked_node_ids.size(); ++r) {
        size_t comp = node_to_comp[ranked_node_ids[r]];
        size_t i, j;
        tie(i, j) = st_node_annotation_idx(GESANode(r, r));
        component_subtree_counts[comp][i][j] = 1;
    }
    
    // TODO: this would probably be more efficient if i figured out how to also
    // emit leaves in this traversal rather than relying on children
    
    // compute the sum of occurences from each component from children
    auto do_dp = [&](const tuple<size_t, size_t, size_t>& record) {
        GESANode node(get<1>(record), get<2>(record));
        size_t i, j;
        tie(i, j) = st_node_annotation_idx(node);
        for (auto child : children(node)) {
            size_t ci, cj;
            tie(ci, cj) = st_node_annotation_idx(child);
            for (size_t c = 0; c < component_subtree_counts.size(); ++c) {
                component_subtree_counts[c][i][j] += component_subtree_counts[c][ci][cj];
            }
        }
    };
    
    // records of (lcp, left, right)
    vector<tuple<size_t, size_t, size_t>> stack;
    stack.emplace_back(0, 0, -1);
    for (size_t i = 1; i < lcp_array.size(); ++i) {
        
        // figure out which internal nodes we're leaving
        size_t left = i - 1;
        while (get<0>(stack.back()) > lcp_array[i]) {
            auto& top = stack.back();
            // emit an internal node
            get<2>(top) = i - 1;
            do_dp(top);
            left = get<1>(top);
            stack.pop_back();
        }
        if (lcp_array[i] > get<0>(stack.back())) {
            stack.emplace_back(lcp_array[i], left, -1);
        }
    }
    // clear the stack
    while (!stack.empty()) {
        auto& top = stack.back();
        // emit an internal node
        get<2>(top) = lcp_array.size() - 1;
        do_dp(top);
        stack.pop_back();
    }
}

void GESA::print(ostream& out) const {
    out << "i" << '\t' << "Nd" << '\t' << "LCP" << '\t' << "Ch" << '\t' << "LL1" << '\t' << "LL2" << '\t' << "IL1" << '\t' << "IL2";
    for (size_t c = 0; c < component_subtree_counts.size(); ++c) {
        out << '\t' << "SC" << c;
    }
    out << '\t' << "Es" << '\n';
    for (size_t i = 0; i < ranked_node_ids.size(); ++i) {
        out << i << '\t' << ranked_node_ids[i] << '\t' << lcp_array[i] << '\t' << (i < child_array.size() ? (int) child_array[i] : -1);
        for (auto j : {0, 1}) {
            if (suffix_links[j][i] != GESANode()) {
                out << '\t' << suffix_links[j][i].begin << '\t' << suffix_links[j][i].end;
            }
            else {
                out << '\t' << -1 << '\t' << -1;
            }
        }
        for (size_t j = 0; j < component_subtree_counts.size(); ++j) {
            out << '\t' << component_subtree_counts[j][1][i];
        }
        auto& node_edges = edges[i];
        out << '\t';
        for (size_t j = 0; j < node_edges.size(); ++j) {
            if (j) {
                out << ',';
            }
            out << node_edges[j];
        }
        out << '\n';
    }
}

}
