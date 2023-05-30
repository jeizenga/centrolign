#include <numeric>
#include <limits>
#include <unordered_set>

#include "centrolign/gesa.hpp"
#include "centrolign/range_unique_query.hpp"

namespace centrolign {

using namespace std;

bool GESA::debug_gesa = false;

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

vector<tuple<GESANode, size_t, vector<uint64_t>>> GESA::minimal_rare_matches(size_t max_count) const {
    
    if (debug_gesa) {
        cerr << "finding minimal rare matches with max count " << max_count << '\n';
    }
    
    logging::log(logging::Debug, "Constructing Range-Unique-Query structures");
    
    // construct range unique queries to compute subtree counts
    vector<RUQ> ruqs;
    ruqs.reserve(component_ranked_ids.size());
    for (const auto& ranked_ids : component_ranked_ids) {
        ruqs.emplace_back(ranked_ids);
    }
    
    vector<tuple<GESANode, size_t, vector<uint64_t>>> matches;
    auto add_matches = [&](const GESANode& parent, const vector<GESANode>& children) {
        
        // in order to have the same count as the children, we need one more character
        // after the parent's depth
        size_t unique_length = depth(parent) + 1;
        if (debug_gesa) {
            cerr << "checking non-leaf children of " << parent.begin << ',' << parent.end << ", which has depth " << depth(parent) << '\n';
        }
        
        if (unique_length == 1) {
            if (debug_gesa) {
                cerr << "parent is the root, handling as a special case\n";
            }
            // we only need to check for the max count on children of root
            for (const auto& child : children) {
                if (debug_gesa) {
                    cerr << "considering match node " << child.begin << ',' << child.end << '\n';
                }
                std::vector<uint64_t> counts(component_size());
                bool above_max = false;
                size_t num_nonzero = 0;
                for (size_t c = 0; c < component_size(); ++c) {
                    uint64_t count = ruqs[c].range_unique(nearest_comp_rank[c][child.begin],
                                                          nearest_comp_rank[c][child.end + 1]);
                    counts[c] = count;
                    if (count > max_count) {
                        // breaks max count on this component
                        if (debug_gesa) {
                            cerr << "count on component " << c << " is " << count << ", which is above max\n";
                        }
                        above_max = true;
                        break;
                    }
                    if (count) {
                        ++num_nonzero;
                    }
                }
                if (!above_max && num_nonzero == component_size()) {
                    if (debug_gesa) {
                        cerr << "is a minimal match with length " << unique_length << " and counts:";
                        for (auto cnt : counts) {
                            cerr << ' ' << cnt;
                        }
                        cerr << '\n';
                    }
                    matches.emplace_back(child, unique_length, std::move(counts));
                }
                else if (debug_gesa) {
                    cerr << "not a minimal rare match\n";
                }
            }
            return;
        }
        
        // align the children of the parent's suffix link to the children we're testing
        vector<GESANode> link_children = this->children(link(parent));
        for (size_t i = 0, j = 0; i < children.size(); ++j) {
            if (label(children[i]) == label(link_children[j])) {
                link_children[i] = link_children[j];
                ++i;
            }
        }
        link_children.resize(children.size());
        
        for (size_t k = 0; k < children.size(); ++k) {
            if (debug_gesa) {
                cerr << "considering match node " << children[k].begin << ',' << children[k].end << '\n';
            }
            
            auto& child = children[k];
            auto& link_child = link_children[k];
            
            std::vector<uint64_t> counts(component_size());
            size_t num_nonzero = 0;
            bool link_more_frequent = false;
            bool above_max = false;
            for (size_t c = 0; c < component_size(); ++c) {
                uint64_t count = ruqs[c].range_unique(nearest_comp_rank[c][child.begin],
                                                      nearest_comp_rank[c][child.end + 1]);
                counts[c] = count;
                if (count > max_count) {
                    // breaks max count on this component
                    if (debug_gesa) {
                        cerr << "count on component " << c << " is " << count << ", which is above max\n";
                    }
                    above_max = true;
                    break;
                }
                if (count) {
                    ++num_nonzero;
                }
                uint64_t link_count = ruqs[c].range_unique(nearest_comp_rank[c][link_child.begin],
                                                           nearest_comp_rank[c][link_child.end + 1]);
                link_more_frequent = (link_more_frequent || count < link_count);
            }
            if (num_nonzero > 1 && link_more_frequent && !above_max) {
                // occurs on more than one component and removing the first character
                // involves introducing more matches
                if (debug_gesa) {
                    cerr << "is a minimal match with length " << unique_length << " and counts:";
                    for (auto cnt : counts) {
                        cerr << ' ' << cnt;
                    }
                    cerr << '\n';
                }
                matches.emplace_back(children[k], unique_length, std::move(counts));
            }
            else if (debug_gesa && num_nonzero <= 1 && !above_max) {
                cerr << "only occurs on one component\n";
            }
            else if (debug_gesa && !link_more_frequent && !above_max) {
                cerr << "suffix link has same occurrences\n";
            }
        }
    };
    
    
    logging::log(logging::Debug, "Traversing the LCP tree");
    
    // records of (lcp, left, children)
    vector<tuple<size_t, size_t, vector<GESANode>>> stack;
    stack.emplace_back(0, 0, vector<GESANode>());
    for (size_t i = 1; i < lcp_array.size(); ++i) {
        
//        if (debug_gesa) {
//            cerr << "iter " << i << " stack state:\n";
//            for (const auto& rec : stack) {
//                cerr << " (" << get<0>(rec) << "; " << get<1>(rec) << "; ";
//                if (!get<2>(rec).empty()) {
//                    for (size_t j = 0; j < get<2>(rec).size(); ++j) {
//                        if (j) {
//                            cerr << ',';
//                        }
//                        const auto& n = get<2>(rec)[j];
//                        cerr << '[' << n.begin << ',' << n.end << ']';
//                    }
//                }
//                else {
//                    cerr << '.';
//                }
//                cerr << ')';
//            }
//            cerr << '\n';
//        }
        
        GESANode last_node(-1, -1);
        
        // figure out which internal nodes we're leaving
        size_t left = i - 1;
        while (get<0>(stack.back()) > lcp_array[i]) {
            
            auto& top = stack.back();
            
            // emit an internal node
            last_node = GESANode(get<1>(top), i - 1);
            add_matches(last_node, get<2>(top));
            
            left = get<1>(top);
            stack.pop_back();
            if (get<0>(stack.back()) >= lcp_array[i]) {
                // record this as a child of the parent
                get<2>(stack.back()).push_back(last_node);
                last_node = GESANode(-1, -1);
            }
        }
        if (get<0>(stack.back()) < lcp_array[i]) {
            stack.emplace_back(lcp_array[i], left, vector<GESANode>());
            if (last_node != GESANode(-1, -1)) {
                get<2>(stack.back()).push_back(last_node);
            }
        }
    }
    // clear the stack
    while (!stack.empty()) {
        auto& top = stack.back();
        // emit an internal node
        GESANode node(get<1>(top), lcp_array.size() - 1);
        add_matches(node, get<2>(top));
        stack.pop_back();
        // record this as a child of the parent
        if (!stack.empty()) {
            get<2>(stack.back()).push_back(node);
        }
    }
    
    return matches;
}

vector<pair<size_t, vector<uint64_t>>> GESA::walk_matches(const GESANode& node,
                                                          size_t length) const {
    
    if (debug_gesa) {
        cerr << "walking node " << node.begin << " " << node.end << '\n';
    }
    
    vector<pair<size_t, vector<uint64_t>>> matches;
    
    unordered_set<pair<size_t, uint64_t>> match_starts;
    for (size_t i = node.begin; i <= node.end; ++i) {
        
        // get the first node and the component information
        size_t idx = i;
        size_t comp = leaf_to_comp[idx];
        const auto& ranked_ids = component_ranked_ids[comp];
        const auto& nearest_rank = nearest_comp_rank[comp];
        uint64_t node_id = ranked_ids[nearest_rank[idx]];
        
        if (match_starts.count(make_pair(comp, node_id))) {
            // this is a duplicated node
            if (debug_gesa) {
                cerr << "walk at " << i << " has duplicated start " << comp << " " << node_id << '\n';
            }
            continue;
        }
        
        match_starts.emplace(comp, node_id);
        
        matches.emplace_back();
        auto& match = matches.back();
        match.second.reserve(length);
        match.first = comp;
        match.second.push_back(node_id);
                
        // walk the rest of the match
        // FIXME: how should i choose from the combinatorially many identical matches
        // that there might be?
        for (size_t j = 1; j < length; ++j) {
            idx = edges[idx].front();
            match.second.push_back(ranked_ids[nearest_rank[idx]]);
        }
        if (debug_gesa) {
            cerr << "walk from " << i << ", comp " << match.first << ":\n";
            for (auto v : match.second) {
                cerr << '\t' << v << '\n';
            }
        }
    }
    
    return matches;
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
        std::cerr << "computing suffix links\n";
    }
    
    // TODO: decide if I want the leaf suffix links
    
    // records of (lcp, left, right)
    vector<tuple<size_t, size_t, size_t>> stack;
    
    // do a top down traversal to construct lists of l-intervals
    vector<vector<GESANode>> l_interval_lists;
    stack.emplace_back(0, 0, -1);
    for (size_t i = 1; i < lcp_array.size(); ++i) {
        if (debug_gesa) {
            std::cerr << "iteration " << i << ", stack state:\n";
            for (auto x : stack) {
                cerr << " (" << get<0>(x) << ',' << get<1>(x) << ',' << (int64_t) get<2>(x) << ')';
            }
            cerr << '\n';
        }
        
        // figure out which internal nodes we're leaving
        size_t left = i - 1;
        while (get<0>(stack.back()) > lcp_array[i]) {
            auto& top = stack.back();
            get<2>(top) = i - 1;
            if (debug_gesa) {
                std::cerr << "emit internal [" << get<1>(top) << ", " << get<2>(top) << "], l = " << get<0>(top) << '\n';
            }
            while (l_interval_lists.size() <= get<0>(top)) {
                l_interval_lists.emplace_back();
            }
            l_interval_lists[get<0>(top)].emplace_back(get<1>(top), get<2>(top));
            left = get<1>(top);
            stack.pop_back();
        }
        if (lcp_array[i] > get<0>(stack.back())) {
            stack.emplace_back(lcp_array[i], left, -1);
        }
    }
    if (debug_gesa) {
        std::cerr << "clearing the stack\n";
    }
    // clear the stack
    while (!stack.empty()) {
        auto& top = stack.back();
        get<2>(top) = lcp_array.size() - 1;
        if (debug_gesa) {
            std::cerr << "emit internal [" << get<1>(top) << ", " << get<2>(top) << "], l = " << get<0>(top) << '\n';
        }
        while (l_interval_lists.size() <= get<0>(top)) {
            l_interval_lists.emplace_back();
        }
        l_interval_lists[get<0>(top)].emplace_back(get<1>(top), get<2>(top));
        stack.pop_back();
    }
    
    if (debug_gesa) {
        std::cerr << "l interval lists:\n";
        for (size_t l = 0; l < l_interval_lists.size(); ++l) {
            std::cerr << l << ':';
            for (auto node : l_interval_lists[l]) {
                std::cerr << " [" << node.begin << ',' << node.end << "]";
            }
            std::cerr << '\n';
        }
    }
    
    suffix_links.resize(lcp_array.size());
    for (size_t l = 1; l < l_interval_lists.size(); ++l) {
        if (debug_gesa) {
            std::cerr << "forming links for nodes at depth " << l << "\n";
        }
        auto& link_interval_list = l_interval_lists[l - 1];
        for (const GESANode& node : l_interval_lists[l]) {
            const auto& node_edges = edges[node.begin];
            size_t i, j;
            tie(i, j) = st_node_annotation_idx(node);
            if (node_edges.empty()) {
                // this is a sink node, so removing a character gives the empty string
                suffix_links[j] = root();
                if (debug_gesa) {
                    std::cerr << "link [" << node.begin << ',' << node.end << "] -> [" << root().begin << ',' << root().end << "], storing at " << j << "\n";
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
                    cerr << "link [" << node.begin << ',' << node.end << "] -> [" << link_interval_list[lo].begin << ',' << link_interval_list[lo].end << "], storing at " << j << "\n";
                }
                
                // record the suffix link
                suffix_links[j] = link_interval_list[lo];
            }
        }
    }
}


void GESA::label_edges(size_t doubling_steps, const BaseGraph& joined, const PathGraph& path_graph) {
    
    if (debug_gesa) {
        cerr << "labeling edges" << endl;
    }
    
    // the first skip edge is just a normal edge
    vector<vector<uint64_t>> skip_edges(edges.size());
    for (size_t i = 0; i < edges.size(); ++i) {
        const auto& node_edges = edges[i];
        if (!node_edges.empty()) {
            skip_edges[i].push_back(node_edges.front());
        }
    }
    

    // construct skip edge lists
    for (size_t step = 0; step < doubling_steps; ++step) {
        for (size_t i = 0; i < skip_edges.size(); ++i) {
            
            auto& node_skip_edges = skip_edges[i];
            if (node_skip_edges.size() > step) {

                size_t prefix_len = lcp_array[i];
                if (i + 1 < lcp_array.size()) {
                    prefix_len = max(prefix_len, lcp_array[i + 1]);
                }
                if ((1 << step) > prefix_len) {
                    // we don't need to skip any further than the prefix required to
                    // sort this node
                    continue;
                }
                
                uint64_t next = node_skip_edges[step];
                auto& next_skip_edges = skip_edges[next];
                if (next_skip_edges.size() > step) {
                    node_skip_edges.push_back(next_skip_edges[step]);
                }
            }
        }
    }
    
    if (debug_gesa) {
        cerr << "constructed skip edge lists:\n";
        for (size_t i = 0; i < skip_edges.size(); ++i) {
            cerr << i << ": ";
            for (size_t j = 0; j < skip_edges[i].size(); ++j) {
                if (j) {
                    cerr << ", ";
                }
                cerr << skip_edges[i][j];
            }
            cerr << '\n';
        }
    }
    
    // init the annotation
    for (int i : {0, 1}) {
        edge_label[i].resize(lcp_array.size(), -1);
    }
    
    auto add_child_labels = [&](const GESANode& parent) {
        
        // we need to walk this far to traverse this edge
        size_t branch_depth = depth(parent);
        
        // get the order of step sizes corresponding to this depth
        vector<uint32_t> step_sizes;
        for (uint32_t i = 0, n = hi_bit(branch_depth); i <= n; ++i) {
            if (branch_depth & (1 << i)) {
                step_sizes.push_back(i);
            }
        }
        // we want the largest step first
        std::reverse(step_sizes.begin(), step_sizes.end());
        
        if (debug_gesa) {
            cerr << "adding labels to children of " << parent.begin << "," << parent.end << " with branch depth " << branch_depth << " and steps:\n";
            for (auto s : step_sizes) {
                cerr << '\t' << (1 << s) << '\n';
            }
        }
        
        for (GESANode& child : children(parent)) {
            if (debug_gesa) {
                cerr << "labeling child " << child.begin << ',' << child.end << '\n';
            }
            // walk the skip edges to the first
            uint64_t here = child.begin;
            for (size_t step_size : step_sizes) {
                here = skip_edges[here][step_size];
                if (debug_gesa) {
                    cerr << "\tskip " << (1 << step_size) << " to " << here << '\n';
                }
            }
            
            
            size_t i, j;
            tie(i, j) = st_node_annotation_idx(child);
            
            if (debug_gesa) {
                auto base = joined.label(path_graph.from(here));
                if (base <= 4) {
                    base = decode_base(base);
                }
                cerr << "end at base " << base << '\n';
            }
            
            edge_label[i][j] = joined.label(path_graph.from(here));
        }
    };
    
    // records of (lcp, left)
    vector<pair<size_t, size_t>> stack;
    stack.emplace_back(0, 0);
    for (size_t i = 1; i < lcp_array.size(); ++i) {
        
        // figure out which internal nodes we're leaving
        size_t left = i - 1;
        while (stack.back().first > lcp_array[i]) {
            auto& top = stack.back();
            // emit an internal node
            add_child_labels(GESANode(top.second, i - 1));
            left = top.second;
            stack.pop_back();
        }
        if (lcp_array[i] > stack.back().first) {
            stack.emplace_back(lcp_array[i], left);
        }
    }
    // clear the stack
    while (!stack.empty()) {
        auto& top = stack.back();
        // emit an internal node
        add_child_labels(GESANode(top.second, lcp_array.size() - 1));
        stack.pop_back();
    }
}

void GESA::print(ostream& out) const {
    auto label = [](unsigned char c) {
        char l;
        if (c <= 4) {
            l = decode_base(c);
        }
        else if (c == numeric_limits<uint8_t>::max()) {
            l = '.';
        }
        else {
            l = c;
        }
        return l;
    };
    
    out << "i" << '\t' << "Cmp" << '\t' << "Nd" << '\t' << "LCP" << '\t' << "Ch" << '\t' << "SL1" << '\t' << "SL2";
    out << '\t' << "IB" << '\t' << "LB";
    out <<  '\t' << "Es" << '\n';
    for (size_t i = 0; i < lcp_array.size(); ++i) {
        const auto& ranked_id = component_ranked_ids[leaf_to_comp[i]];
        const auto& nearest_rank = nearest_comp_rank[leaf_to_comp[i]];
        out << i << '\t' << leaf_to_comp[i] << '\t' << ranked_id[nearest_rank[i]] << '\t' << lcp_array[i] << '\t' << (i < child_array.size() ? (int) child_array[i] : -1);
        out << " (";
        if (child_array_is_down(i)) {
            out << 'D';
        }
        else if (child_array_is_l_index(i)) {
            out << 'N';
        }
        else if (child_array_is_up(i)) {
            out << 'U';
        }
        else {
            out << '.';
        }
        out << ')';
        if (suffix_links[i] != GESANode()) {
            out << '\t' << suffix_links[i].begin << '\t' << suffix_links[i].end;
        }
        else {
            out << '\t' << -1 << '\t' << -1;
        }
        out << '\t' << label(edge_label[0][i]) << '\t' << label(edge_label[1][i]) << '\t';
        auto& node_edges = edges[i];
        if (node_edges.empty()) {
            out << '.';
        }
        else {
            for (size_t j = 0; j < node_edges.size(); ++j) {
                if (j) {
                    out << ',';
                }
                out << node_edges[j];
            }
        }
        out << '\n';
    }
}

}
