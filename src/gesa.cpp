#include <numeric>
#include <limits>
#include <unordered_set>

#include "centrolign/gesa.hpp"
#include "centrolign/range_unique_query.hpp"

namespace centrolign {

using namespace std;

bool GESA::debug_gesa = false;

const char* GESASizeException::what() const noexcept {
    return "Size limit exceeded while constructing GESA";
}

const std::vector<std::vector<uint64_t>>& GESASizeException::from_counts() const {
    return curr_counts;
}

const std::vector<std::vector<uint64_t>>& GESASizeException::previous_from_counts() const {
    return prev_counts;
}

size_t GESASizeException::doubling_step() const {
    return step;
}

GESASizeException::GESASizeException(const PathGraphSizeException& ex,
                                     const std::vector<uint16_t>& node_to_comp,
                                     const std::vector<size_t>& component_ranges) noexcept {
    
    step = ex.doubling_step();
    
    curr_counts.resize(component_ranges.size() - 1);
    prev_counts.resize(component_ranges.size() - 1);
    
    assert(ex.from_count().size() == ex.previous_from_count().size());
    
    for (uint64_t node_id = 0; node_id < ex.from_count().size(); ++node_id) {
        auto comp = node_to_comp[node_id];
        auto orig_id = node_id - component_ranges[comp];
        while (curr_counts[comp].size() <= orig_id) {
            curr_counts[comp].emplace_back(0);
            prev_counts[comp].emplace_back(0);
        }
        curr_counts[comp][orig_id] = ex.from_count()[node_id];
        prev_counts[comp][orig_id] = ex.previous_from_count()[node_id];
    }
    
}


vector<SANode> GESA::children(const SANode& parent) const {
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

vector<tuple<SANode, size_t, vector<uint64_t>>> GESA::minimal_rare_matches(size_t max_count) const {
    auto label_getter = [&](const SANode& parent, const SANode& child) {
        return label(child);
    };
    return ESA::minimal_rare_matches_internal(max_count, label_getter);
}

vector<pair<size_t, vector<uint64_t>>> GESA::walk_matches(const SANode& node,
                                                          size_t length) const {
    auto advance = [&](size_t node) -> size_t {
        return edges[node].front();
    };
    return ESA::walk_matches_internal(node, length, advance);
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
    
    auto add_child_labels = [&](const SANode& parent) {
        
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
        
        for (SANode& child : children(parent)) {
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
            add_child_labels(SANode(top.second, i - 1));
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
        add_child_labels(SANode(top.second, lcp_array.size() - 1));
        stack.pop_back();
    }
}

size_t GESA::memory_size() const {
    size_t edge_size = 0;
    for (const auto& edge_list : edges) {
        edge_size += edge_list.capacity();
    }
    return (esa_struct_size()
            + sizeof(edges) + edges.capacity() * sizeof(decltype(edges)::value_type)
            + edge_size * sizeof(decltype(edges)::value_type::value_type));
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
        if (suffix_links[i] != SANode()) {
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
