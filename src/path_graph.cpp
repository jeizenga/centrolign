#include "centrolign/path_graph.hpp"

#include <algorithm>
#include <tuple>
#include <limits>

#include "centrolign/range_min_query.hpp"
#include "centrolign/utility.hpp"

namespace centrolign {

using namespace std;

PathGraphSizeException::PathGraphSizeException(const PathGraph& current_graph,
                                               const PathGraph& previous_graph, size_t step) noexcept
    : curr_count(make_from_count(current_graph)), prev_count(make_from_count(previous_graph)), step(step)
{
    while (curr_count.size() < prev_count.size()) {
        curr_count.push_back(0);
    }
}

std::vector<uint64_t> PathGraphSizeException::make_from_count(const PathGraph& path_graph) noexcept {
    std::vector<uint64_t> count;
    for (uint64_t node_id = 0; node_id < path_graph.node_size(); ++node_id) {
        while (count.size() <= path_graph.from(node_id)) {
            count.push_back(0);
        }
        ++count[path_graph.from(node_id)];
    }
    return count;
}


const char* PathGraphSizeException::what() const noexcept {
    return "Exceeded PathGraph size limit";
}

const vector<uint64_t>& PathGraphSizeException::from_count() const {
    return curr_count;
}

const vector<uint64_t>& PathGraphSizeException::previous_from_count() const {
    return prev_count;
}

size_t PathGraphSizeException::doubling_step() const {
    return step;
}

const bool PathGraph::debug_path_graph = false;
const bool PathGraph::instrument_path_graph = true;
const size_t PathGraph::from_doubling_step = 13;
const size_t PathGraph::min_count = 8;

uint64_t PathGraph::null_id = -1;

PathGraph::PathGraph(const PathGraph& graph, size_t size_limit) {
    
    // step 1: do the relational join to generate the new nodes
    doubling_step = graph.doubling_step + 1;
    
    {
        // compute the sorted order by the join variables for a fast merge join
        vector<size_t> order_by_from, order_by_to;
        {
            auto indexes = range_vector(graph.node_size());
            order_by_from = integer_sort(indexes, [&graph](size_t i) { return graph.from(i); });
            // note: we adjust by +1 to wrap the -1 null ID into positive integers
            order_by_to = integer_sort(indexes, [&graph](size_t i) { return graph.to(i) + 1; });
            if (debug_path_graph) {
                cerr << "order by to:\n";
                for (auto i : order_by_to) {
                    cerr << ' ' << i;
                }
                cerr << '\n';
                cerr << "order by from:\n";
                for (auto i : order_by_from) {
                    cerr << ' ' << i;
                }
                cerr << '\n';
            }
        }
        
        // we only need to know uniqueness, so count up the occurrences of each rank up to 2
        vector<uint8_t> rank_count;
        for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
            size_t rank = graph.rank(node_id);
            while (rank_count.size() <= rank) {
                rank_count.push_back(0);
            }
            rank_count[rank] = min(rank_count[rank] + 1, 2);
        }
        
        for (size_t i = 0, j = 0; i < order_by_to.size(); ) {
            
            // find range of equal "to" in left relation
            size_t i_end = i + 1;
            while (i_end < order_by_to.size() && graph.to(order_by_to[i_end]) == graph.to(order_by_to[i])) {
                ++i_end;
            }
            
            // walk up to/past matched "from" in right relation
            // note: we have to use the +1 offset to make the -1 null id wrap to the beginning
            while (j < order_by_from.size() && graph.from(order_by_from[j]) + 1 < graph.to(order_by_to[i]) + 1) {
                ++j;
            }
            // find range of equal "from" in right relation
            size_t j_end = j;
            while (j_end < order_by_from.size() && graph.from(order_by_from[j_end]) == graph.to(order_by_to[i])) {
                ++j_end;
            }
            
            if (debug_path_graph) {
                cerr << "match ranges " << i << ":" << i_end << " and " << j << ":" << j_end << " on to " << (int) graph.to(order_by_to[i]) << ", from " << (j < order_by_from.size() ? (int) order_by_from[j] : -1) << '\n';
                for (size_t i_inner = i; i_inner < i_end; ++i_inner) {
                    if (graph.rank(order_by_to[i_inner]) != 1) {
                        cerr << ' ' << order_by_to[i_inner];
                    }
                }
                cerr << '\n';
                if (j < order_by_from.size()) {
                    for (size_t j_inner = j; j_inner < j_end; ++j_inner) {
                        cerr << ' ' << order_by_from[j_inner];
                    }
                }
                cerr << '\n';
            }
            
            // generate all pairs as nodes
            for (size_t i_inner = i; i_inner < i_end; ++i_inner) {
                uint64_t node_id = order_by_to[i_inner];
                if (rank_count[graph.rank(node_id)] == 1) {
                    // this node has a unique rank, can be copied directly
                    if (rank_count[graph.rank(node_id)] == 1) {
                        nodes.emplace_back(graph.from(node_id), graph.to(node_id),
                                           graph.rank(node_id), 0);
                        // check if we broke the limit
                        if (nodes.size() > size_limit) {
                            throw PathGraphSizeException(*this, graph, doubling_step);
                        }
                    }
                }
                else {
                    // generate all pairs to double the prefix length for this node
                    for (size_t j_inner = j; j_inner < j_end; ++j_inner) {
                        nodes.emplace_back(graph.from(node_id), graph.to(order_by_from[j_inner]),
                                           graph.rank(node_id), graph.rank(order_by_from[j_inner]));
                        // check if we broke the limit
                        if (nodes.size() > size_limit) {
                            throw PathGraphSizeException(*this, graph, doubling_step);
                        }
                    }
                    continue;
                }
            }
            
            i = i_end;
            j = j_end;
        }
    }
    
    if (debug_path_graph) {
        cerr << "peformed relational join, current graph state:\n";
        print_graph(cerr);
    }
    
    // step 2: convert from pair ranks to integer ranks and merge redundant nodes
    
    // sort the indexes of the new nodes, first by the second rank then by the first
    vector<size_t> indexes = integer_sort(range_vector(nodes.size()), [&](size_t i) { return nodes[i].join_rank; });
    indexes = integer_sort(indexes, [&](size_t i) { return nodes[i].rank; });
    
    // build an RMQ over the previous LCP array to answer general LCP(i, j) queries
    RMQ<size_t> lcp_rmq(graph.lcp_array);
    
    vector<bool> remove(nodes.size(), false);
    size_t next_rank = 0;
    pair<size_t, size_t> prev_pre_rank(0, 0);
    for (size_t i = 0; i < indexes.size();) {
        // find the range of equal-ranked nodes
        size_t j = i + 1;
        while (j < indexes.size() && nodes[indexes[i]].rank == nodes[indexes[j]].rank
               && nodes[indexes[i]].join_rank == nodes[indexes[j]].join_rank) {
            ++j;
        }
        
        // figure out the LCP between this rank and the previous rank
        if (next_rank != 0) {
            auto& node = nodes[indexes[i]];
            size_t lcp;
            if (node.rank == prev_pre_rank.first) {
                if (node.join_rank == prev_pre_rank.second) {
                    // match over the full prefix out to this doubling length
                    lcp = (1 << doubling_step);
                }
                else {
                    // full match over the first half of the prefix, followed by LCP over second half
                    lcp = ((1 << graph.doubling_step) +
                           graph.lcp_array[lcp_rmq.range_arg_min(min(node.join_rank, prev_pre_rank.second),
                                                                 max(node.join_rank, prev_pre_rank.second))]);
                }
            }
            else {
                // LCP over the first half of the prefix
                lcp = graph.lcp_array[lcp_rmq.range_arg_min(min(node.rank, prev_pre_rank.first),
                                                            max(node.rank, prev_pre_rank.first))];
            }
            lcp_array.push_back(lcp);
        }
        
        // remember the pair rank for the next iteration
        prev_pre_rank.first = nodes[indexes[i]].rank;
        prev_pre_rank.second = nodes[indexes[i]].join_rank;
        
        // assign this entire range the next rank
        bool shared_from = true;
        for (size_t k = i; k < j; ++k) {
            auto& node = nodes[indexes[k]];
            node.rank = next_rank;
            node.join_rank = 0;
            shared_from = shared_from && (node.from == nodes[indexes[i]].from);
        }
        ++next_rank;
        
        if (shared_from) {
            // the nodes are redundant and all but one can be removed
            for (size_t k = i + 1; k < j; ++k) {
                remove[indexes[k]] = true;
            }
        }
        // advance to next range
        i = j;
    }
    
    // filter out the nodes that we merged
    size_t removed_so_far = 0;
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (remove[i]) {
            ++removed_so_far;
        }
        else {
            nodes[i - removed_so_far] = nodes[i];
        }
    }
    nodes.resize(nodes.size() - removed_so_far);
    
    if (debug_path_graph) {
        cerr << "merged redundant nodes, current graph state:\n";
        print_graph(cerr);
        
        // just to make this easier to debug
        auto order = integer_sort(range_vector(nodes.size()), [&](size_t i) { return nodes[i].rank; });
        auto idx = invert(order);
        for (size_t i = 0; i < idx.size(); ++i) {
            while (idx[i] != i) {
                std::swap(nodes[i], nodes[idx[i]]);
                std::swap(idx[i], idx[idx[i]]);
            }
        }
        
        cerr << "reorder nodes, current graph state:\n";
        print_graph(cerr);
    }
    
    if (instrument_path_graph) {
        if (doubling_step >= from_doubling_step) {
            
            vector<size_t> from_count;
            for (uint64_t node_id = 0; node_id < node_size(); ++node_id) {
                while (from_count.size() <= from(node_id)) {
                    from_count.push_back(0);
                }
                ++from_count[from(node_id)];
            }
            vector<size_t> histogram(*max_element(from_count.begin(), from_count.end()) + 1, 0);
            for (size_t i = 0; i < from_count.size(); ++i) {
                ++histogram[from_count[i]];
            }
            
            cerr << "node size: " << node_size() << '\n';
            cerr << "high from-count histogram (min " << min_count << "):\n";
            for (size_t i = min_count; i < histogram.size(); ++i) {
                if (histogram[i] != 0) {
                    cerr << i << ": " << histogram[i] << '\n';
                }
            }
        }
    }
}

uint64_t PathGraph::from(uint64_t node_id) const {
    return nodes[node_id].from;
}

uint64_t PathGraph::to(uint64_t node_id) const {
    return nodes[node_id].to;
}

size_t PathGraph::rank(uint64_t node_id) const {
    return nodes[node_id].rank;
}

size_t PathGraph::lcp(size_t rank) const {
    return lcp_array[rank];
}

size_t PathGraph::node_size() const {
    return nodes.size();
}

const vector<uint64_t>& PathGraph::next(uint64_t node_id) const {
    return edges[node_id].next;
}

const vector<uint64_t>& PathGraph::previous(uint64_t node_id) const {
    return edges[node_id].prev;
}

size_t PathGraph::next_size(uint64_t node_id) const {
    return edges[node_id].next.size();
}

size_t PathGraph::previous_size(uint64_t node_id) const {
    return edges[node_id].prev.size();
}

bool PathGraph::is_prefix_sorted() const {
    // do all nodes have unique ranks?
    return nodes.empty() || lcp_array.size() + 1 == nodes.size();
}

void PathGraph::order_by_rank() {
    
    // re-order the nodes (nodes shouldn't exist yet)
    for (uint64_t i = 0; i < nodes.size(); i++) {
        while (rank(i) != i) {
            std::swap(nodes[rank(i)], nodes[i]);
        }
    }
    
    if (debug_path_graph) {
        cerr << "reordered nodes by rank\n";
        print_graph(cerr);
    }
}

void PathGraph::prefix_range_sort() {
    
    if (debug_path_graph) {
        cerr << "prefix range sorting\n";
    }
    
    size_t removed_total = 0;
    for (uint64_t i = 0, j = 0; i < nodes.size(); i = j) {
        
        // figure out how many to remove this round
        size_t removed = 0;
        j = i + 1;
        while (j < nodes.size() && from(j) == from(i)) {
            if (debug_path_graph) {
                cerr << "remove " << j << " in favor of " << i << "\n";
            }
            if (j < lcp_array.size()) {
                // accumulate LCP information in the corresponding overlap index
                lcp_array[i] = min(lcp_array[i], lcp_array[j]);
            }
            ++j;
            ++removed;
        }
        
        // move the nodes and lcp entries into the prefix
        nodes[i - removed_total] = nodes[i];
        if (i < lcp_array.size()) {
            lcp_array[i - removed_total] = lcp_array[i];
        }
        removed_total += removed;
    }
    nodes.resize(nodes.size() - removed_total);
    lcp_array.resize(lcp_array.size() - removed_total);
    
    if (debug_path_graph) {
        cerr << "obtained graph:\n";
        print_graph(cerr);
    }
}

void PathGraph::merge_overexpanded_nodes() {
    
    if (debug_path_graph) {
        cerr << "looking for over-expanded nodes\n";
    }
    
    // find the subtrees of the suffix tree that all have the same from()
    // value, which can then be merged
    
    using frame_t = tuple<int64_t, int64_t, int64_t, vector<pair<int64_t, int64_t>>, bool, uint64_t>;
    
    static const frame_t null(-1, -1, -1, vector<pair<int64_t, int64_t>>(), false, -1);
    
    vector<pair<size_t, size_t>> to_merge;
    
    // records of (lcp, lb, rb, children, all from are equal, from)
    vector<frame_t> stack;
    
    auto process = [&](frame_t& frame) {
        if (debug_path_graph) {
            cerr << "processing LCP interval " << get<1>(frame) << ":" << get<2>(frame) << '\n';
        }
        if (!get<4>(frame)) {
            // we've already identified non-uniformities
            return;
        }
        if (get<5>(frame) == -1) {
            // we haven't communicated up from any children, decide what the from value will be here
            get<5>(frame) = from(get<1>(frame));
            if (debug_path_graph) {
                cerr << "all children are leaves, choosing " << from(get<1>(frame)) << " as from\n";
            }
        }
        
        // check the leaf children to see if they match the from value
        auto& children = get<3>(frame);
        for (int64_t i = 0, n = children.size() + 1; i < n; ++i) {
            int64_t begin = i == 0 ? get<1>(frame) : children[i - 1].second + 1;
            int64_t end = i == children.size() ? get<2>(frame) + 1 : children[i].first;
            for (size_t j = begin; j < end; ++j) {
                get<4>(frame) = get<4>(frame) && from(j) == get<5>(frame);
                if (debug_path_graph) {
                    cerr << '\t' << "leaf child " << j << " has from " << from(j) << ", still equivalent? " << get<4>(frame) << '\n';
                }
            }
        }
        
        if (get<4>(frame)) {
            // all of the children of this node have the same from value
            
            if (debug_path_graph) {
                cerr << "identified " << get<1>(frame) << "," << get<2>(frame) << " as mergeable\n";
            }
            
            // clear out children of this node from the return value
            while (!to_merge.empty()
                   && get<1>(frame) <= to_merge.back().first
                   && get<2>(frame) >= to_merge.back().second) {
                if (debug_path_graph) {
                    cerr << "clearing child " << to_merge.back().first << "," << to_merge.back().second << '\n';
                }
                to_merge.pop_back();
            }
            
            to_merge.emplace_back(get<1>(frame), get<2>(frame));
        }
    };
    
    auto communicate_to_parent = [&](frame_t& frame,
                                     frame_t& parent_frame) {
        if (debug_path_graph) {
            cerr << "communicating LCP interval " << get<1>(frame) << ":" << get<2>(frame) << " to parent " << get<1>(parent_frame) << ":" << get<2>(parent_frame) << '\n';
        }
        if (get<5>(parent_frame) == -1) {
            // this is the first child that has communicated a from to the parent
            get<5>(parent_frame) = get<5>(frame);
            if (debug_path_graph) {
                cerr << "\tpassing from value " << get<5>(parent_frame) << " to parent\n";
            }
        }
        
        // check that the child from from agrees with the parent from
        get<4>(parent_frame) = (get<4>(parent_frame) && get<4>(frame)
                                && get<5>(frame) == get<5>(parent_frame));
        
        if (debug_path_graph) {
            cerr << "\tstill consistent? " << get<4>(parent_frame) << '\n';
        }
        
        // remember it as a child
        get<3>(parent_frame).emplace_back(get<1>(frame), get<2>(frame));
    };
    
    // top down traversal of the LCP interval tree
    auto last_frame = null;
    stack.emplace_back(0, 0, 0, vector<pair<int64_t, int64_t>>(), true, -1);
    for (size_t i = 0; i < lcp_array.size(); ++i) {
        
        int64_t lb = i;
        while (get<0>(stack.back()) > lcp_array[i]) {
            // LCPs are falling, must be ending LCP intervals
            get<2>(stack.back()) = i;
            
            last_frame = move(stack.back());
            stack.pop_back();
            process(last_frame);
            
            lb = get<1>(last_frame);
            if (get<0>(stack.back()) >= lcp_array[i]) {
                communicate_to_parent(last_frame, stack.back());
                last_frame = null;
            }
        }
        if (get<0>(stack.back()) < lcp_array[i]) {
            // LCP is increasing, indicates start of LCP intervals
            stack.emplace_back(lcp_array[i], lb, -1, vector<pair<int64_t, int64_t>>(), true, -1);
            if (last_frame != null) {
                communicate_to_parent(last_frame, stack.back());
                last_frame = null;
            }
        }
    }
    
    // clear the satck
    while (!stack.empty()) {
        get<2>(stack.back()) = lcp_array.size();
        
        last_frame = move(stack.back());
        stack.pop_back();
        
        process(last_frame);
        
        if (!stack.empty()) {
            communicate_to_parent(last_frame, stack.back());
        }
    }
    
    if (to_merge.empty()) {
        // we didn't find any subtrees to remove
        return;
    }
    
    if (debug_path_graph) {
        cerr << "found mergeable LCP intervals:\n";
        for (auto interval : to_merge) {
            cerr << '\t' << interval.first << ":" << interval.second << '\n';
        }
    }
    
    // note: nodes before the first merge interval are already where we want them
    size_t removed = 0;
    for (size_t i = 0; i < to_merge.size(); ++i) {
        // remove all but one (note closed interval arithmetic)
        removed += to_merge[i].second - to_merge[i].first;
        
        // keep the final node of this previous interval
        int64_t begin = to_merge[i].second;
        // stop at the start of the next interval (or the end)
        int64_t end = i + 1 == to_merge.size() ? nodes.size() : to_merge[i + 1].first;
        
        // note: because the LCP interval has higher LCP internal than to its neighbors, the
        // LCP onto it from above "falls through", and the final LCP out of it applies to
        // whichever node from the interval that we keep
        
        for (size_t j = begin; j < end; ++j) {
            nodes[j - removed] = nodes[j];
            if (j < lcp_array.size()) {
                lcp_array[j - removed] = lcp_array[j];
            }
        }
    }
    nodes.resize(nodes.size() - removed);
    lcp_array.resize(lcp_array.size() - removed);
    
    if (debug_path_graph) {
        cerr << "finished merging:\n";
        print_graph(cerr);
    }
}

void PathGraph::print_graph(std::ostream& out) const {
    for (size_t i = 0; i < nodes.size(); ++i) {
        const auto& node = nodes[i];
        out << i << ": {from=" << (int64_t) node.from << ", to=" << (int64_t) node.to << ", r1=" << node.rank << ", r2=" << node.join_rank << "}\n";
        if (!edges.empty()) {
            for (auto next : edges[i].next) {
                out << "\t-> " << next << '\n';
            }
        }
    }
    out << "LCP:\n";
    for (size_t i = 0; i < lcp_array.size(); ++i) {
        out << i << ": " << lcp_array[i] << '\n';
    }
}

}
