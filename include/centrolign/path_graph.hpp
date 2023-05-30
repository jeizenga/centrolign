#ifndef centrolign_path_graph_hpp
#define centrolign_path_graph_hpp

#include <vector>
#include <cstdint>
#include <utility>
#include <iostream>
#include <stdexcept>
#include <tuple>
#include <functional>

#include "centrolign/utility.hpp"
#include "centrolign/integer_sort.hpp"
#include "centrolign/topological_order.hpp"

namespace centrolign {

class GESA; // forward declaration

/*
 * Graph for intermediate indexing steps representing the product of
 * a phase of prefix doubling
 */
class PathGraph {
public:
    
    // copy from an initial base graph
    template<class BGraph>
    PathGraph(const BGraph& graph);
    // construct by a prefix-doubling step from another PathGraph
    PathGraph(const PathGraph& graph);
    
    PathGraph() = default;
    ~PathGraph() = default;
    
    size_t node_size() const;
    uint64_t from(uint64_t node_id) const;
    uint64_t to(uint64_t node_id) const;
    size_t rank(uint64_t node_id) const;
    
    // LCP (possibly partial if not prefix sorted) between rank and rank + 1
    size_t lcp(size_t rank) const;
    
    bool is_prefix_sorted() const;
    
    // called after prefix sorted to remove redundancy, order by rank,
    // and identify the edges
    template<class BGraph>
    void finish(const BGraph& graph);
    
    // these functions are only valid after calling construct_edges
    const std::vector<uint64_t>& next(uint64_t node_id) const;
    const std::vector<uint64_t>& previous(uint64_t node_id) const;
    size_t next_size(uint64_t node_id) const;
    size_t previous_size(uint64_t node_id) const;
    
private:
    
    static bool debug_path_graph;
    
    
    // reassign node IDs so that they coincide with rank and prefix-range-sort
    // note: only valid if prefix sorted, else causes infinite loop
    void order_by_rank();
    
    void merge_overexpanded_nodes();
    
    void prefix_range_sort();
    
    // construct edges using the original graph you copied from
    template<class BGraph>
    void construct_edges(const BGraph& graph);
    
    
    void print_graph(std::ostream& out) const;
    
    struct PathGraphNode;
    struct PathGraphEdges;
    
    // used to represent the "next node" after a sink
    static uint64_t null_id;
    
    std::vector<PathGraphNode> nodes;
    std::vector<PathGraphEdges> edges;
    
    // LCP array over the shared prefixes of the unique ranks, up to the doubling length
    std::vector<size_t> lcp_array;
    
    size_t doubling_step = 0;
    
    struct PathGraphNode {
        
        PathGraphNode(uint64_t from, uint64_t to) : from(from), to(to) { }
        PathGraphNode(uint64_t from, uint64_t to,
                      size_t rank1, size_t rank2) : from(from), to(to),
                                                    rank(rank1), join_rank(rank2) { }
        PathGraphNode() = default;
        ~PathGraphNode() = default;
        
        uint64_t from = null_id;
        uint64_t to = null_id;
        size_t rank = 0;
        size_t join_rank = 0;
    };
    
    // these are held separately since the nodes frequently exist without edges
    struct PathGraphEdges {
        PathGraphEdges() = default;
        ~PathGraphEdges() = default;
        
        std::vector<uint64_t> next;
        std::vector<uint64_t> prev;
    };
    
    friend class GESA;
};

/*
 * Template implementations
 */

template<class BGraph>
PathGraph::PathGraph(const BGraph& graph) {
    
    // TODO: probably over-engineering the initial sort
    std::vector<size_t> cumul;
    
    nodes.reserve(graph.node_size());
    for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
        if (graph.next_size(node_id) == 0) {
            nodes.emplace_back(node_id, null_id);
        }
        else {
            for (uint64_t next_id : graph.next(node_id)) {
                nodes.emplace_back(node_id, next_id);
            }
        }
        
        // record this base as 'seen'
        char base = graph.label(node_id);
        while (cumul.size() <= base) {
            cumul.emplace_back(0);
        }
        cumul[base] = 1;
    }
    
    // convert the 'seen' vector to a cumulative sum vector
    for (size_t i = 1; i < cumul.size(); ++i) {
        cumul[i] += cumul[i - 1];
    }
    // remove offset to convert to 0-based ranks
    for (size_t i = 0; i < cumul.size(); ++i) {
        // note: this goes to -1 for 0s, but we never will read them anyway
        --cumul[i];
    }
    
    // give rank to node based on its base
    for (PathGraphNode& node : nodes) {
        node.rank = cumul[graph.label(node.from)];
    }
    
    if (debug_path_graph) {
        std::cerr << "initialized graph, current state:\n";
        print_graph(std::cerr);
    }
    
    // merge ranks
    // TODO: largely copied from doubling step, could merge implementations...
    {
        auto indexes = integer_sort(range_vector(nodes.size()), [&](size_t i) {
            return nodes[i].rank;
        });

        std::vector<bool> remove(nodes.size(), false);
        for (size_t i = 0; i < indexes.size();) {
            // find the range of equal-ranked nodes
            size_t j = i + 1;
            bool shared_from = true;
            while (j < indexes.size() && nodes[indexes[j]].rank == nodes[indexes[i]].rank) {
                shared_from = shared_from && (nodes[indexes[j]].from == nodes[indexes[i]].from);
                ++j;
            }
            if (shared_from) {
                // the nodes are redundant and all but one can be removed
                for (size_t k = i + 1; k < j; ++k) {
                    remove[indexes[k]] = true;
                }
            }
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
    }

    if (debug_path_graph) {
        std::cerr << "merged nodes, path graph state:\n";
        print_graph(std::cerr);
    }
    
    // at this point ranks correspond to distinct characters, so the LCP between
    // distinct ranks is always 0
    lcp_array.resize(cumul.back(), 0);
}

template<class BGraph>
void PathGraph::finish(const BGraph& graph) {
    order_by_rank();
    merge_overexpanded_nodes();
    construct_edges(graph);
}

template<class BGraph>
void PathGraph::construct_edges(const BGraph& graph) {
    
    // get incomplete edges as (parent_prev_id, child_next_id)
    std::vector<std::pair<uint64_t, uint64_t>> pre_edges;
    pre_edges.reserve(node_size()); // heuristic to speed up allocation
    for (uint64_t node_id = 0; node_id < node_size(); ++node_id) {
        for (uint64_t parent_prev_id : graph.previous(from(node_id))) {
            
            pre_edges.emplace_back(parent_prev_id, node_id);
        }
    }
    
    // get the edges sorted into blocks sharing the same first node
    std::vector<size_t> indexes = integer_sort(range_vector(pre_edges.size()),
                                               [&](size_t i) { return this->rank(pre_edges[i].second); });
    indexes = integer_sort(indexes, [&](size_t i) { return graph.label(pre_edges[i].first); });
    
    if (debug_path_graph) {
        std::cerr << "pre-edges:\n";
        size_t j = 0;
        for (auto i : indexes) {
            std::cerr << j++ << ": " << pre_edges[i].first << '\t' << pre_edges[i].second << '\n';
        }
    }
    
    // records of (node id begin, node id end, pre edge index begin, pre edge index end)
    std::vector<std::tuple<uint64_t, uint64_t, size_t, size_t>> unresolved_intervals;
    
    // construct the adjacency lists (assumes rank ordering for node IDs)
    edges.resize(node_size());
    for (uint64_t node_id = 0, i = 0; node_id < node_size();) {
        
        // find the interval of path node IDs that match the current "from"
        uint64_t node_id_end = node_id + 1;
        while (node_id_end < node_size() && from(node_id_end) == from(node_id)) {
            ++node_id_end;
        }
        // find the interval of edge indexes that match the current "from"
        size_t j = i;
        while (j < indexes.size() && pre_edges[indexes[j]].first == from(node_id)) {
            ++j;
        }
        
        if (node_id_end == node_id + 1 || i == j) {
            // there is a unique "from" node, so there is no confusion over which
            // path node should get the edges
            
            // iterate through edges that come from this node
            auto& edge_lists = edges[node_id];
            for (; i < j; ++i) {
                // create an edge
                uint64_t next_id = pre_edges[indexes[i]].second;
                edge_lists.next.emplace_back(next_id);
                edges[next_id].prev.emplace_back(node_id);
                if (debug_path_graph) {
                    std::cerr << "create edge " << node_id << " -> " << next_id << '\n';
                }
            }
        }
        else if (node_id_end - node_id == j - i) {
            // every node except for a sink sentinel must have an outgoing edge, and
            // a sink sentinel can't end up on this path because it is fully prefix sorted
            // by its own character. that let's us assign the edges one-to-one with nodes
            // by the pigeonhole principle
            if (debug_path_graph) {
                std::cerr << "trivially subdividing ambiguous interval " << node_id << ":" << node_id_end << " and " << i << ":" << j << '\n';
            }
            while (node_id != node_id_end) {
                uint64_t next_id = pre_edges[indexes[i]].second;
                edges[node_id].next.emplace_back(next_id);
                edges[next_id].prev.emplace_back(node_id);
                if (debug_path_graph) {
                    std::cerr << "create edge " << node_id << " -> " << next_id << '\n';
                }
                ++node_id;
                ++i;
            }
        }
        else {
            // we cannot immediately tell which path node should get the edges
            unresolved_intervals.emplace_back(node_id, node_id_end, i, j);
        }
        
        node_id = node_id_end;
        i = j;
    }
    
    
    
    if (debug_path_graph) {
        std::cerr << "failed to resolve " << unresolved_intervals.size() << " intervals:\n";
        for (auto& interval : unresolved_intervals) {
            uint64_t path_node_begin, path_node_end;
            size_t idx_begin, idx_end;
            std::tie(path_node_begin, path_node_end, idx_begin, idx_end) = interval;
            
            
            std::cerr << '\t' << (path_node_end - path_node_begin) << " nodes " << path_node_begin << ":" << path_node_end << ", and " << (idx_end - idx_begin) << " edges " << idx_begin << ":" << idx_end << '\n';
            std::cerr << "\tlcps:\n";
            for (auto n = path_node_begin; n + 1 < path_node_end; ++n) {
                std::cerr << "\t\t" << lcp_array[n] << '\n';
            }
            std::cerr << "\ttargets:\n";
            for (auto i = idx_begin; i < idx_end; ++i) {
                std::cerr << "\t\t" << pre_edges[i].second << '\n';
            }
        }
    }
    
    if (!unresolved_intervals.empty()) {
        // we need to walk out some edges to disambiguate where they come from
        
        if (debug_path_graph) {
            std::cerr << "succesfully completed edges:\n";
            for (uint64_t node_id = 0; node_id < node_size(); ++node_id) {
                std::cerr << node_id << ":\n";
                for (auto next_id : next(node_id)) {
                    std::cerr << "\t-> " << next_id << '\n';
                }
            }
        }
        
        // we will lazily construct these
        std::vector<std::vector<uint64_t>> skip_edges(node_size());
        
        // lazily build skip edges and return a node that is 2^skip_power steps away
        // from the indicated node
        std::function<uint64_t(uint64_t, uint32_t)> get_skip = [&](uint64_t node_id, uint32_t skip_power) {
            auto& node_skip_edges = skip_edges[node_id];
            if (skip_power < node_skip_edges.size()) {
                // we've already constructed this skip edge, return it
                return node_skip_edges[skip_power];
            }
            else if (skip_power == 0) {
                // base case, a single edge
                node_skip_edges.push_back(edges[node_id].next.front());
                return node_skip_edges.front();
            }
            else {
                // build up earlier skips
                while (node_skip_edges.size() < skip_power) {
                    get_skip(node_id, node_skip_edges.size());
                }
                // take the first hop
                auto next_id = node_skip_edges[skip_power - 1];
                // take the second hop and memoize it
                node_skip_edges.push_back(get_skip(next_id, skip_power - 1));
                return node_skip_edges.back();
            }
        };
        
        // use skip edges to walk some path through the graph of the given length
        // in log(length) time, assumes that all necessary skip edges are present
        auto skip_walk = [&](uint64_t node_id, size_t length) {
            // figure out which size steps we need to take
            std::vector<uint32_t> step_sizes;
            for (uint32_t i = 0, n = hi_bit(length); i <= n; ++i) {
                if (length & (1 << i)) {
                    step_sizes.push_back(i);
                }
            }
            // walk the skip edges to the first
            uint64_t here = node_id;
            for (size_t step_size : ReverseForEachAdapter<std::vector<uint32_t>>(step_sizes)) {
                here = get_skip(here, step_size);
            }
            return here;
        };
        
        // we resolve the edges in reverse topological order of the original graph, since the last one
        // must be resolvable by the currently constructed skip edges (no edges further than it are
        // missing)
        // TODO: is this reasoning solid? can some of the other nodes in the interval be missing skip edges?
        std::vector<std::vector<size_t>> queue(graph.node_size());
        {
            auto top_index = invert(topological_order(graph));
            for (size_t i = 0; i < unresolved_intervals.size(); ++i) {
                uint64_t path_node_begin, path_node_end;
                size_t idx_begin, idx_end;
                std::tie(path_node_begin, path_node_end, idx_begin, idx_end) = unresolved_intervals[i];
                
                size_t max_idx = 0;
                uint64_t which_max_id = -1; // this should get optimized away without debug
                for (uint64_t node_id = path_node_begin; node_id < path_node_end; ++node_id) {
                    auto idx = top_index[from(node_id)];
                    if (idx >= max_idx) {
                        max_idx = idx;
                        which_max_id = node_id;
                    }
                    max_idx = std::max<size_t>(max_idx, top_index[from(node_id)]);
                }
                if (debug_path_graph) {
                    std::cerr << "interval " << i << " has max topological index of " << max_idx << ", achieved by " << which_max_id << " with from() value " << from(which_max_id) << '\n';
                }
                
                queue[max_idx].push_back(i);
            }
        }
        for (int64_t k = queue.size() - 1; k >= 0; --k) {
            for (auto j : queue[k]) {
                uint64_t path_node_begin, path_node_end;
                size_t idx_begin, idx_end;
                std::tie(path_node_begin, path_node_end, idx_begin, idx_end) = unresolved_intervals[j];
                
                if (debug_path_graph) {
                    std::cerr << "directly resolving intervals " << path_node_begin << ":" << path_node_end << " and " <<  idx_begin << ":" << idx_end << " from topo index bucket " << k << '\n';
                }
                
                uint64_t curr_node = path_node_begin;
                
                for (size_t i = idx_begin; i < idx_end; ++i) {
                    if (debug_path_graph) {
                        std::cerr << "\tassigning edge " << i << '\n';
                    }
                    // the node we will link to
                    uint64_t tail_id = pre_edges[indexes[i]].second;
                    
                    if (idx_end - i < path_node_end - curr_node) {
                        // one edge each for the rest of the nodes by pigeonhole principle
                        
                        if (debug_path_graph) {
                            std::cerr << "\tassign by pigeonhole principle" << '\n';
                        }
                        ++curr_node;
                    }
                    else if (i != idx_begin                       // first edge always goes to first node
                             && curr_node + 1 != path_node_end) { // after moving to last node, it gets all the edges
                        
                        // this is a case where we actually just need to check whether the current edge belongs
                        // to the current node
                        
                        // walk forward to the first mismatch with the next node on this and the previous edge
                        size_t walk_len = lcp_array[curr_node] - 1;
                        if (debug_path_graph) {
                            std::cerr << "\twalking " << walk_len << " from adjacent edge tails " << pre_edges[indexes[i - 1]].second << " and " << tail_id << '\n';
                        }
                        uint64_t prev_walked_id = skip_walk(pre_edges[indexes[i - 1]].second, walk_len);
                        uint64_t curr_walked_id = skip_walk(tail_id, walk_len);
                        if (debug_path_graph) {
                            auto c1 = graph.label(from(prev_walked_id));
                            auto c2 = graph.label(from(curr_walked_id));
                            if (c1 <= 4) {
                                c1 = decode_base(c1);
                            }
                            if (c2 <= 4) {
                                c2 = decode_base(c2);
                            }
                            std::cerr << "\tgot nodes " << prev_walked_id << " and " << curr_walked_id << " with labels " << c1 << " and " << c2 << '\n';
                        }
                        if (graph.label(from(prev_walked_id)) != graph.label(from(curr_walked_id))) {
                            // we found a mismatch where we would expect to at the boundary of the
                            // two nodes' edges, so we've actually moved into the next nodes block of edges
                            ++curr_node;
                        }
                    }
                    
                    // create the edge
                    edges[curr_node].next.emplace_back(tail_id);
                    edges[tail_id].prev.emplace_back(curr_node);
                    if (debug_path_graph) {
                        std::cerr << "\tcreate edge " << curr_node << " -> " << tail_id << '\n';
                    }
                }
                // make sure that we cleared the entire list of nodes
                assert(curr_node + 1 == path_node_end);
            }
        }
    }
    
    if (debug_path_graph) {
        std::cerr << "constructed edges:\n";
        print_graph(std::cerr);
    }
}

}

#endif /* centrolign_path_graph_hpp */
