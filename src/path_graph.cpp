#include "path_graph.hpp"

#include <algorithm>

#include "range_min_query.hpp"
#include "utility.hpp"
#include "topological_order.hpp"

namespace centrolign {

using namespace std;

PathGraph::PathGraph(const PathGraph& graph) {
    
    // step 1: do the relational join to generate the new nodes
    doubling_step = graph.doubling_step + 1;
    
    {
        // compute the sorted order by the join variables for a fast merge join
        vector<size_t> order_by_from, order_by_to;
        {
            auto indexes = range_vector(graph.node_size());
            order_by_from = integer_sort(indexes, [&graph](size_t i) { return graph.from(i); });
            order_by_to = integer_sort(indexes, [&graph](size_t i) { return graph.to(i); });
        }
        
        for (size_t i = 0, j = 0; i < order_by_to.size() && j < order_by_from.size(); ) {
            // find range of equal "to" in left relation
            size_t i_end = i + 1;
            while (i_end < order_by_to.size() && graph.to(order_by_to[i_end]) == graph.to(order_by_to[i])) {
                ++i_end;
            }
            
            // walk up to/past matched "from" in right relation
            while (j < order_by_from.size() && graph.from(order_by_from[j]) == graph.to(order_by_to[i])) {
                ++j;
            }
            // find range of equal "from" in right relation
            size_t j_end = j;
            while (j_end < order_by_from.size() && graph.from(order_by_from[j_end]) == graph.to(order_by_from[j])) {
                ++j_end;
            }
            // generate all pairs as nodes
            for (size_t j_inner = j; j_inner < j_end; ++j_inner) {
                for (size_t i_inner = i; i_inner < i_end; ++i_inner) {
                    nodes.emplace_back(graph.from(order_by_to[i_inner]), graph.to(order_by_from[j_inner])
                                       graph.rank(order_by_to[i_inner]), graph.rank(order_by_from[j_inner]));
                }
            }
            
            i = i_end;
            j = j_end;
        }
    }
    
    
    // step 2: convert from pair ranks to integer ranks and merge redundant nodes
    
    // sort the indexes of the new nodes, first by the second rank then by the first
    vector<size_t> indexes = integer_sort(range_vector(nodes.size()), [&](size_t i) { return nodes[i].join_rank; });
    indexes = integer_sort(indexes, [&](size_t i) { return nodes[i].rank; });
    
    // build a RMQ over the previous LCP array to answer general LCP(i, j) queries
    RMQ<size_t> lcp_rmq(graph.lcp_array);
    auto get_lcp = [&](size_t rank1, size_t rank2) {
        if (rank1 == rank2) {
            // they match along the whole 2^k length prefix
            return 1 << graph.doubling_step;
        }
        else {
            // the shortest match between adjacent prefixes is the LCP
            return graph.lcp_array[lcp_rmq.range_arg_min(min(rank1, rank2), max(rank1, rank2))];
        }
    };
    
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
                // the first 2^k bases are the same, the rest is the LCP
                lcp = (1 << graph.doubling_step) + get_lcp(node.join_rank, prev_pre_rank.second);
            }
            else {
                // the match is just the LCP of the first ranks
                lcp = get_lcp(node.rank, prev_pre_rank.first);
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
        if (to_remove[i]) {
            ++removed_so_far;
        }
        else {
            nodes[i - removed_so_far] = nodes[i];
        }
    }
    nodes.resize(nodes.size() - removed_so_far);
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
    
    if (nodes.empty()) {
        return true;
    }
    
    size_t max_rank = 0;
    for (const PathGraphNode& node : nodes) {
        max_rank = max(max_rank, node.rank);
    }
    
    // do all nodes have unique ranks?
    return max_rank + 1 == nodes.size();
}

void PathGraph::order_by_rank() {
    
    // update the edges (if any)
    for (auto& node_edges : edges) {
        for (auto& node_id : node_edges.next) {
            node_id = rank(node_id);
        }
        for (auto& node_id : node_edges.prev) {
            node_id = rank(node_id);
        }
    }
    
    // re-order the nodes
    for (uint64_t i = 0; i < nodes.size(); i++) {
        while (rank(i) != i) {
            std::swap(nodes[rank(i)], nodes[i]);
        }
    }
}

}
