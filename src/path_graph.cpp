#include "centrolign/path_graph.hpp"

#include <algorithm>

#include "centrolign/range_min_query.hpp"
#include "centrolign/utility.hpp"

namespace centrolign {

using namespace std;

bool PathGraph::debug_path_graph = false;

uint64_t PathGraph::null_id = -1;

PathGraph::PathGraph(const PathGraph& graph) {
    
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
                        if (graph.rank(order_by_from[j_inner]) != 1) {
                            cerr << ' ' << order_by_from[j_inner];
                        }
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
                    }
                }
                else {
                    // generate all pairs to double the prefix length for this node
                    for (size_t j_inner = j; j_inner < j_end; ++j_inner) {
                        nodes.emplace_back(graph.from(node_id), graph.to(order_by_from[j_inner]),
                                           graph.rank(node_id), graph.rank(order_by_from[j_inner]));
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
    
    if (debug_path_graph) {
        cerr << "max rank is " << max_rank << " in graph of size " << nodes.size() << '\n';
    }
    
    // do all nodes have unique ranks?
    return max_rank + 1 == nodes.size();
}

void PathGraph::order_by_rank() {
    
    // re-order the nodes (nodes shouldn't exist yet)
    for (uint64_t i = 0; i < nodes.size(); i++) {
        while (rank(i) != i) {
            std::swap(nodes[rank(i)], nodes[i]);
        }
    }
    
    if (debug_path_graph) {
        cerr << "reordered graph by rank:\n";
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
