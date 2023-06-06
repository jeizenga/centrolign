#include "centrolign/simplifier.hpp"

#include <limits>
#include <unordered_set>
#include <queue>
#include <set>

#include "centrolign/superbubbles.hpp"
#include "centrolign/superbubble_distances.hpp"
#include "centrolign/step_index.hpp"
#include "centrolign/count_walks.hpp"
#include "centrolign/trie.hpp"
#include "centrolign/utility.hpp"
#include "centrolign/compacted_graph.hpp"

namespace centrolign {

using namespace std;

ExpandedGraph Simplifier::simplify(const BaseGraph& graph, const SentinelTableau& tableau) const {
        
    if (debug) {
        std::cerr << "beginning simplification algorithm\n";
    }
    
    SuperbubbleTree bub_tree(graph, tableau);
    SuperbubbleDistances bub_dists(bub_tree, graph);
    StepIndex step_index(graph);
    
    // the number of walks through a chain in the simplified graph
    std::vector<uint64_t> chain_subwalks(bub_tree.chain_size(), 0);
    
    // a bank of tries that we will replace complicated segments with, along with the
    // node ID for the entry point of the snarl that they include
    // note: they are reverse so that we maintain reverse determinism
    std::vector<std::pair<Trie, uint64_t>> interval_rev_tries;
    // for nodes that have been replaced, the trie that they were replaced with
    std::vector<size_t> node_to_trie(graph.node_size(), -1);
    
    for (const auto& feature : bub_tree.postorder()) {
        if (!feature.second) {
            // we only iterate over chains in this postorder
            continue;
        }
        
        uint64_t chain_id = feature.first;
        
        if (debug) {
            std::cerr << "simplifying chain " << chain_id << '\n';
        }
        
        const auto& chain = bub_tree.superbubbles_inside(chain_id);
        
        // count walks in the child snarls
        std::vector<uint64_t> walk_sub_counts(chain.size());
        
        std::vector<bool> do_split(chain.size(), false);
        window_prod_t prod;
        
        size_t window_width = 0;
        size_t window_begin = 0;
        for (size_t i = 0; i < chain.size(); ++i) {
            
            if (debug) {
                auto bs = bub_tree.superbubble_boundaries(chain[i]);
                cerr << "traversing superbubble " << chain[i] << " with boundaries " << bs.first << " " << bs.second << '\n';
            }
            
            NetGraph netgraph(graph, bub_tree, chain[i]);
            walk_sub_counts[i] = count_walks_hierarchical(netgraph,
                                                          [&](uint64_t net_id) -> uint64_t {
                uint64_t feature_id;
                bool is_chain;
                std::tie(feature_id, is_chain) = netgraph.label(net_id);
                return is_chain ? chain_subwalks[feature_id] : 1;
            });
            
            prod *= walk_sub_counts[i];
            
            size_t min_bub_dist, max_bub_dist;
            std::tie(min_bub_dist, max_bub_dist) = bub_dists.superbubble_min_max_dist(chain[i]);
            
            if (max_bub_dist >= preserve_bubble_size) {
                //  we hit a bubble with an allele we want to preserve
                if (debug) {
                    cerr << "skipping bubble of size " << max_bub_dist << '\n';
                }
                window_begin = i + 1;
                window_width = 0;
                prod = 1;
                continue;
            }
            
            // extend the window rightward to include this index
            window_width += min_bub_dist;
            if (window_begin != i) {
                // account for the overlap by 1 base
                --window_width;
            }
            if (debug) {
                std::cerr << "extend window to interval " << window_begin << ":" << i << ", width " << window_width << ", walk count " << prod.str() << '\n';
            }
            // shrink the window from the left to be under the window size
            while (window_width > min_dist_window) {
                window_width -= bub_dists.superbubble_min_max_dist(chain[window_begin]).first;
                if (window_begin != i) {
                    // account for the overlap by 1 base
                    ++window_width;
                }
                prod /= walk_sub_counts[window_begin];
                ++window_begin;
                
                if (debug) {
                    std::cerr << "retract window to interval " << window_begin << ":" << i << ", width " << window_width << ", walk count " << prod.str() << '\n';
                }
            }
            
            if (prod > max_walks) {
                // mark yet-unmarked superbubbles in the window to be split
                for (int64_t j = i; j >= (int64_t) window_begin && !do_split[j]; --j) {
                    if (debug) {
                        std::cerr << "marking " << j << " for splitting\n";
                    }
                    do_split[j] = true;
                }
            }
        }
        
        // the number of walks in the simplified chain
        uint64_t simp_count_walks = 1;
        
        // identify the intervals that we need to split
        std::vector<std::pair<size_t, size_t>> split_intervals;
        for (size_t i = 0; i < do_split.size(); ) {
            if (do_split[i]) {
                size_t j = i + 1;
                while (j < do_split.size() && do_split[j]) {
                    ++j;
                }
                split_intervals.emplace_back(i, j);
                i = j;
            }
            else {
                // we can already add walks through this bubble to the tally
                simp_count_walks = sat_mult(simp_count_walks, walk_sub_counts[i]);
                ++i;
            }
        }
        
        // split the graph into a trie at each split interval
        for (const auto& interval : split_intervals) {
            
            simplify_chain_interval(graph, step_index, bub_tree,
                                    interval_rev_tries, node_to_trie,
                                    chain_id, interval.first, interval.second);
            
            // update the walk count
            simp_count_walks = sat_mult(simp_count_walks, count_walks(interval_rev_tries.back().first));
        }
        
        chain_subwalks[chain_id] = simp_count_walks;
    }
    
    return perform_simplification(graph, tableau, step_index, interval_rev_tries, node_to_trie);
}

void Simplifier::simplify_chain_interval(const BaseGraph& graph, const StepIndex& step_index,
                                         const SuperbubbleTree& bub_tree,
                                         std::vector<std::pair<Trie, uint64_t>>& interval_rev_tries,
                                         std::vector<size_t>& node_to_trie,
                                         uint64_t chain_id, size_t begin, size_t end) const {
    
    
    const auto& chain = bub_tree.superbubbles_inside(chain_id);
    
    uint64_t start_id = bub_tree.superbubble_boundaries(chain[begin]).first;
    uint64_t end_id = bub_tree.superbubble_boundaries(chain[end - 1]).second;
    
    if (debug) {
        std::cerr << "splitting interval of bubbles from node " << start_id << " to " << end_id << '\n';
    }
    
    size_t trie_idx = interval_rev_tries.size();
    interval_rev_tries.emplace_back();
    interval_rev_tries.back().second = start_id;
    auto& interval_rev_trie = interval_rev_tries.back().first;
    
    // iterate over paths on the start node
    // note: all paths must traverse the entire chain interval or else they would be connected
    // to the sentinel, which would disrupt a superbubble
    for (const auto& occurrence : step_index.path_steps(end_id)) {
        
        if (debug) {
            std::cerr << "walking across bubble on path " << occurrence.first << " starting from step " << occurrence.second << '\n';
        }
        
        // gather the path interval sequence and record the node->trie correspondence
        // note: we leave the exit node alone
        const auto& path = graph.path(occurrence.first);
        std::vector<uint64_t> rev_path_seq;
        for (size_t i = occurrence.second; path[i] != start_id; --i) {
            auto node_id = path[i];
            
            node_to_trie[node_id] = trie_idx;
            
            rev_path_seq.push_back(node_id);
        }
        
        if (debug) {
            std::cerr << "inserting into trie reverse sequence ";
            for (size_t i = 0; i < rev_path_seq.size(); ++i) {
                if (i) {
                    std::cerr << ',';
                }
                std::cerr << rev_path_seq[i];
            }
            std::cerr << " from path " << occurrence.first << '\n';
        }
        
        // add to the trie
        interval_rev_trie.insert_sequence(graph.path_name(occurrence.first), rev_path_seq);
    }
}

ExpandedGraph Simplifier::perform_simplification(const BaseGraph& graph,
                                                 const SentinelTableau& tableau,
                                                 const StepIndex& step_index,
                                                 std::vector<std::pair<Trie, uint64_t>>& interval_rev_tries,
                                                 const std::vector<size_t>& node_to_trie) const {
        
    // construct the simplified graph
    ExpandedGraph simplified;
    
    // init empty paths
    for (uint64_t path_id = 0; path_id < graph.path_size(); ++path_id) {
        simplified.graph.add_path(graph.path_name(path_id));
    }
    
    // the forward translation for nodes outside of tries and for the sink of a trie
    std::vector<uint64_t> forward_translation(graph.node_size(), -1);
    
    // we iterate in topological order to guarantee that edges to predecessors can be
    // contructed in the current iteration and paths can be incrementally extended
    for (uint64_t node_id : topological_order(graph)) {
        if (node_to_trie[node_id] == -1) {
            // we haven't simplified this node
            
            // copy it over
            auto new_node_id = simplified.graph.add_node(graph.label(node_id));
            simplified.back_translation.push_back(node_id);
            forward_translation[node_id] = new_node_id;
            
            if (debug) {
                std::cerr << "directly copying node " << node_id << " as " << new_node_id << '\n';
            }
            
            for (auto prev_id : graph.previous(node_id)) {
                simplified.graph.add_edge(forward_translation[prev_id], new_node_id);
            }
            
            // add the new copies onto its paths
            for (const auto& occurrence : step_index.path_steps(node_id)) {
                simplified.graph.extend_path(occurrence.first, new_node_id);
            }
        }
        else if (interval_rev_tries[node_to_trie[node_id]].first.node_size() != 0) {
            // this is part of a trie that we haven't inserted yet
            
            size_t trie_idx = node_to_trie[node_id];
            auto& trie = interval_rev_tries[trie_idx].first;
            uint64_t entry_id = interval_rev_tries[trie_idx].second;
            std::vector<uint64_t> trie_forward_translation(trie.node_size(), -1);
                        
            if (debug) {
                std::cerr << "encountered trie " << trie_idx << " at " << node_id << '\n';
            }
            
            // this should be true for tries constructed from a chain
            assert(trie.next_size(trie.get_root()) == 1);
            auto trie_sink_id = trie.next(trie.get_root()).front();
            
            // make nodes for the mergeable leaves
            for (auto& mergeable_node_group : mergeable_nodes(trie)) {
                auto original_node_id = trie.label(mergeable_node_group.front());
                auto new_node_id = simplified.graph.add_node(graph.label(original_node_id));
                for (auto trie_id : mergeable_node_group) {
                    trie_forward_translation[trie_id] = new_node_id;
                }
                simplified.back_translation.push_back(original_node_id);
                if (debug) {
                    cerr << "make merged node " << new_node_id << " from set:\n";
                    for (auto n : mergeable_node_group) {
                        cerr << ' ' << n;
                    }
                    cerr << '\n';
                }
            }
                        
            // copy the other nodes nodes over
            for (uint64_t trie_id = 0; trie_id < trie.node_size(); ++trie_id) {
                if (trie_id == trie.get_root() || trie_forward_translation[trie_id] != -1) {
                    // the root is a placeholder
                    continue;
                }
                
                auto original_node_id = trie.label(trie_id);
                auto new_node_id = simplified.graph.add_node(graph.label(original_node_id));
                trie_forward_translation[trie_id] = new_node_id;
                simplified.back_translation.push_back(original_node_id);
                
                if (debug) {
                    cerr << "make new node " << new_node_id << " from trie node " << trie_id << " representing " << original_node_id << '\n';
                }
            }
            
            // identify the edges internal to the trie
            std::unordered_set<std::pair<uint64_t, uint64_t>> edges;
            for (uint64_t trie_id = 0; trie_id < trie.node_size(); ++trie_id) {
                if (trie_id == trie.get_root() || trie_id == trie_sink_id) {
                    // parent edge is nonexistent or irrelevant
                    continue;
                }
                
                // parent edges are internal
                edges.emplace(trie_forward_translation[trie_id],
                              trie_forward_translation[trie.get_parent(trie_id)]);
            }
            
            // copy the paths and identify edges out of the trie
            for (uint64_t trie_path_id = 0; trie_path_id < trie.path_size(); ++trie_path_id) {
                auto path_id = graph.path_id(trie.path_name(trie_path_id));
                const auto& trie_path = trie.path(trie_path_id);
                for (int64_t i = trie_path.size() - 1; i >= 0; --i) {
                    simplified.graph.extend_path(path_id, trie_forward_translation[trie_path[i]]);
                }
                // path ends connect to the chain start ID outside the trie
                edges.emplace(forward_translation[entry_id],
                              trie_forward_translation[trie_path.back()]);
            }
            // actually add the edges
            for (const auto& edge : edges) {
                simplified.graph.add_edge(edge.first, edge.second);
            }
            
            // this node alone goes into the main forward translation
            forward_translation[trie.label(trie_sink_id)] = trie_forward_translation[trie_sink_id];
            
            // clean out the trie to mark it as having been added
            trie.clear();
        }
    }
    // the source and sink are not part of bubbles, so they should not have been simplified
    simplified.tableau.src_id = forward_translation[tableau.src_id];
    simplified.tableau.src_sentinel = tableau.src_sentinel;
    simplified.tableau.snk_id = forward_translation[tableau.snk_id];
    simplified.tableau.snk_sentinel = tableau.snk_sentinel;
    
    return simplified;
}

ExpandedGraph Simplifier::targeted_simplify(const BaseGraph& graph, const SentinelTableau& tableau,
                                            const vector<uint64_t>& node_ids, size_t distance) const {
    
    // make a unipath graph
    CompactedGraph compacted(graph);
    
    // the queue we'll use for dijkstra
    priority_queue<pair<size_t, uint64_t>,
                   vector<pair<size_t, uint64_t>>,
                   greater<pair<size_t, uint64_t>>> pqueue;
    unordered_set<uint64_t> dequeued;
    
    vector<tuple<uint64_t, size_t, size_t>> traversed_intervals;
    
    // find where the source nodes fit into the compacted nodes
    {
        unordered_set<uint64_t> node_id_set;
        for (size_t i = 0; i < node_ids.size(); ++i) {
            node_id_set.insert(node_ids[i]);
        }
        for (uint64_t comp_id = 0; comp_id < compacted.node_size(); ++comp_id) {
            uint64_t here = compacted.front(comp_id);
            size_t pos = 0;
            if (node_id_set.count(here)) {
                if (pos + distance <= compacted.label_size(here)) {
                    traversed_intervals.emplace_back(here, pos, pos + distance);
                }
                else {
                    pqueue.emplace(compacted.label_size(here) - pos, here);
                }
            }
            while (here != compacted.back(comp_id)) {
                ++pos;
                here = graph.next(here).front();
                if (pos + distance <= compacted.label_size(here)) {
                    traversed_intervals.emplace_back(here, pos, pos + distance);
                }
                else {
                    pqueue.emplace(compacted.label_size(here) - pos, here);
                }
            }
        }
    }
    
    // dijkstra
    while (!pqueue.empty()) {
        auto top = pqueue.top();
        pqueue.pop();
        if (dequeued.count(top.second)) {
            continue;
        }
        dequeued.insert(top.second);
        
        auto dist_thru = top.first + compacted.label_size(top.second);
        if (dist_thru < distance) {
            // we can reach the next node
            traversed_intervals.emplace_back(top.second, 0, compacted.label_size(top.second));
            for (auto next_id : compacted.next(top.second)) {
                pqueue.emplace(dist_thru, next_id);
            }
        }
        else {
            // we cannot reach the next node
            traversed_intervals.emplace_back(top.second, 0, distance - top.first);
        }
    }
    
    // merge together any of these that ended up on the same compacted node
    sort(traversed_intervals.begin(), traversed_intervals.end());
    size_t removed = 0;
    for (size_t i = 0; i < traversed_intervals.size(); ) {
        size_t j = i + 1;
        while (j < traversed_intervals.size() && get<0>(traversed_intervals[j]) == get<0>(traversed_intervals[i])) {
            ++j;
        }
        get<2>(traversed_intervals[i]) = get<2>(traversed_intervals[j - 1]);
        traversed_intervals[i - removed] = traversed_intervals[i];
        removed += (j - i - 1);
        i = j;
    }
    traversed_intervals.resize(traversed_intervals.size() - removed);
    
    // convert traversals to the original node IDs
    vector<uint64_t> simplify_node_ids;
    for (const auto& interval : traversed_intervals) {
        auto here = compacted.front(get<0>(interval));
        if (get<1>(interval) == 0) {
            simplify_node_ids.push_back(here);
        }
        size_t pos = 0;
        while (++pos < get<2>(interval)) {
            here = graph.next(here).front();
            if (pos >= get<0>(interval)) {
                simplify_node_ids.push_back(here);
            }
        }
    }
    
    SuperbubbleTree superbubbles(graph, tableau);
    
    // identify the nearest containing superbubble
    std::vector<bool> simplify_superbubble(false, superbubbles.superbubble_size());
    std::vector<bool> traversed(false, graph.node_size());
    for (auto node_id : simplify_node_ids) {
        
        int rel_depth = 0;
        while (graph.next_size(node_id) != 0) {
            node_id = graph.next(node_id).front();
            if (traversed[node_id]) {
                // we're already found a superbubble by walking this node
                break;
            }
            traversed[node_id] = true;
            if (superbubbles.superbubble_ending_at(node_id) != -1) {
                if (rel_depth == 0) {
                    // reached the end of the superbubble
                    simplify_superbubble[superbubbles.superbubble_ending_at(node_id)] = true;
                    break;
                }
                else {
                    // exiting a contained superbubble
                    --rel_depth;
                }
            }
            if (superbubbles.superbubble_beginning_at(node_id) != -1) {
                // entering a contained superbubble
                ++rel_depth;
            }
        }
    }
    
    // a bank of tries that we will replace complicated segments with, along with the
    // node ID for the entry point of the snarl that they include
    // note: they are reverse so that we maintain reverse determinism
    std::vector<std::pair<Trie, uint64_t>> interval_rev_tries;
    // for nodes that have been replaced, the trie that they were replaced with
    std::vector<size_t> node_to_trie(graph.node_size(), -1);
    
    StepIndex step_index(graph);
    
    for (const auto& feature : superbubbles.postorder()) {
        if (!feature.second) {
            // we only iterate over chains in this postorder
            continue;
        }
        
        uint64_t chain_id = feature.first;
        
        if (debug) {
            std::cerr << "simplifying chain " << chain_id << '\n';
        }
        
        const auto& chain = superbubbles.superbubbles_inside(chain_id);
        
        for (size_t i = 0; i < chain.size(); ) {
            if (simplify_superbubble[chain[i]]) {
                size_t j = i;
                while (j < chain.size() && simplify_superbubble[chain[j]]) {
                    ++j;
                }
                
                simplify_chain_interval(graph, step_index, superbubbles,
                                        interval_rev_tries, node_to_trie,
                                        chain_id, i, j);
                
                i = j;
            }
            else {
                ++i;
            }
        }
    }
    
    return perform_simplification(graph, tableau, step_index, interval_rev_tries, node_to_trie);
}

vector<vector<uint64_t>> Simplifier::mergeable_nodes(const Trie& trie) const {
    
    vector<vector<uint64_t>> mergeable_sets;
    // the core recursive algorithm to find nodes that can be merged
    function<void(const vector<uint64_t>&)> find_mergeable = [&](const vector<uint64_t>& node_set) {
        
        unordered_map<uint64_t, vector<uint64_t>> sets;
        for (auto node_id : node_set) {
            sets[trie.label(node_id)].push_back(node_id);
        }
        for (pair<const uint64_t, vector<uint64_t>>& subset : sets) {
            if (subset.second.size() > 1) {
                // this is a set of leaves (or leaf parents that are all the same node)
                
                vector<uint64_t> parents;
                parents.reserve(subset.second.size());
                for (auto node_id : subset.second) {
                    auto parent_id = trie.get_parent(node_id);
                    if (trie.next_size(parent_id) == 1) {
                        parents.push_back(parent_id);
                    }
                }
                // record this mergeable set
                mergeable_sets.emplace_back(move(subset.second));
                if (parents.size() > 1) {
                    // recurse into the parents of these nodes
                    find_mergeable(parents);
                }
            }
        }
    };
    
    // init recursion on the leaves of the trie
    std::vector<uint64_t> leaves;
    for (uint64_t node_id = 0; node_id < trie.node_size(); ++node_id) {
        if (trie.next_size(node_id) == 0) {
            leaves.push_back(node_id);
        }
    }
    find_mergeable(leaves);
    
    return mergeable_sets;
}


}
