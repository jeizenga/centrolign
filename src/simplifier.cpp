#include "centrolign/simplifier.hpp"

#include <limits>
#include <unordered_set>
#include <queue>
#include <set>
#include <functional>

#include "centrolign/superbubbles.hpp"
#include "centrolign/structure_distances.hpp"
#include "centrolign/step_index.hpp"
#include "centrolign/count_walks.hpp"
#include "centrolign/trie.hpp"
#include "centrolign/utility.hpp"
#include "centrolign/compacted_graph.hpp"
#include "centrolign/logging.hpp"

namespace centrolign {

using namespace std;

const bool Simplifier::debug = false;

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
        
        const auto& chain = bub_tree.structures_inside(chain_id);
        
        // count walks in the child snarls
        std::vector<uint64_t> walk_sub_counts(chain.size());
        
        std::vector<bool> do_split(chain.size(), false);
        window_prod_t prod;
        
        size_t window_width = 0;
        size_t window_begin = 0;
        for (size_t i = 0; i < chain.size(); ++i) {
            
            if (debug) {
                auto bs = bub_tree.structure_boundaries(chain[i]);
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
            std::tie(min_bub_dist, max_bub_dist) = bub_dists.structure_min_max_dist(chain[i]);
            
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
                window_width -= bub_dists.structure_min_max_dist(chain[window_begin]).first;
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
        for (size_t i = 0; i < do_split.size(); ) {
            if (do_split[i]) {
                size_t j = i + 1;
                while (j < do_split.size() && do_split[j]) {
                    ++j;
                }
                
                simplify_chain_interval(graph, step_index, bub_tree,
                                        interval_rev_tries, node_to_trie,
                                        chain_id, i, j);
                
                // update the walk count
                simp_count_walks = sat_mult(simp_count_walks, count_walks(interval_rev_tries.back().first));
                
                i = j;
            }
            else {
                // we can already add walks through this bubble to the tally
                simp_count_walks = sat_mult(simp_count_walks, walk_sub_counts[i]);
                ++i;
            }
        }
        
        chain_subwalks[chain_id] = simp_count_walks;
    }
    
    return std::move(perform_simplification(graph, tableau, step_index, interval_rev_tries, node_to_trie));
}

void Simplifier::simplify_chain_interval(const BaseGraph& graph, const StepIndex& step_index,
                                         const SuperbubbleTree& bub_tree,
                                         std::vector<std::pair<Trie, uint64_t>>& interval_rev_tries,
                                         std::vector<size_t>& node_to_trie,
                                         uint64_t chain_id, size_t begin, size_t end) const {
    
    
    const auto& chain = bub_tree.structures_inside(chain_id);
    
    uint64_t start_id = bub_tree.structure_boundaries(chain[begin]).first;
    uint64_t end_id = bub_tree.structure_boundaries(chain[end - 1]).second;
    
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
    
    // iterate in topological order to guarantee that edges to predecessors can be
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
    
    logging::log(logging::Debug, "Node count increased from " + to_string(graph.node_size()) + " to " + to_string(simplified.graph.node_size()) + " in simplified graph");
    
    return simplified;
}

ExpandedGraph Simplifier::targeted_simplify(const BaseGraph& graph, const SentinelTableau& tableau,
                                            const vector<uint64_t>& node_ids, size_t distance) const {
        
    if (debug) {
        std::cerr << "beginning targeted simplification algorithm at distance " << distance << " from nodes:\n";
        for (auto n : node_ids) {
            cerr << '\t' << n << '\n';
        }
    }
    
    logging::log(logging::Debug, "Targeted simplify using " + to_string(node_ids.size()) + " of " + to_string(graph.node_size()) + " with walk distance " + to_string(distance));
    
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
            while (true) {
                if (node_id_set.count(here)) {
                    if (pos + 1 + distance <= compacted.label_size(comp_id)) {
                        traversed_intervals.emplace_back(comp_id, pos, pos + 1 + distance);
                    }
                    else {
                        traversed_intervals.emplace_back(comp_id, pos, compacted.label_size(comp_id));
                        for (auto next_id : compacted.next(comp_id)) {
                            pqueue.emplace(compacted.label_size(comp_id) - pos - 1, next_id);
                        }
                    }
                }
                
                if (here == compacted.back(comp_id)) {
                    break;
                }
                
                ++pos;
                here = graph.next(here).front();
            }
        }
    }
    
    if (debug) {
        cerr << "traversed intervals before dijkstra:\n";
        for (auto interval : traversed_intervals) {
            std::cerr << get<0>(interval) << " (" << compacted.front(get<0>(interval)) << " to " << compacted.back(get<0>(interval)) << ") " << get<1>(interval) << ":" << get<2>(interval) << "\n";
        }
    }
    
    // dijkstra
    // distances correspond to the path traversed before the beginning of the node
    while (!pqueue.empty()) {
        auto top = pqueue.top();
        pqueue.pop();
        if (debug) {
            std::cerr << "dequeue " << top.second << " (" << compacted.front(top.second) << " to " << compacted.back(top.second) << ") at distance " << top.first << "\n";
        }
        if (dequeued.count(top.second)) {
            if (debug) {
                cerr << "filtered out, already traversed\n";
            }
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
    
    if (debug) {
        cerr << "deduplicating traversed intervals\n";
    }
    
    // merge together any of these that ended up on the same compacted node
    sort(traversed_intervals.begin(), traversed_intervals.end());
    size_t removed = 0;
    for (size_t i = 0; i < traversed_intervals.size(); ) {
        size_t j = i + 1;
        while (j < traversed_intervals.size() &&
               get<0>(traversed_intervals[j]) == get<0>(traversed_intervals[i]) &&
               get<1>(traversed_intervals[j]) <= get<2>(traversed_intervals[j - 1])) {
            ++j;
        }
        get<2>(traversed_intervals[i]) = get<2>(traversed_intervals[j - 1]);
        traversed_intervals[i - removed] = traversed_intervals[i];
        removed += (j - i - 1);
        i = j;
    }
    traversed_intervals.resize(traversed_intervals.size() - removed);
    
    if (debug) {
        cerr << "traversed intervals:\n";
        for (auto interval : traversed_intervals) {
            std::cerr << get<0>(interval) << " (" << compacted.front(get<0>(interval)) << " to " << compacted.back(get<0>(interval)) << ") " << get<1>(interval) << ":" << get<2>(interval) << "\n";
        }
    }
    
    // convert traversals to the original node IDs
    vector<uint64_t> simplify_node_ids;
    for (const auto& interval : traversed_intervals) {
        auto here = compacted.front(get<0>(interval));
        if (get<1>(interval) == 0 && here != tableau.src_id && here != tableau.snk_id) {
            simplify_node_ids.push_back(here);
        }
        size_t pos = 0;
        while (++pos < get<2>(interval)) {
            here = graph.next(here).front();
            if (pos >= get<1>(interval) && here != tableau.src_id && here != tableau.snk_id) {
                simplify_node_ids.push_back(here);
            }
        }
    }
    
    SuperbubbleTree superbubbles(graph, tableau);
    
    if (debug) {
        cerr << "converted back into node IDs:\n";
        for (auto node_id : simplify_node_ids) {
            cerr << '\t' << node_id << '\n';
        }
    }
    
    // identify the nearest containing superbubble
    std::vector<bool> simplify_superbubble(superbubbles.structure_size(), false);
    std::vector<bool> traversed(graph.node_size(), false);
    for (auto node_id : simplify_node_ids) {
        
        if (debug) {
            cerr << "beginning bubble-finding traversal from " << node_id << '\n';
        }
        
        // TODO: should I move this into the bubble tree?
        
        // handle this as a special case so that afterwards we can focus on just endings
        if (superbubbles.structure_beginning_at(node_id) != -1) {
            if (debug) {
                cerr << "hit a bubble boundary on first node\n";
            }
            simplify_superbubble[superbubbles.structure_beginning_at(node_id)] = true;
            continue;
        }
        
        std::vector<uint64_t> stack;
        if (!traversed[node_id]) {
            stack.push_back(node_id);
        }
        while (!stack.empty()) {
            
            auto id_here = stack.back();
            stack.pop_back();
            
            if (traversed[id_here]) {
                continue;
            }
            
            if (debug) {
                cerr << "traverse to " << id_here << "\n";
            }
            
            if (superbubbles.structure_ending_at(id_here) != -1) {
                // reached the end of the superbubble, and it's not a chain we skipped
                if (debug) {
                    cerr << "identify superbubble\n";
                }
                simplify_superbubble[superbubbles.structure_ending_at(id_here)] = true;
                break;
            }
            
            traversed[id_here] = true;
            
            for (auto next_id : graph.next(id_here)) {
                if (superbubbles.structure_beginning_at(next_id) != -1 &&
                    superbubbles.structure_ending_at(next_id) == -1) {
                    // skip to the end of this chain, and don't mark anything traversed
                    auto bub_id = superbubbles.structure_beginning_at(next_id);
                    auto chain_id = superbubbles.chain_containing(bub_id);
                    auto final_bub_id = superbubbles.structures_inside(chain_id).back();
                    stack.push_back(superbubbles.structure_boundaries(final_bub_id).second);
                }
                else {
                    stack.push_back(next_id);
                }
            }
        }
    }
    
    if (debug) {
        cerr << "superbubbles being simplified:\n";
        for (uint64_t bub_id = 0; bub_id < superbubbles.structure_size(); ++bub_id) {
            auto b = superbubbles.structure_boundaries(bub_id);
            cerr << b.first << "," << b.second << " ? " << simplify_superbubble[bub_id] << '\n';
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
        
        const auto& chain = superbubbles.structures_inside(chain_id);
        
        for (size_t i = 0; i < chain.size(); ) {
            if (simplify_superbubble[chain[i]]) {
                size_t j = i + 1;
                while (j < chain.size() && simplify_superbubble[chain[j]]) {
                    ++j;
                }
                
                simplify_chain_interval(graph, step_index, superbubbles, interval_rev_tries,
                                        node_to_trie, chain_id, i, j);
                
                i = j;
            }
            else {
                ++i;
            }
        }
    }
    
    return std::move(perform_simplification(graph, tableau, step_index, interval_rev_tries, node_to_trie));
}

vector<vector<uint64_t>> Simplifier::mergeable_nodes(const Trie& trie) const {
    
    vector<vector<uint64_t>> mergeable_sets;
    
    std::vector<std::vector<uint64_t>> stack(1);
    // init recursion on the leaves of the trie
    for (uint64_t node_id = 0; node_id < trie.node_size(); ++node_id) {
        if (trie.next_size(node_id) == 0) {
            stack.front().push_back(node_id);
        }
    }
    
    while (!stack.empty()) {
        auto node_set = std::move(stack.back());
        stack.pop_back();
        // group up the nodes by their origin node
        unordered_map<uint64_t, vector<uint64_t>> sets;
        for (auto node_id : node_set) {
            sets[trie.label(node_id)].push_back(node_id);
        }
        for (pair<const uint64_t, vector<uint64_t>>& subset : sets) {
            if (subset.second.size() > 1) {
                // this is a non-trivial set of nodes with the same origin
                
                // get their parents
                vector<uint64_t> parents;
                parents.reserve(subset.second.size());
                for (auto node_id : subset.second) {
                    auto parent_id = trie.get_parent(node_id);
                    if (trie.next_size(parent_id) == 1) {
                        parents.push_back(parent_id);
                    }
                }
                // record this mergeable set
                mergeable_sets.emplace_back(std::move(subset.second));
                if (parents.size() > 1) {
                    // recurse into the parents of these nodes
                    stack.push_back(std::move(parents));
                }
            }
        }
    }
    
    return mergeable_sets;
}

vector<vector<uint64_t>> Simplifier::identify_target_nodes(const vector<vector<uint64_t>>& node_counts) const {
    
    // make structures to translate from joined to unjoined node space
    vector<size_t> ranges(node_counts.size() + 1, 0);
    for (size_t i = 0; i < node_counts.size(); ++i) {
        ranges[i + 1] = ranges[i] + node_counts[i].size();
    }
    vector<uint16_t> idx_to_comp(ranges.back(), 0);
    for (uint16_t i = 1; i < node_counts.size(); ++i) {
        for (size_t j = ranges[i], n = ranges[i + 1]; j < n; ++j) {
            idx_to_comp[j] = i;
        }
    }
    
    // find the cutoff based on the proportion we want to simplify
    auto indexes = range_vector(ranges.back());
    auto cutoff_it = indexes.begin() + min_resimplify_fraction * indexes.size();
    nth_element(indexes.begin(), cutoff_it, indexes.end(),
                [&](size_t i, size_t j) {
        auto ci = idx_to_comp[i];
        auto cj = idx_to_comp[j];
        return node_counts[ci][i - ranges[ci]] < node_counts[cj][j - ranges[cj]];
    });
    
    // choose the final cutoff by also considering the max count
    uint64_t cutoff = min(node_counts[idx_to_comp[*cutoff_it]][*cutoff_it - ranges[idx_to_comp[*cutoff_it]]],
                          max_resimplify_count);
    
    // identify the targets
    vector<vector<uint64_t>> targets(node_counts.size());
    for (size_t i = 0; i < node_counts.size(); ++i) {
        auto& comp_counts = node_counts[i];
        auto& comp_targets = targets[i];
        for (size_t j = 0; j < comp_counts.size(); ++j) {
            if (comp_counts[j] > cutoff) {
                comp_targets.push_back(j);
            }
        }
    }
    
    return targets;
}


}
