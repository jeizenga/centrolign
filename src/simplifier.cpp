#include "centrolign/simplifier.hpp"

#include <limits>

#include "centrolign/superbubbles.hpp"
#include "centrolign/superbubble_distances.hpp"
#include "centrolign/step_index.hpp"
#include "centrolign/count_walks.hpp"
#include "centrolign/trie.hpp"
#include "centrolign/utility.hpp"

namespace centrolign {

using namespace std;

ExpandedGraph Simplifier::simplify(const BaseGraph& graph, const SentinelTableau& tableau) const{
    
    static const bool debug = false;
    
    if (debug) {
        std::cerr << "beginning simplification algorithm\n";
    }
    
    SuperbubbleTree bub_tree(graph, tableau);
    SuperbubbleDistances bub_dists(bub_tree, graph);
    StepIndex step_index(graph);
    
    // the number of walks through a chain in the simplified graph
    std::vector<uint64_t> chain_subwalks(bub_tree.chain_size(), 0);
    
    // a bank of tries that we will replace complicated segments with
    std::vector<Trie> interval_tries;
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
        long_product_t prod;
        
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
            
            uint64_t start_id = bub_tree.superbubble_boundaries(chain[interval.first]).first;
            uint64_t end_id = bub_tree.superbubble_boundaries(chain[interval.second - 1]).second;
            
            if (debug) {
                std::cerr << "splitting interval of bubbles from node " << start_id << " to " << end_id << '\n';
            }
            
            size_t trie_idx = interval_tries.size();
            interval_tries.emplace_back();
            auto& interval_trie = interval_tries.back();
            
            // iterate over paths on the start node
            // note: all paths must traverse the entire chain interval or else they would be connected
            // to the sentinel, which would disrupt a superbubble
            for (const auto& occurrence : step_index.path_steps(start_id)) {
                
                if (debug) {
                    std::cerr << "walking across bubble on path " << occurrence.first << " starting from step " << occurrence.second << '\n';
                }
                
                // gather the path interval sequence and record the node->trie correspondence
                // note: we leave the exit node alone
                const auto& path = graph.path(occurrence.first);
                std::string path_seq;
                for (size_t i = occurrence.second; path[i] != end_id; ++i) {
                    auto node_id = path[i];
                    
                    node_to_trie[node_id] = trie_idx;
                    
                    path_seq.push_back(graph.label(node_id));
                }
                
                if (debug) {
                    std::cerr << "inserting into trie sequence " << path_seq << " from path " << occurrence.first << '\n';
                }
                
                // add to the trie
                interval_trie.insert_sequence(graph.path_name(occurrence.first), path_seq);
                
                // update the walk count
                simp_count_walks = sat_mult(simp_count_walks, count_walks(interval_trie));
            }
            if (debug) {
                std::cerr << "final trie:\n";
                print_graph(interval_trie, std::cerr);
            }
        }
        
        chain_subwalks[chain_id] = simp_count_walks;
    }
    
    // construct the simplified graph
    ExpandedGraph simplified;
    
    // init empty paths
    for (uint64_t path_id = 0; path_id < graph.path_size(); ++path_id) {
        simplified.graph.add_path(graph.path_name(path_id));
    }
    
    // for nodes that weren't entered into tries, their corresponding node in
    // the simplified traph
    std::vector<uint64_t> non_trie_forward_translation(graph.node_size(), -1);
    // when we finish adding a trie, its leaves
    std::vector<std::vector<uint64_t>> trie_leaves(interval_tries.size());
    
    for (uint64_t node_id : topological_order(graph)) {
        if (node_to_trie[node_id] == -1) {
            // we haven't simplified this node
            
            
            // copy it over
            auto new_node_id = simplified.graph.add_node(graph.label(node_id));
            simplified.back_translation.push_back(node_id);
            non_trie_forward_translation[node_id] = new_node_id;
            
            if (debug) {
                std::cerr << "directly copying node " << node_id << " as " << new_node_id << '\n';
            }
            
            if (graph.previous_size(node_id) != 0) {
                auto trie_idx = node_to_trie[graph.previous(node_id).front()];
                if (trie_idx == -1) {
                    // the predecessors are simple copied nodes
                    for (auto prev_id : graph.previous(node_id)) {
                        simplified.graph.add_edge(non_trie_forward_translation[prev_id], new_node_id);
                    }
                }
                else {
                    // the predecessors are leaves of a trie
                    for (auto leaf_id : trie_leaves[trie_idx]) {
                        simplified.graph.add_edge(leaf_id, new_node_id);
                    }
                }
            }
            
            // add the new copies onto its paths
            for (const auto& occurrence : step_index.path_steps(node_id)) {
                simplified.graph.extend_path(occurrence.first, new_node_id);
            }
        }
        else if (trie_leaves[node_to_trie[node_id]].empty()) {
            // this is part of a trie that we haven't inserted yet
            
            size_t trie_idx = node_to_trie[node_id];
            const auto& trie = interval_tries[trie_idx];
            auto& leaves = trie_leaves[trie_idx];
            std::vector<uint64_t> trie_forward_translation(trie.node_size(), -1);
            
            if (debug) {
                std::cerr << "encountered trie " << trie_idx << " at " << node_id << '\n';
            }
            
            // make space for the back translation
            for (size_t i = 0, n = trie.node_size() - 1; i < n; ++i) {
                simplified.back_translation.emplace_back(-1);
            }
            
            // this should be true for tries constructed from a chain
            assert(trie.next_size(trie.get_root()) == 1);
            auto entry_point_id = trie.next(trie.get_root()).front();
            
            for (uint64_t trie_id = 0; trie_id < trie.node_size(); ++trie_id) {
                if (trie_id == trie.get_root()) {
                    // the root is a placeholder
                    continue;
                }
                auto new_node_id = simplified.graph.add_node(trie.label(trie_id));
                trie_forward_translation[trie_id] = new_node_id;
                if (trie.next_size(trie_id) == 0) {
                    // record leaves to connect to later
                    leaves.push_back(new_node_id);
                }
                if (trie_id != entry_point_id) {
                    // we can connect this within the trie
                    simplified.graph.add_edge(trie_forward_translation[trie.previous(trie_id).front()],
                                              new_node_id);
                }
                else {
                    // we have to connect this outside the trie
                    // note: because we iterate in topological order, the node we encounter
                    // the trie on must also be the source node
                    for (auto prev_id : graph.previous(node_id)) {
                        simplified.graph.add_edge(non_trie_forward_translation[prev_id], new_node_id);
                    }
                }
            }
            
            // find the ordinal offset value of the step for each path on the entry point
            std::unordered_map<uint64_t, size_t> offset;
            for (auto& occurrence : step_index.path_steps(node_id)) {
                offset[occurrence.first] = occurrence.second;
            }
            
            // copy the paths, and also use them to get the correspondence back to the
            // original graph for the back translation
            for (uint64_t trie_path_id = 0; trie_path_id < trie.path_size(); ++trie_path_id) {
                auto path_id = graph.path_id(trie.path_name(trie_path_id));
                auto trie_path = trie.path(trie_path_id);
                const auto& path = graph.path(path_id);
                for (size_t i = 0, j = offset.at(path_id); i < trie_path.size(); ++i, ++j) {
                    auto new_id = trie_forward_translation[trie_path[i]];
                    auto old_id = path[j];
                    simplified.graph.extend_path(path_id, new_id);
                    simplified.back_translation[new_id] = old_id;
                }
            }
        }
    }
    // the source and sink are not part of bubbles, so they should not have been simplified
    simplified.tableau.src_id = non_trie_forward_translation[tableau.src_id];
    simplified.tableau.src_sentinel = tableau.src_sentinel;
    simplified.tableau.snk_id = non_trie_forward_translation[tableau.snk_id];
    simplified.tableau.snk_sentinel = tableau.snk_sentinel;
    
    return simplified;
}

}
