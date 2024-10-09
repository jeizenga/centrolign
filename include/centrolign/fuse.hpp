#ifndef centrolign_fuse_hpp
#define centrolign_fuse_hpp

#include <cstdint>
#include <unordered_set>
#include <map>
#include <iostream>

#include "centrolign/graph.hpp"
#include "centrolign/modify_graph.hpp"
#include "centrolign/alignment.hpp"
#include "centrolign/union_find.hpp"

namespace centrolign {

// merge two graphs along an alignment (occurs in place for the first
// of the graphs), the in-place graph should be node_id1 in the alignment
// merge two graphs along an alignment (occurs in place for the first
// of the graphs), the in-place graph should be node_id1 in the alignment
template<class Graph1, class Graph2>
void fuse(Graph1& dest, const Graph2& source,
          const SentinelTableau& dest_table,
          const SentinelTableau& source_table,
          const Alignment& alignment);


// create a new graph that fuses together portions of the original graph
// that are aligned together according to a set of alignments
template<class BGraph>
BaseGraph internal_fuse(const BGraph& graph,
                        const std::vector<Alignment>& alignments,
                        const SentinelTableau* tableau_in = nullptr,
                        SentinelTableau* tableau_out = nullptr,
                        const Alignment* alignment_in = nullptr,
                        Alignment* alignment_out = nullptr);





/*****************************
 * Template implementations
 *****************************/


template<class Graph1, class Graph2>
void fuse(Graph1& dest, const Graph2& source,
          const SentinelTableau& dest_table,
          const SentinelTableau& source_table,
          const Alignment& alignment) {
    
    static const bool debug_fuse = false;
    
    if (debug_fuse) {
        std::cerr << "beginning fuse algorithm\n";
    }
    
    std::vector<uint64_t> trans(source.node_size(), -1);
    
    // record match nodes
    for (const auto& aln_pair : alignment) {
        if (aln_pair.node_id1 != AlignedPair::gap &&
            aln_pair.node_id2 != AlignedPair::gap &&
            dest.label(aln_pair.node_id1) == source.label(aln_pair.node_id2)) {
            trans[aln_pair.node_id2] = aln_pair.node_id1;
        }
    }
    
    // join the sentinels
    trans[source_table.src_id] = dest_table.src_id;
    trans[source_table.snk_id] = dest_table.snk_id;
            
    // add any unmatched nodes
    for (uint64_t node_id2 = 0; node_id2 < source.node_size(); ++node_id2) {
        if (trans[node_id2] == -1) {
            trans[node_id2] = dest.add_node(source.label(node_id2));
        }
    }
    
    if (debug_fuse) {
        std::cerr << "adding substitution edges\n";
    }
    
    // add substitution edges from the destination graph
    for (size_t i = 0; i < alignment.size(); ++i) {
        const auto& aln_pair = alignment[i];
        if (aln_pair.node_id1 != AlignedPair::gap &&
            aln_pair.node_id2 != AlignedPair::gap &&
            dest.label(aln_pair.node_id1) != source.label(aln_pair.node_id2)) {
            
            // look right for the adjacent node
            for (size_t j = i + 1; j < alignment.size(); ++j) {
                if (alignment[j].node_id1 != AlignedPair::gap) {
                    dest.add_edge(trans[alignment[i].node_id2], alignment[j].node_id1);
                    break;
                }
            }
            
            // look left for the adjacent node
            for (int64_t j = i - 1; j >= 0; --j) {
                if (alignment[j].node_id1 != AlignedPair::gap) {
                    dest.add_edge(alignment[j].node_id1, trans[alignment[i].node_id2]);
                    break;
                }
            }
        }
    }
    
    if (debug_fuse) {
        std::cerr << "adding in uncovered edges\n";
    }
    
    // add in any missing edges
    for (uint64_t node_id2 = 0; node_id2 < source.node_size(); ++node_id2) {
        
        auto new_id = trans[node_id2];
        std::unordered_set<uint64_t> next_ids;
        for (auto next_id : dest.next(new_id)) {
            next_ids.insert(next_id);
        }
        
        for (auto next_id : source.next(node_id2)) {
            auto new_next_id = trans[next_id];
            if (!next_ids.count(new_next_id)) {
                dest.add_edge(new_id, new_next_id);
            }
        }
    }
    
    
    if (debug_fuse) {
        std::cerr << "copying paths\n";
    }
    
    // add the paths
    for (uint64_t new_path_id = dest.path_size(), path_id = 0; path_id < source.path_size(); ++path_id, ++new_path_id) {
        dest.add_path(source.path_name(path_id));
        for (auto node_id : source.path(path_id)) {
            dest.extend_path(new_path_id, trans[node_id]);
        }
    }
}

template<class BGraph>
BaseGraph internal_fuse(const BGraph& graph,
                        const std::vector<Alignment>& alignments,
                        const SentinelTableau* tableau_in,
                        SentinelTableau* tableau_out,
                        const Alignment* alignment_in,
                        Alignment* alignment_out) {

    if (tableau_out && !tableau_in) {
        throw std::runtime_error("Outputting a sentinel tableau during internal fuse requires an input tableau");
    }
    if (alignment_out && !alignment_in) {
        throw std::runtime_error("Outputting an updated alignment during internal fuse requires an input alignment");
    }
        
    // figure out the set of transitive merges implied by the alignments
    UnionFind transitive_merges(graph.node_size());
    for (const auto& alignment : alignments) {
        for (const auto& aln_pair : alignment) {
            if (aln_pair.node_id1 != AlignedPair::gap && aln_pair.node_id2 != AlignedPair::gap) {
                transitive_merges.union_groups(aln_pair.node_id1, aln_pair.node_id2);
            }
        }
    }
    
    // create the nodes for the fused graph
    BaseGraph fused;
    std::vector<uint64_t> trans(graph.node_size(), -1);
    for (const auto& group : transitive_merges.get_groups()) {
        
        // ordered map to remove machine-dependency on hash order
        std::map<char, std::vector<uint64_t>> label_subgroups;
        for (auto node_id : group) {
            label_subgroups[graph.label(node_id)].push_back(node_id);
        }
        
        for (const auto& subgroup : label_subgroups) {
            
            auto new_node_id = fused.add_node(subgroup.first);
            for (auto node_id : subgroup.second) {
                trans[node_id] = new_node_id;
            }
        }
    }
    
    {
        // we no longer need this, can clear the memory
        auto dummy = std::move(transitive_merges);
    }
    
    // update the sentinel info
    if (tableau_out) {
        tableau_out->src_id = trans[tableau_in->src_id];
        tableau_out->snk_id = trans[tableau_in->snk_id];
        tableau_out->src_sentinel = tableau_in->src_sentinel;
        tableau_out->snk_sentinel = tableau_in->snk_sentinel;
    }
    
    // recreate the edges
    // note: recreate the edges lists as unordered set for O(1) membership tests TODO: maybe overkill
    std::vector<std::unordered_set<uint64_t>> translated_edges(fused.node_size());
    for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
        auto fused_id = trans[node_id];
        auto& edges = translated_edges[fused_id];
        for (auto next_id : graph.next(node_id)) {
            auto fused_next_id = trans[next_id];
            if (edges.insert(fused_next_id).second) {
                // this is the first time we've seen this edge
                fused.add_edge(fused_id, fused_next_id);
            }
        }
    }
    
    
    // recreate the paths
    for (uint64_t path_id = 0; path_id < graph.path_size(); ++path_id) {
        auto fused_path_id = fused.add_path(graph.path_name(path_id));
        for (auto node_id : graph.path(path_id)) {
            fused.extend_path(fused_path_id, trans[node_id]);
        }
    }
    
    if (alignment_out) {
        alignment_out->clear();
        alignment_out->reserve(alignment_in->size());
        for (const auto& aln_pair_in : *alignment_in) {
            alignment_out->emplace_back(aln_pair_in);
            auto& aln_pair = alignment_out->back();
            if (aln_pair.node_id1 != AlignedPair::gap) {
                aln_pair.node_id1 = trans[aln_pair.node_id1];
            }
            if (aln_pair.node_id2 != AlignedPair::gap) {
                aln_pair.node_id2 = trans[aln_pair.node_id2];
            }
        }
    }
    
    return fused;
}

    
}

#endif /* centrolign_fuse_hpp */
