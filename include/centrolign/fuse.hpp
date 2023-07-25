#ifndef centrolign_fuse_hpp
#define centrolign_fuse_hpp

#include <cstdint>
#include <unordered_set>
#include <iostream>

#include "centrolign/graph.hpp"
#include "centrolign/modify_graph.hpp"
#include "centrolign/alignment.hpp"

namespace centrolign {

// merge two graphs along an alignment (occurs in place for the first
// of the graphs), the in-place graph should be node_id1 in the alignment
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
    
//    print_graph(dest, std::cerr);
}

}

#endif /* centrolign_fuse_hpp */
