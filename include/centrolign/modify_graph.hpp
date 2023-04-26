#ifndef centrolign_modify_graph_hpp
#define centrolign_modify_graph_hpp

#include <cstdint>
#include <unordered_set>

#include "centrolign/graph.hpp"
#include "centrolign/alignment.hpp"

namespace centrolign {

// make simple graphs with encoded chars and an embedded path, but no sentinels
SequenceGraph make_sequence_graph(const std::string& name,
                                  const std::string& sequence);
BaseGraph make_base_graph(const std::string& name,
                          const std::string& sequence);

// append a graph as a separate connected component onto another, including
// embedded paths
template<class Graph1, class Graph2>
void append_component(Graph1& appending, const Graph2& component);

/*
 * A convenience struct to keep track of sentinel node information
 */
struct SentinelTableau {
    uint64_t src_id = -1;
    uint64_t snk_id = -1;
    char src_sentinel = 0;
    char snk_sentinel = 0;
};

// add sentinel nodes and construct a tableau
SentinelTableau add_sentinels(BaseGraph& graph,
                              char src_sentinel, char snk_sentinel);

// change the sentinel node characters
void reassign_sentinels(BaseGraph& graph, SentinelTableau& tableau,
                        char src_sentinel, char snk_sentinel);

// merge two graphs along an alignment (occurs in place for the first
// of the graphs), the in-place graph should be node_id1 in the alignment
template<class Graph1, class Graph2>
void fuse(Graph1& dest, const Graph2& source,
          const SentinelTableau& dest_table,
          const SentinelTableau& source_table,
          const Alignment& alignment);




/*
 * Template implementations
 */

template<class Graph1, class Graph2>
void append_component(Graph1& appending, const Graph2& component) {
    
    uint64_t prev_node_size = appending.node_size();
    
    // add the nodes
    for (uint64_t node_id = 0; node_id < component.node_size(); ++node_id) {
        appending.add_node(component.label(node_id));
    }
    
    // add the edges
    for (uint64_t node_id = 0; node_id < component.node_size(); ++node_id) {
        for (uint64_t next_id : component.next(node_id)) {
            // IDs are offset by a constant amount
            appending.add_edge(node_id + prev_node_size, next_id + prev_node_size);
        }
    }
    
    // add the paths
    for (uint64_t path_id = 0; path_id < component.path_size(); ++path_id) {
        uint64_t new_path_id = appending.add_path(component.path_name(path_id));
        for (uint64_t node_id : component.path(path_id)) {
            appending.extend_path(new_path_id, prev_node_size + node_id);
        }
    }
}

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

}

#endif /* centrolign_modify_graph_hpp */
