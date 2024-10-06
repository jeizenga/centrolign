#ifndef centrolign_gfa_hpp
#define centrolign_gfa_hpp

#include <iostream>
#include <vector>
#include <deque>

#include "centrolign/graph.hpp"
#include "centrolign/utility.hpp"
#include "centrolign/modify_graph.hpp"

namespace centrolign {

// write as a maximally node-compacted GFA
template<class BGraph>
void write_gfa(const BGraph& graph, std::ostream& out, bool decode = true);

// write as a maximally node-compacted GFA, removing sentinel nodes
template<class BGraph>
void write_gfa(const BGraph& graph, const SentinelTableau& tableau,
               std::ostream& out, bool decode = true);

// read in a GFA and break into a base-level graph
// assumes GFA v1.1, and that segment IDs are integers, and that lines are in the order H, S, L, P
BaseGraph read_gfa(std::istream& in, bool encode = true);

/*
 * Template implementations
 */

template<class BGraph>
void write_gfa(const BGraph& graph, std::ostream& out, bool decode) {
    write_gfa_internal(graph, nullptr, out, decode);
}


template<class BGraph>
void write_gfa(const BGraph& graph, const SentinelTableau& tableau,
               std::ostream& out, bool decode) {
    write_gfa_internal(graph, &tableau, out, decode);
}

template<class BGraph>
void write_gfa_internal(const BGraph& graph, const SentinelTableau* tableau,
                        std::ostream& out, bool decode) {
    
    static const bool debug = false;
    
    std::vector<bool> path_begin(graph.node_size(), false);
    std::vector<bool> path_end(graph.node_size(), false);
    std::vector<bool> compacted_end(graph.node_size(), false);
    std::vector<uint64_t> compacted_id(graph.node_size(), -1);
    
    for (uint64_t path_id = 0; path_id < graph.path_size(); ++path_id) {
        path_begin[graph.path(path_id).front()] = true;
        path_end[graph.path(path_id).back()] = true;
    }
    
    
    // header
    out << "H\tVN:Z:1.0\n";
    
    // construct the nodes
    uint64_t next_compacted_id = 1;
    for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
        if (compacted_id[node_id] != -1) {
            continue;
        }
        if (tableau && (node_id == tableau->src_id || node_id == tableau->snk_id)) {
            continue;
        }
        if (debug) {
            std::cerr << "forming a compacted node seeded from " << node_id << '\n';
        }
        
        std::vector<uint64_t> compacted(1, node_id);
        // walk left
        while (!path_begin[compacted.back()]
               && graph.previous_size(compacted.back()) == 1
               && !path_end[graph.previous(compacted.back()).front()]
               && graph.next_size(graph.previous(compacted.back()).front()) == 1
               && (!tableau || (graph.previous(compacted.back()).front() != tableau->src_id &&
                                graph.previous(compacted.back()).front() != tableau->snk_id))) {
            compacted.push_back(graph.previous(compacted.back()).front());
        }
        std::reverse(compacted.begin(), compacted.end());
        // walk right
        while (!path_end[compacted.back()]
               && graph.next_size(compacted.back()) == 1
               && !path_begin[graph.next(compacted.back()).front()]
               && graph.previous_size(graph.next(compacted.back()).front()) == 1
               && (!tableau || (graph.next(compacted.back()).front() != tableau->src_id &&
                                graph.next(compacted.back()).front() != tableau->snk_id))) {
            compacted.push_back(graph.next(compacted.back()).front());
        }
        
        if (debug) {
            std::cerr << "walked compacted node:\n";
            for (auto n : compacted) {
                std::cerr << ' ' << n;
            }
            std::cerr << '\n';
        }
        
        std::string seq;
        for (auto n : compacted) {
            compacted_id[n] = next_compacted_id;
            seq.push_back(decode ? decode_base(graph.label(n)) : graph.label(n));
        }
        
        compacted_end[compacted.back()] = true;
        
        // write the S line
        out << "S\t" << next_compacted_id++ << '\t' << seq << '\n';
    }
    
    // edges
    for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
        if (!compacted_end[node_id] || (tableau && (node_id == tableau->src_id || node_id == tableau->snk_id))) {
            continue;
        }
        for (auto next_id : graph.next(node_id)) {
            if (tableau && (next_id == tableau->src_id || next_id == tableau->snk_id)) {
                continue;
            }
            // write the L line
            out << "L\t" << compacted_id[node_id] << "\t+\t" << compacted_id[next_id] << "\t+\t*\n";
        }
    }
    
    // paths
    for (uint64_t path_id = 0; path_id < graph.path_size(); ++path_id) {
        out << "P\t" << graph.path_name(path_id) << '\t';
        bool write_next = true;
        bool first = true;
        for (auto node_id : graph.path(path_id)) {
            if (tableau && (node_id == tableau->src_id || node_id == tableau->snk_id)) {
                // note: as is, this shouldn't ever happen, but hopefully this future-proofs
                // any future decisions to include sentinel nodes in paths
                continue;
            }
            if (write_next) {
                if (!first) {
                    out << ',';
                }
                out << compacted_id[node_id] << '+';
                first = false;
            }
            write_next = compacted_end[node_id];
        }
        out << "\t*\n";
    }
}

}

#endif /* centrolign_gfa_hpp */
