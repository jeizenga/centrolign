#ifndef centrolign_gfa_hpp
#define centrolign_gfa_hpp

#include <iostream>
#include <vector>
#include <deque>

#include "centrolign/graph.hpp"

namespace centrolign {

// write as a maximally node-compacted GFA
template<class BGraph>
void write_gfa(const BGraph& graph, std::ostream& out);


/*
 * Template implemenatations
 */

template<class BGraph>
void write_gfa(const BGraph& graph, std::ostream& out) {
    
    static const bool debug = false;
    
    std::vector<bool> path_begin(graph.node_size(), false);
    std::vector<bool> path_end(graph.node_size(), false);
    std::vector<bool> written(graph.node_size(), false);
    std::vector<uint64_t> compacted_id(graph.node_size(), -1);
    
    for (uint64_t path_id = 0; path_id < graph.path_size(); ++path_id) {
        path_begin[graph.path(path_id).front()] = true;
        path_end[graph.path(path_id).back()] = true;
    }
    
    
    // header
    out << "H\tVN:Z:1.2\n";
    
    // construct the nodes
    uint64_t next_compacted_id = 1;
    for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
        if (written[node_id]) {
            continue;
        }
        if (debug) {
            std::cerr << "forming a compacted node seeded from " << node_id << '\n';
        }
        
        std::deque<uint64_t> compacted(1, node_id);
        // walk left
        while (!path_begin[compacted.front()]
               && graph.previous_size(compacted.front()) == 1
               && !path_end[graph.previous(compacted.front()).front()]
               && graph.next_size(graph.previous(compacted.front()).front()) == 1) {
            compacted.push_front(graph.previous(compacted.front()).front());
        }
        // walk right
        while (!path_end[compacted.back()]
               && graph.next_size(compacted.back()) == 1
               && !path_begin[graph.next(compacted.back()).front()]
               && graph.previous_size(graph.next(compacted.back()).front()) == 1) {
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
            written[n] = true;
            seq.push_back(graph.label(n));
        }
        
        // write the S line
        out << "S\t" << next_compacted_id++ << '\t' << seq << '\n';
    }
    
    // edges
    for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
        for (auto next_id : graph.next(node_id)) {
            if (compacted_id[next_id] != compacted_id[node_id]) {
                // write the L line
                out << "L\t" << compacted_id[node_id] << "\t+\t" << compacted_id[next_id] << "\t+\t*\n";
            }
        }
    }
    
    // paths
    for (uint64_t path_id = 0; path_id < graph.path_size(); ++path_id) {
        out << "P\t" << graph.path_name(path_id) << '\t';
        uint64_t curr_id = -1;
        for (auto node_id : graph.path(path_id)) {
            if (compacted_id[node_id] != curr_id) {
                if (curr_id != -1) {
                    out << ',';
                }
                curr_id = compacted_id[node_id];
                out << curr_id;
            }
        }
        out << "\t*\n";
    }
}

}

#endif /* centrolign_gfa_hpp */
