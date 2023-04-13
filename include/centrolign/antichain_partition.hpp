#ifndef centrolign_antichain_partition_hpp
#define centrolign_antichain_partition_hpp

#include <vector>
#include <iostream>
#include <cassert>

namespace centrolign {

// Modified Kahn's algorithm, returns an antichain partition consisting
// of an antichain ID indexed by node ID. This partition additionally
// has the property that x \pred y ==> AC(x) < AC(y)
template<class Graph>
std::vector<size_t> antichain_partition(const Graph& graph) {
    
    static const bool debug_partition = false;
    
    std::vector<size_t> partition(graph.node_size(), -1);
    
    // find sources
    std::vector<uint64_t> next_antichain;
    std::vector<size_t> in_degree(graph.node_size());
    for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
        in_degree[node_id] = graph.previous_size(node_id);
        if (in_degree[node_id] == 0) {
            if (debug_partition) {
                std::cerr << "found node in minimum set " << node_id << '\n';
            }
            next_antichain.push_back(node_id);
        }
    }
    
    size_t antichain_idx = 0;
    while (!next_antichain.empty()) {
        
        if (debug_partition) {
            std::cerr << "antichain " << antichain_idx << " contains nodes:\n";
            for (auto node_id : next_antichain) {
                std::cerr << ' ' << node_id;
            }
            std::cerr << '\n';
        }
        
        // assign the next antichain and identify the next min set
        std::vector<uint64_t> antichain = std::move(next_antichain);
        next_antichain.clear();
        for (auto node_id : antichain) {
            partition[node_id] = antichain_idx;
            for (uint64_t next_id : graph.next(node_id)) {
                --in_degree[next_id];
                if (in_degree[next_id] == 0) {
                    next_antichain.push_back(next_id);
                }
            }
        }
        ++antichain_idx;
    }
    
    return partition;
}

}

#endif /* centrolign_antichain_partition_hpp */
