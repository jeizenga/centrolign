#include "centrolign/chain_merge.hpp"

namespace centrolign {

using namespace std;


vector<vector<pair<uint64_t, uint64_t>>> ChainMerge::chain_forward_edges() const {
    
    // identify the forward links
    vector<vector<pair<uint64_t, uint64_t>>> forward(table.size());
    for (uint64_t node_id = 0; node_id < table.size(); ++node_id) {
        const auto& row = table[node_id];
        for (size_t c = 0; c < row.size(); ++c) {
            if (row[c] != -1) {
                forward[chains[c][row[c]]].emplace_back(node_id, c);
            }
        }
    }
    return forward;
}

size_t ChainMerge::memory_size() const {
    
    size_t table_size = 0;
    for (const auto& row : table) {
        table_size += row.capacity();
    }
    size_t chain_size = 0;
    for (const auto& chain : chains) {
        chain_size += chain.capacity();
    }
    
    return (sizeof(table) + sizeof(node_to_chain) + sizeof(chains)
            + table.capacity() * sizeof(std::vector<size_t>)
            + table_size * sizeof(size_t)
            + node_to_chain.capacity() * sizeof(std::pair<uint64_t, size_t>)
            + chains.capacity() * sizeof(std::vector<uint64_t>)
            + chain_size * sizeof(uint64_t));
}

}
