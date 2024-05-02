#include "centrolign/chain_merge.hpp"

namespace centrolign {

using namespace std;


vector<vector<pair<uint64_t, uint64_t>>> ChainMerge::chain_forward_edges(const std::vector<bool>* to_nodes,
                                                                         const std::vector<bool>* from_nodes) const {
    
    // identify the forward links
    vector<vector<pair<uint64_t, uint64_t>>> forward(table.size());
    for (uint64_t node_id = 0; node_id < table.size(); ++node_id) {
        if (to_nodes && !(*to_nodes)[node_id]) {
            continue;
        }
        const auto& row = table[node_id];
        for (size_t c = 0; c < row.size(); ++c) {
            if (row[c] != -1) {
                auto from_id = chains[c][row[c]];
                if (!from_nodes || (*from_nodes)[from_id]) {
                    forward[from_id].emplace_back(node_id, c);
                }
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
            + table.capacity() * sizeof(decltype(table)::value_type)
            + table_size * sizeof(decltype(table)::value_type::value_type)
            + node_to_chain.capacity() * sizeof(decltype(node_to_chain)::value_type)
            + chains.capacity() * sizeof(decltype(chains)::value_type)
            + chain_size * sizeof(decltype(chains)::value_type::value_type));
}

}
