#include "centrolign/chain_merge.hpp"

namespace centrolign {

using namespace std;

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
