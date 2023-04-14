#include "centrolign/chain_merge.hpp"

namespace centrolign {

using namespace std;


vector<vector<uint64_t>> ChainMerge::chain_forward_edges() const {
    
    // intantiate the chains
    std::vector<std::vector<uint64_t>> chains(chain_size());
    for (uint64_t node_id = 0; node_id < table.size(); ++node_id) {
        uint64_t c;
        size_t idx;
        tie(c, idx) = node_to_chain[node_id];
        auto& chain = chains[c];
        while (chain.size() <= idx) {
            chain.emplace_back();
        }
        chain[idx] = node_id;
    }
    
    // identify the forward links
    std::vector<std::vector<uint64_t>> forward(table.size());
    for (uint64_t node_id = 0; node_id < table.size(); ++node_id) {
        const auto& row = table[node_id];
        for (size_t c = 0; c < row.size(); ++c) {
            if (row[c] != -1) {
                forward[chains[c][row[c]]].push_back(node_id);
            }
        }
    }
    return forward;
}

}
