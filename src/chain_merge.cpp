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

}
