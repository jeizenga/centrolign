#include "centrolign/path_merge.hpp"

namespace centrolign {

using namespace std;


vector<vector<uint64_t>> PathMerge::chain_forward_edges() const {
    
    // TODO: this would be simpler if i took the graph as input
    
    // reconstruct the paths
    std::vector<std::vector<uint64_t>> paths(chain_size());
    for (uint64_t node_id = 0; node_id < table.size(); ++node_id) {
        uint64_t p = path_head[node_id];
        while (p != -1) {
            auto& path = paths[p];
            size_t idx;
            tie(idx, p) = index_on_path[node_id][p];
            while (path.size() <= idx) {
                path.emplace_back();
            }
            path[idx] = node_id;
        }
    }
    
    // identify the forward links
    std::vector<std::vector<uint64_t>> forward(table.size());
    for (uint64_t node_id = 0; node_id < table.size(); ++node_id) {
        const auto& row = table[node_id];
        for (size_t p = 0; p < row.size(); ++p) {
            if (row[p] != -1) {
                forward[paths[p][row[p]]].push_back(node_id);
            }
        }
    }
    return forward;
}

}
