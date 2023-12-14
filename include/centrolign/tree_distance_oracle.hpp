#ifndef centrolign_tree_distance_oracle_hpp
#define centrolign_tree_distance_oracle_hpp

#include <vector>
#include <cstdint>
#include <algorithm>

#include "centrolign/tree.hpp"
#include "centrolign/range_min_query.hpp"

namespace centrolign {

/*
 * Data structure providing O(1) distance queries between nodes of a tree
 */
class TreeDistanceOracle {
public:
    TreeDistanceOracle(const Tree& tree);
    TreeDistanceOracle() = default;
    ~TreeDistanceOracle() = default;
    
    inline double distance(uint64_t node_id1, uint64_t node_id2) const;
    
private:
    
    std::vector<uint64_t> euler_nodes;
    std::vector<size_t> euler_depths;
    RMQ<size_t> euler_rmq;
    std::vector<size_t> position;
    std::vector<double> depths;
};


/*
 * Inline implementations
 */

inline double TreeDistanceOracle::distance(uint64_t node_id1, uint64_t node_id2) const {
    size_t lo = position[node_id1];
    size_t hi = position[node_id2];
    if (hi < lo) {
        std::swap(lo, hi);
    }
    return depths[node_id1] + depths[node_id2] - 2 * depths[euler_nodes[euler_rmq.range_arg_min(lo, hi + 1)]];
}

}
#endif /* centrolign_tree_distance_oracle_hpp */
