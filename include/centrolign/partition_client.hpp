#ifndef centrolign_partition_client_hpp
#define centrolign_partition_client_hpp

#include <vector>
#include <utility>
#include <algorithm>

namespace centrolign {


/*
 * Class to share part of the algorithms with other partitioning code
 */
class PartitionClient {
public:
    PartitionClient() = default;
    ~PartitionClient() = default;
protected:
    template<class T>
    std::vector<std::pair<size_t, size_t>> traceback(const std::vector<std::pair<T, T>>& dp,
                                                     const std::vector<size_t>& backpointer,
                                                     size_t tb_idx) const;
};

/*
 * Template implementations
 */

template<class T>
std::vector<std::pair<size_t, size_t>> PartitionClient::traceback(const std::vector<std::pair<T, T>>& dp,
                                                                  const std::vector<size_t>& backpointer, size_t tb_idx) const {
    
    std::vector<std::pair<size_t, size_t>> partition;
    
    bool in_interval = true;
    while (tb_idx > 0) {
        if (in_interval) {
            size_t prev = backpointer[tb_idx];
            partition.emplace_back(prev, tb_idx);
            tb_idx = prev;
            in_interval = false;
        }
        else {
            in_interval = (dp[tb_idx].first == dp[tb_idx - 1].second);
            --tb_idx;
        }
    }
    
    // put in forward order
    std::reverse(partition.begin(), partition.end());
    
    return partition;
}

}

#endif /* centrolign_partition_client_hpp */
