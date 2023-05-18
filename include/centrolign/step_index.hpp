#ifndef centrolign_step_index_hpp
#define centrolign_step_index_hpp

#include <vector>
#include <utility>
#include <cstdint>

namespace centrolign {

/*
 * Memoizes node ID to path step queries
 */
class StepIndex {
public:
    template<class Graph>
    StepIndex(const Graph& graph);
    StepIndex() = default;
    ~StepIndex() = default;
    
    // returns a list of (path ID, step index) pairs that are on the node
    inline const std::vector<std::pair<uint64_t, size_t>>& path_steps(uint64_t node_id) const;
    
private:
    
    std::vector<std::vector<std::pair<uint64_t, size_t>>> steps;
};

/*
 * Template and inline implementations
 */

template<class Graph>
StepIndex::StepIndex(const Graph& graph) : steps(graph.node_size()) {
    for (uint64_t path_id = 0; path_id < graph.path_size(); ++path_id) {
        size_t i = 0;
        for (auto node_id : graph.path(path_id)) {
            steps[node_id].emplace_back(path_id, i++);
        }
    }
}

inline const std::vector<std::pair<uint64_t, size_t>>& StepIndex::path_steps(uint64_t node_id) const {
    return steps[node_id];
}

}

#endif /* centrolign_step_index_hpp */
