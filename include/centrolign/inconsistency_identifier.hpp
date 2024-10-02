#ifndef centrolign_inconsistency_identifier_hpp
#define centrolign_inconsistency_identifier_hpp

#include <vector>
#include <utility>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <list>

#include "centrolign/snarls.hpp"
#include "centrolign/step_index.hpp"

namespace centrolign {

/*
 * Class to identify poorly normalized graph regions
 */
class InconsistencyIdentifier {
    
public:
    InconsistencyIdentifier() = default;
    ~InconsistencyIdentifier() = default;

    template<class BGraph>
    std::vector<std::pair<uint64_t, uint64_t>> identify_inconsistencies(const BGraph& graph, const SentinelTableau& tableau) const;
    
    size_t max_tight_cycle_size = 5000;
    
    size_t max_bond_inconsistency_window = 100;
    size_t min_inconsistency_disjoint_length = 8;
    size_t min_inconsistency_total_length = 50;
    
protected:
    
    std::vector<std::pair<uint64_t, uint64_t>> identify_tight_cycles(const SnarlTree& snarls, const StepIndex step_index,
                                                                     const std::vector<bool>& nontrivial_left_boundary,
                                                                     const std::vector<bool>& nontrivial_right_boundary) const;
    
    std::vector<std::pair<uint64_t, uint64_t>> identify_inconsistent_bonds(const SnarlTree& snarls, const StepIndex step_index,
                                                                           const std::vector<bool>& nontrivial_left_boundary,
                                                                           const std::vector<bool>& nontrivial_right_boundary) const;
};





/**
 * Template and inline implentations
 */
template<class BGraph>
std::vector<std::pair<uint64_t, uint64_t>> InconsistencyIdentifier::identify_inconsistencies(const BGraph& graph, const SentinelTableau& tableau) const {

    SnarlTree snarls(graph, tableau);
    StepIndex step_index(graph);
    
    // label the nodes that can be boundaries of non-trivial snarls
    std::vector<bool> nontrivial_left_boundary(graph.node_size(), false);
    std::vector<bool> nontrivial_right_boundary(graph.node_size(), false);
    {
        CompactedGraph compacted_graph(graph);
        for (uint64_t node_id = 0; node_id < compacted_graph.node_size(); ++node_id) {
            nontrivial_right_boundary[compacted_graph.front(node_id)] = true;
            nontrivial_left_boundary[compacted_graph.back(node_id)] = true;
        }
    }
    // TODO: do i need to do anything special to handle the sentinel nodes?
    // it seems to me that they could never be non-trivial boundaries on their respective ends
    
    auto inconsistent_bonds = identify_inconsistent_bonds(snarls, step_index, nontrivial_left_boundary, nontrivial_right_boundary);
    auto tight_cycles = identify_tight_cycles(snarls, step_index, nontrivial_left_boundary, nontrivial_right_boundary);
    
    // TODO: merge and return
    
    return std::vector<std::pair<uint64_t, uint64_t>>();
}



}

#endif /* centrolign_inconsistency_identifier_hpp */
