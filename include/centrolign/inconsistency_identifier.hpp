#ifndef centrolign_inconsistency_identifier_hpp
#define centrolign_inconsistency_identifier_hpp

#include <vector>
#include <utility>
#include <deque>

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
    
    // maximum size cycle that will be identified as "tight" and flagged
    size_t max_tight_cycle_size = 5000;
    
    // maximum distance along a chain that where bubbles might be identified as a bond-induced inconsistency
    size_t max_bond_inconsistency_window = 100;
    // length of sequence that is disjoint between two passes of a sequence through a chain window to be identified as
    // a bond-induced inconsistency
    size_t min_inconsistency_disjoint_length = 8;
    // length of total inserted/deleted sequence between two passes of a sequence through a chain window to be identified as
    // a bond-induced inconsistency
    size_t min_inconsistency_total_length = 50;
    
protected:
    
    std::vector<std::pair<uint64_t, uint64_t>> identify_tight_cycles(const SnarlTree& snarls, const StepIndex step_index,
                                                                     const std::vector<bool>& nontrivial_left_boundary) const;
    
    std::vector<std::pair<uint64_t, uint64_t>> identify_inconsistent_bonds(const SnarlTree& snarls, const StepIndex step_index,
                                                                           const std::vector<bool>& nontrivial_left_boundary) const;
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
    {
        CompactedGraph compacted_graph(graph);
        for (uint64_t node_id = 0; node_id < compacted_graph.node_size(); ++node_id) {
            nontrivial_left_boundary[compacted_graph.back(node_id)] = true;
        }
    }
    // TODO: do i need to do anything special to handle the sentinel nodes?
    // it seems to me that they could never be non-trivial boundaries on their respective ends
    
    auto tight_cycles = identify_tight_cycles(snarls, step_index, nontrivial_left_boundary);
    auto inconsistent_bonds = identify_inconsistent_bonds(snarls, step_index, nontrivial_left_boundary);
    
    std::vector<std::pair<uint64_t, uint64_t>> snarl_inconsistencies(snarls.structure_size(), std::pair<uint64_t, uint64_t>(-1, -1));
    // note: we add these first so that they can possible be overwritten by inconsistent bonds, which possibly
    // encompass multiple snarls
    for (const auto& tight_cycle : tight_cycles) {
        snarl_inconsistencies[snarls.structure_beginning_at(tight_cycle.first)] = tight_cycle;
    }
    for (const auto& bond_inconsistency : inconsistent_bonds) {
        snarl_inconsistencies[snarls.structure_beginning_at(bond_inconsistency.first)] = bond_inconsistency;
    }
    
    std::vector<std::pair<uint64_t, uint64_t>> merged_inconsistencies;
    
    std::deque<std::pair<uint64_t, bool>> queue;
    for (uint64_t chain_id = 0; chain_id < snarls.chain_size(); ++chain_id) {
        if (snarls.structure_containing(chain_id) == -1) {
            queue.emplace_back(chain_id, true);
        }
    }
    while (!queue.empty()) {
        auto feature = queue.front();
        queue.pop_front();
        if (feature.second) {
            // chain
            const auto& chain = snarls.structures_inside(feature.first);
            for (size_t i = 0; i < chain.size(); ++i) {
                if (snarl_inconsistencies[chain[i]].first != -1) {
                    // this is the first node of an inconsistency that we've identified
                    merged_inconsistencies.emplace_back(snarl_inconsistencies[chain[i]]);
                    // advance to the end of the window, if necessary
                    while (chain[i] != snarls.structure_ending_at(merged_inconsistencies.back().second)) {
                        ++i;
                    }
                }
                else {
                    // we can still find inconsistencies in this snarl's children
                    queue.emplace_back(chain[i], false);
                }
            }
        }
        else {
            // snarl
            for (auto chain_id : snarls.chains_inside(feature.first)) {
                queue.emplace_back(chain_id, true);
            }
        }
    }
    
    return merged_inconsistencies;
}



}

#endif /* centrolign_inconsistency_identifier_hpp */
