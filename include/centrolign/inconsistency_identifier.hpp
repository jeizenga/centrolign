#ifndef centrolign_inconsistency_identifier_hpp
#define centrolign_inconsistency_identifier_hpp

#include <vector>
#include <utility>
#include <deque>

#include "centrolign/snarls.hpp"
#include "centrolign/step_index.hpp"
#include "centrolign/structure_distances.hpp"

namespace centrolign {

/*
 * Class to identify poorly normalized graph regions
 */
class InconsistencyIdentifier {
    
public:
    InconsistencyIdentifier() = default;
    ~InconsistencyIdentifier() = default;

    // identify regions that might contain alignment inconsistencies
    // returns pairs of node IDs that represent a separable subgraph
    // the separable subgraphs are all mutually disjoint with each other
    template<class BGraph>
    std::vector<std::pair<uint64_t, uint64_t>> identify_inconsistencies(const BGraph& graph, const SentinelTableau& tableau) const;
    
    // maximum size cycle that will be identified as "tight" and flagged as an inconsistency
    size_t max_tight_cycle_size = 5000;
    
    // maximum distance along a chain that where bubbles might be identified as a bond-induced inconsistency
    size_t max_bond_inconsistency_window = 100;
    // length of sequence that is disjoint between two passes of a sequence through a chain window to be identified as
    // a bond-induced inconsistency
    size_t min_inconsistency_disjoint_length = 8;
    // length of total inserted/deleted sequence between two passes of a sequence through a chain window to be identified as
    // a bond-induced inconsistency
    size_t min_inconsistency_total_length = 50;
        
    // try to pad inconsistent regions with this much extra sequence on each side, measured in minimum distance
    size_t padding_target_min_length = 1000;
    // stop padding inconsistent regions if it would require a maximum length of this much
    size_t padding_max_length_limit = 10000;
    
protected:
    
    std::vector<std::pair<uint64_t, uint64_t>> identify_tight_cycles(const SnarlTree& snarls, const StepIndex step_index,
                                                                     const std::vector<bool>& nontrivial_left_boundary) const;
    
    std::vector<std::pair<uint64_t, uint64_t>> identify_inconsistent_bonds(const SnarlTree& snarls, const StepIndex step_index,
                                                                           const std::vector<bool>& nontrivial_left_boundary) const;
    
    template<class BGraph>
    void expand_inconsistencies(std::vector<std::pair<uint64_t, uint64_t>>& inconsistencies,
                                const BGraph& graph, const SnarlTree& snarls) const;
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
            // TODO: do i need to do anything special to handle the sentinel nodes?
            // it seems to me that they could never be non-trivial boundaries on their respective ends
            nontrivial_left_boundary[compacted_graph.back(node_id)] = true;
        }
    }
    
    auto tight_cycles = identify_tight_cycles(snarls, step_index, nontrivial_left_boundary);
    auto inconsistent_bonds = identify_inconsistent_bonds(snarls, step_index, nontrivial_left_boundary);
    {
        auto dummy = std::move(nontrivial_left_boundary);
    }
    
    // record the position of snarls in their chain so we can tell which intervals are larger than others
    std::vector<size_t> position_in_chain(snarls.structure_size());
    for (uint64_t chain_id = 0; chain_id < snarls.chain_size(); ++chain_id) {
        const auto& chain = snarls.structures_inside(chain_id);
        for (size_t i = 0; i < chain.size(); ++i) {
            position_in_chain[chain[i]] = i;
        }
    }
    
    // the furthest way snarl including which, there's an inconsistency starting at ths snarl
    std::vector<uint64_t> snarl_inconsistencies(snarls.structure_size(), -1);
    for (const auto& tight_cycle : tight_cycles) {
        snarl_inconsistencies[snarls.structure_beginning_at(tight_cycle.first)] = snarls.structure_ending_at(tight_cycle.second);
    }
    for (const auto& bond_inconsistency : inconsistent_bonds) {
        uint64_t snarl_id = snarls.structure_beginning_at(bond_inconsistency.first);
        uint64_t other_snarl_id = snarls.structure_ending_at(bond_inconsistency.second);
        if (snarl_inconsistencies[snarl_id] == -1 || position_in_chain[snarl_inconsistencies[snarl_id]] < position_in_chain[other_snarl_id]) {
            snarl_inconsistencies[snarl_id] = other_snarl_id;
        }
    }
    {
        auto dummy = std::move(position_in_chain);
    }
    
    std::vector<std::pair<uint64_t, uint64_t>> merged_inconsistencies;
    
    // top down traversal in which we stop traversing downward into snarls that have been flagged
    // inconsistent (this effectively deduplicates in favor of larger structures)
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
                if (snarl_inconsistencies[chain[i]] != -1) {
                    // this is the first node of an inconsistency that we've identified
                    if (!merged_inconsistencies.empty() &&
                        merged_inconsistencies.back().second == snarls.structure_boundaries(chain[i]).first) {
                        // it's immediately adjacent to the last snarl inteval we've flagged, so we can merge them
                        merged_inconsistencies.back().second = snarls.structure_boundaries(snarl_inconsistencies[chain[i]]).second;
                    }
                    else {
                        // add a new interval to the merged list
                        merged_inconsistencies.emplace_back(snarls.structure_boundaries(chain[i]).first,
                                                            snarls.structure_boundaries(snarl_inconsistencies[chain[i]]).second);
                    }
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
    
    // add padding for the alignment
    expand_inconsistencies(merged_inconsistencies, graph, snarls);
    
    return merged_inconsistencies;
}

template<class BGraph>
void InconsistencyIdentifier::expand_inconsistencies(std::vector<std::pair<uint64_t, uint64_t>>& inconsistencies,
                                                     const BGraph& graph, const SnarlTree& snarls) const {
    
    
    
    // all of the info that we need to track how much we've expanded each inconsitency
    struct HeapFrame {
        // minimum distance expanded
        size_t expanded_left_min = 0;
        size_t expanded_right_min = 0;
        // maximum distance expanded
        size_t expanded_left_max = 0;
        size_t expanded_right_max = 0;
        // have we hit a limit or the end of a region we can expand to
        size_t can_expand_left = true;
        size_t can_expand_right = true;
        // where have we expanded to
        uint64_t left_boundary = -1;
        uint64_t right_boundary = -1;
        // what inconsistency did this originate from
        size_t origin = -1;
        // the side we're expanding next and how much its been expanded
        std::pair<size_t, bool> frontier() const {
            assert(can_expand_left || can_expand_right);
            if ((can_expand_left && expanded_left_min < expanded_right_min) || !can_expand_right) {
                return std::make_pair(expanded_left_min, true);
            }
            else {
                return std::make_pair(expanded_right_min, false);
            }
        }
        // prioritize by the least expanded
        bool operator<(const HeapFrame& other) const {
            return frontier() > other.frontier();
        }
    };
    
    // for tracking how far we've moved along the chain
    SnarlDistances snarl_distances(snarls, graph);
    
    // initialize the heap
    std::vector<HeapFrame> heap;
    heap.reserve(inconsistencies.size());
    std::unordered_set<uint64_t> is_boundary;
    for (size_t i = 0; i < inconsistencies.size(); ++i) {
        const auto& boundaries = inconsistencies[i];
        is_boundary.insert(boundaries.first);
        is_boundary.insert(boundaries.second);
        heap.emplace_back();
        heap.back().left_boundary = boundaries.first;
        heap.back().right_boundary = boundaries.second;
        heap.back().origin = i;
    }
    std::make_heap(heap.begin(), heap.end());
    
    // expand the inconsistencies until we can't anymore
    while (!heap.empty()) {
        
        // select the least-expanded inconsistency
        std::pop_heap(heap.begin(), heap.end());
        auto& next = heap.back();
        
        if (next.frontier().second) {
            // expand left
            uint64_t next_snarl = snarls.structure_ending_at(next.left_boundary);
            if (next_snarl == -1) {
                // end of the chain
                next.can_expand_left = false;
            }
            else {
                uint64_t next_boundary = snarls.structure_boundaries(next_snarl).first;
                if (is_boundary.count(next_boundary)) {
                    // abutting another expanded region
                    next.can_expand_left = false;
                }
                else {
                    auto dists = snarl_distances.structure_min_max_dist(next_snarl);
                    if (dists.second == -1) {
                        // the next snarl contains a cycle
                        next.can_expand_left = false;
                    }
                    else {
                        size_t next_min_dist = next.expanded_left_min + dists.first - graph.label_size(next_boundary);
                        size_t next_max_dist = next.expanded_left_max + dists.second - graph.label_size(next_boundary);
                        if (next_min_dist > padding_target_min_length || next_max_dist > padding_max_length_limit) {
                            // snarl would exceed the limits
                            next.can_expand_left = false;
                        }
                        else {
                            // we can expand past this snarl
                            next.expanded_left_min = next_min_dist;
                            next.expanded_left_max = next_max_dist;
                            is_boundary.erase(next.left_boundary);
                            next.left_boundary = next_boundary;
                            is_boundary.insert(next_boundary);
                        }
                        
                    }
                }
            }
        }
        else {
            // expand right
            // TODO: extremely repetitive
            uint64_t next_snarl = snarls.structure_beginning_at(next.right_boundary);
            if (next_snarl == -1) {
                // end of the chain
                next.can_expand_right = false;
            }
            else {
                uint64_t next_boundary = snarls.structure_boundaries(next_snarl).second;
                if (is_boundary.count(next_boundary)) {
                    // abutting another expanded region
                    next.can_expand_right = false;
                }
                else {
                    auto dists = snarl_distances.structure_min_max_dist(next_snarl);
                    if (dists.second == -1) {
                        // the next snarl contains a cycle
                        next.can_expand_right = false;
                    }
                    else {
                        size_t next_min_dist = next.expanded_right_min + dists.first - graph.label_size(next_boundary);
                        size_t next_max_dist = next.expanded_right_max + dists.second - graph.label_size(next_boundary);
                        if (next_min_dist > padding_target_min_length || next_max_dist > padding_max_length_limit) {
                            // snarl would exceed the limits
                            next.can_expand_right = false;
                        }
                        else {
                            // we can expand past this snarl
                            next.expanded_right_min = next_min_dist;
                            next.expanded_right_max = next_max_dist;
                            is_boundary.erase(next.right_boundary);
                            next.right_boundary = next_boundary;
                            is_boundary.insert(next_boundary);
                        }
                        
                    }
                }
            }
        }
        
        if (!next.can_expand_left && !next.can_expand_right) {
            // we can no longer expand here, record the final result
            inconsistencies[next.origin].first = next.left_boundary;
            inconsistencies[next.origin].second = next.right_boundary;
            heap.pop_back();
        }
        else {
            // set this region up to be expanded again
            std::push_heap(heap.begin(), heap.end());
        }
    }
}

}

#endif /* centrolign_inconsistency_identifier_hpp */
