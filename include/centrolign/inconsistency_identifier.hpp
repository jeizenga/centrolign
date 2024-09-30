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
    std::vector<std::pair<uint64_t, uint64_t>> identify_tight_cycles(const BGraph& graph, const SentinelTableau& tableau) const;
    
//    template<class BGraph>
//    std::vector<std::pair<uint64_t, uint64_t>> identify_inconsistent_bonds(const BGraph& graph, const SentinelTableau& tableau) const;
    
    size_t max_tight_cycle_size = 5000;
};





/**
 * Template and inline implentations
 */
//template<class BGraph>
//void InconsistencyIdentifier::identify_inconsistencies(const BGraph& graph, const SentinelTableau& tableau) const {
//
//    CactusGraph cactus_graph(graph, tableau);
//    CactusTree cactus_tree(cactus_graph);
//
//    // do a postorder traversal
//    std::vector<std::pair<uint64_t>> stack(1, std::make_pair(cactus_tree.get_root(), false));
//    while (!stack.empty()) {
//        if (!stack.back().second) {
//            stack.back().second = true;
//            for (auto next_id : cactus_tree.get_children(node_id)) {
//                stack.emplace_back(next_id, false);
//            }
//        }
//        else {
//            stack.pop_back();
//        }
//    }
//}


template<class BGraph>
std::vector<std::pair<uint64_t, uint64_t>> InconsistencyIdentifier::identify_tight_cycles(const BGraph& graph, const SentinelTableau& tableau) const {
    
    static const bool debug = false;
    if (debug) {
        std::cerr << "identifying tight cycles with size limit " << max_tight_cycle_size << '\n';
    }

    SnarlTree snarls(graph, tableau);
    StepIndex step_index(graph);
    
    std::vector<bool> chain_blocked(snarls.chain_size(), false);
    std::vector<bool> snarl_blocked(snarls.structure_size(), false);
    
    std::vector<std::list<uint64_t>> chain_cyclic_descendents(snarls.chain_size());
    std::vector<std::list<uint64_t>> snarl_cyclic_descendents(snarls.structure_size());
    
    std::unordered_set<std::pair<uint64_t, bool>> tight_cycle_features;
    
    for (auto feature : snarls.postorder()) {
        
        if (debug) {
            std::cerr << "at feature " << feature.first << " " << feature.second << " in postorder\n";
        }
        
        uint64_t start_node, end_node;
        if (feature.second) {
            // is a chain
            if (chain_blocked[feature.first]) {
                // already too big in a descendent
                uint64_t snarl_id = snarls.structure_containing(feature.first);
                if (snarl_id != -1) {
                    snarl_blocked[snarl_id] = true;
                }
                if (debug) {
                    std::cerr << "blocked\n";
                }
                continue;
            }
            
            start_node = snarls.structure_boundaries(snarls.structures_inside(feature.first).front()).first;
            end_node = snarls.structure_boundaries(snarls.structures_inside(feature.first).back()).second;
        }
        else {
            // is a snarl
            if (snarl_blocked[feature.first]) {
                // already too big in a descendent
                chain_blocked[snarls.chain_containing(feature.first)] = true;
                if (debug) {
                    std::cerr << "blocked\n";
                }
                continue;
            }
            std::tie(start_node, end_node) = snarls.structure_boundaries(feature.first);
        }
        
        std::unordered_map<uint64_t, std::pair<std::vector<size_t>, std::vector<size_t>>> path_positions;
        for (const auto& step : step_index.path_steps(start_node)) {
            path_positions[step.first].first.push_back(step.second);
        }
        for (const auto& step : step_index.path_steps(end_node)) {
            path_positions[step.first].second.push_back(step.second);
        }
        
        size_t max_path_size = 0;
        for (auto it = path_positions.begin(); it != path_positions.end(); ++it) {
            // make sure that multiple passes are in increasing order, although i think the object already provides this
            if (!std::is_sorted(it->second.first.begin(), it->second.first.end())) {
                std::sort(it->second.first.begin(), it->second.first.end());
            }
            if (!std::is_sorted(it->second.second.begin(), it->second.second.end())) {
                std::sort(it->second.second.begin(), it->second.second.end());
            }
            
            for (size_t i = 0; i < it->second.first.size(); ++i) {
                max_path_size = std::max(max_path_size, it->second.second[i] - it->second.first[i]);
            }
        }
        
        if (debug) {
            std::cerr << "max path size " << max_path_size << " between boundary nodes " << start_node << " and " << end_node << "\n";
        }
        
        if (max_path_size > max_tight_cycle_size) {
            // too big, block the parent
            if (feature.second) {
                uint64_t snarl_id = snarls.structure_containing(feature.first);
                if (snarl_id != -1) {
                    snarl_blocked[snarl_id] = true;
                }
            }
            else {
                chain_blocked[feature.first] = true;
            }
        }
        else if (!feature.second) {
            // a snarl, so we need to check for cycles
            if (!snarls.net_graph_is_acyclic(feature.first)) {
                // this cycle subsumes any contained cycles
                if (debug) {
                    std::cerr << "snarl is cyclic, clearing out existing cyclic descendents:\n";
                    for (auto snarl_id : snarl_cyclic_descendents[feature.first]) {
                        std::cerr << '\t' << snarl_id << '\n';
                    }
                }
                snarl_cyclic_descendents[feature.first].clear();
                snarl_cyclic_descendents[feature.first].emplace_back(feature.first);
            }
        }
        
        // pass the list of cyclic components up to the parent
        if (feature.second) {
            // chain
            uint64_t snarl_id = snarls.structure_containing(feature.first);
            if (snarl_id != -1) {
                auto& parent_list = snarl_cyclic_descendents[snarl_id];
                parent_list.splice(parent_list.end(), chain_cyclic_descendents[feature.first]);
            }
        }
        else {
            // snarl
            uint64_t chain_id = snarls.chain_containing(feature.first);
            auto& parent_list = chain_cyclic_descendents[chain_id];
            parent_list.splice(parent_list.end(), snarl_cyclic_descendents[feature.first]);
        }
    }
    
    std::vector<std::pair<uint64_t, uint64_t>> tight_cycles;
    // the lists stop propagating upward when they get blocked for being too big, so we collect them here
    for (auto result_lists_ptr : {&chain_cyclic_descendents, &snarl_cyclic_descendents}) {
        for (auto result_list : *result_lists_ptr) {
            for (uint64_t snarl_id : result_list) {
                tight_cycles.emplace_back(snarls.structure_boundaries(snarl_id));
            }
        }
    }
    
    return tight_cycles;
}

}

#endif /* centrolign_inconsistency_identifier_hpp */
