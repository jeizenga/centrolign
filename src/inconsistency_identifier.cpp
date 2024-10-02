#include "centrolign/inconsistency_identifier.hpp"

#include "centrolign/utility.hpp"

namespace centrolign {

std::vector<std::pair<uint64_t, uint64_t>> InconsistencyIdentifier::identify_inconsistent_bonds(const SnarlTree& snarls, const StepIndex step_index,
                                                                                                const std::vector<bool>& nontrivial_left_boundary,
                                                                                                const std::vector<bool>& nontrivial_right_boundary) const {
    
    auto query_path_positions = [&](uint64_t node_id) -> std::unordered_map<uint64_t, std::vector<size_t>> {
        std::unordered_map<uint64_t, std::vector<size_t>> path_positions;
        for (const auto& step : step_index.path_steps(node_id)) {
            path_positions[step.first].push_back(step.second);
        }
        // make sure that they're in ascending order (even though i think the data structure guarantees this)
        for (auto it = path_positions.begin(); it != path_positions.end(); ++it) {
            if (it->second.size() != 1 && !std::is_sorted(it->second.begin(), it->second.end())) {
                std::sort(it->second.begin(), it->second.end());
            }
        }
        return path_positions;
    };
    
    // TODO: I think I actually want something more like topological min distance, and a rather tight window...
    // I mostly want indels that overlap or are very close to each other
    // maybe measure it in just snarl boundary nodes
    // minimum distance over all two-pass paths
    auto min_path_distance = [&](const std::unordered_map<uint64_t, std::vector<size_t>>& path_positions_left,
                                 const std::unordered_map<uint64_t, std::vector<size_t>>& path_positions_right,
                                 const std::unordered_map<uint64_t, std::vector<size_t>>& multipass_intervals) -> size_t {
        size_t dist = -1;
        for (const auto& path_record : multipass_intervals) {
            uint64_t path_id = path_record.first;
            auto it1 = path_positions_left.find(pos_record.first);
            auto it2 = path_positions_right.find(pos_record.first);
            assert(it1->second.size() == it2->second.size());
            for (size_t i = 0; i < it1->second.size(); ++i) {
                dist = std:min(dist, it2->second[i] - it1->second[i]);
            }
        }
        return dist;
    };
    
    
    // maximum distance over any pass of a specific path
    auto max_path_distance = [&](const std::unordered_map<uint64_t, std::vector<size_t>>& path_positions_left,
                                 const std::unordered_map<uint64_t, std::vector<size_t>>& path_positions_right,
                                 uint64_t path_id) -> size_t {
        size_t dist = 0;
        auto it1 = path_positions_left.find(path_id);
        auto it2 = path_positions_right.find(path_id);
        assert(it1->second.size() == it2->second.size());
        for (size_t i = 0; i < it1->second.size(); ++i) {
            dist = std:max(dist, it2->second[i] - it1->second[i]);
        }
        return dist;
    };
    
    // median distance across all paths
    auto median_path_distance = [&](const std::unordered_map<uint64_t, std::vector<size_t>>& path_positions_left,
                                    const std::unordered_map<uint64_t, std::vector<size_t>>& path_positions_right) -> size_t {
        std::vector<size_t> dists;
        dists.reserve(path_positions_left.size());
        for (const auto& pos_record : path_positions_left) {
            auto it = path_positions_right.find(pos_record.first);
            assert(pos_record.second.size() == it->second.size());
            for (size_t i = 0; i < pos_record.second.size(); ++i) {
                dists.push_back(it->second[i] - pos_record.second[i]);
            }
        }
        return median(dists.begin(), dists.end());
    };
    
    // for each chain, the closed intervals of the snarls identified as supicious
    std::vector<std::vector<std::pair<size_t, size_t>>> chain_suspicious_intervals(snarls.chain_size());
    
    for (auto feature : snarls.postorder()) {
        if (feature.second) {
            // chain
            
            // look for a non-trivial snarl in this chain
            const auto& chain = snarls.structures_inside(feature.first);
            size_t window_start = 0;
            while (window_start < chain.size()) {
                if (nontrivial_left_boundary[snarls.structure_boundaries(chain[window_start]).first] &&
                    nontrivial_right_boundary[snarls.structure_boundaries(chain[window_start]).second]) {
                    break;
                }
                ++window_start;
            }
            
            if (window_start == chain.size()) {
                // there were no non-trivial snarls on this chain
                continue;
            }
            auto left_path_positions = query_path_positions(snarls.structure_boundaries(chain[window_start]).first);
            
            // create a look-up structure for which pass of a path a position is on
            auto multipass_intervals = query_path_positions(snarls.structure_boundaries(chain.front()).first);
            for (auto it = multipass_intervals.begin() it != multipass_intervals.end(); ) {
                auto prev = it++;
                // remove paths with a single interval
                if (prev->second.size() == 1) {
                    multipass_intervals.erase(prev);
                }
            }
            
            if (multipass_intervals.empty()) {
                // there are no paths that traverse this chain multiple times, so nothing to identify
                continue;
            }
            
            // shrink the buckets, since we may have deleted many paths and we will iterate over this structure
            multipass_intervals.rehash(multipass_intervals.size());
            // add the right ends of the intervals
            for (const auto& pos_record : query_path_positions(snarls.structure_boundaries(chain.back()).second)) {
                auto it = multipass_intervals.find(pos_record.first)
                if (it != multipass_intervals.end()) {
                    for (auto pos : pos_record.second) {
                        it->second.push_back(pos);
                    }
                }
                // sort for fast lookup
                std::sort(it->second.begin(), it->second.end());
            }
            
            // return the ordinal value of the path interval
            auto identify_pass = [&](uint64_t path_id, size_t pos) -> size_t {
                auto it = multipass_intervals.find(path_id);
                if (it == multipass_intervals.end()) {
                    return -1;
                }
                size_t i = std::upper_bound(it->second.begin(), it->second.end(), pos) - it->second.begin();
                assert(i % 2 == 1); // make sure that the position is inside the intervals
                return i / 2;
            };
            
            size_t window_end = window_start;
            size_t last_nontrivial = window_end;
            auto right_path_positions = query_path_positions(snarls.structure_boundaries(chain[window_end]).second);
            size_t dist = min_path_distance(left_path_positions, right_path_positions, multipass_intervals);
            
            while (window_end < chain.size()) {
                
                // walk the right end as far out as we can
                while (window_end < chain.size() && dist + (window_end - last_nontrivial) < max_bond_inconsistency_window) {
                    // TODO: this can lead to querying the same path positions multiple times
                    if (window_end + 1 < chain.size() &&
                        nontrivial_right_boundary[snarls.structure_boundaries(chain[window_end + 1]).second]) {
                        // proceeding to the next entry includes passing through an entire non-trivial snarl
                        
                        auto next_right_path_positions = query_path_positions(snarls.structure_boundaries(chain[window_end + 1]).second);
                        size_t next_dist = min_path_distance(left_path_positions, next_right_path_positions, multipass_intervals);
                        
                        if (next_dist < max_bond_inconsistency_window) {
                            // we can still cross this snarl without breaking the distance limit
                            ++window_end;
                            last_nontrivial = window_end;
                            dist = next_dist;
                            right_path_positions = std::move(next_right_path_positions);
                        }
                        else {
                            // we can't go through this snarl within the distance limit
                            break;
                        }
                    }
                    else {
                        // extend by one base
                        ++window_end;
                    }
                }
                
                if (dist < max_bond_inconsistency_window) {
                    // the window is within the limit
                    // note: it's necessary to check this for the case when a single snarl is too large on its own
                    
                    // count up the length of child chains, subdivivided by which
                    std::unordered_map<uint64_t, std::map<std::vector<bool>, size_t>> pass_set_length;
                    for (size_t i = window_start; i <= window_end; ++i) {
                        for (uint64_t chain_id : snarls.chains_inside(chain[i])) {
                            const auto& child_chain = snarls.structures_inside(chain_id);
                            uint64_t start_node_id = snarls.structure_boundaries(child_chain.front()).first;
                            uint64_t end_node_id = snarls.structure_boundaries(child_chain.back()).second;
                            
                            auto chain_start_path_positions = query_path_positions(start_node_id);
                            auto chain_end_path_positions = query_path_positions(end_node_id);
                            
                            for (const auto& pass_record : multipass_intervals) {
                                size_t length;
                                std::vector<bool> which_passes(pass_record.second.size() / 2, false); // note: contains both begin and ends of intervals
                                auto it = chain_start_path_positions.find(pass_record.first);
                                // TODO: are these good definitions for the chain lengths?
                                if (it == chain_start_path_positions.end()) {
                                    // get the typical length of the non-duplicated
                                    length = median_path_distance(chain_start_path_positions, chain_end_path_positions);
                                }
                                else {
                                    // we take the longest allele among the duplicated
                                    length = max_path_distance(chain_start_path_positions, chain_end_path_positions, pass_record.first);
                                    for (size_t position : it->second) {
                                        which_passes[identify_pass(it->first, position)] = true;
                                    }
                                }
                                pass_set_length[pass_record.first][std::move(which_passes)] += length;
                            }
                        }
                    }
                    
                    bool is_suspicious = false;
                    for (auto it = pass_set_length.being(), end = pass_set_length.end(); it != end && !is_suspicious; ++it) {
                        const auto& path_pass_set_length = *it;
                        // try to find a pair of passes that looks like it might be represented inconsistently over a bond
                        // TODO: could possibly speed this up by indexing pass -> pass set iterator
                        size_t num_passes = path_pass_set_length.second.begin()->first.size();
                        for (size_t pass1 = 0; pass1 < num_passes && !is_suspicious; ++pass1) {
                            for (size_t pass2 = pass1 + 1; pass2 < num_passes; ++pass2) {
                                
                                // count up the sequence that belongs to only one of this pair of passes and the sequence
                                // that belongs to neither
                                size_t length_disjoint1 = 0, length_disjoint2 = 0, length_nonoverlapping = 0;
                                for (const auto& pass_set_record : path_pass_set_length.second) {
                                    if (pass_set_record.first[pass1] && !pass_set_record.first[pass2]) {
                                        length_disjoint1 += pass_set_record.second;
                                    }
                                    else if (!pass_set_record.first[pass1] && pass_set_record.first[pass2]) {
                                        length_disjoint2 += pass_set_record.second;
                                    }
                                    else if (!pass_set_record.first[pass1] && !pass_set_record.first[pass2]) {
                                        length_nonoverlapping += pass_set_record.second;
                                    }
                                }
                                
                                if (length_disjoint1 > min_inconsistency_disjoint_length &&
                                    length_disjoint2 > min_inconsistency_disjoint_length &&
                                    (length_disjoint1 + length_disjoint2) / 2 + length_nonoverlapping > min_inconsistency_total_length) {
                                    // this looks like it could have been the result of an indel that got inconsistent placement over a bond
                                    is_suspicious = true;
                                    break;
                                }
                            }
                        }
                    }
                    
                    if (is_suspicious) {
                        // record that we found a suspicious interval and move past it
                        chain_suspicious_intervals[feature.first].emplace_back(window_start, last_nontrivial);
                        window_start = last_nontrivial - 1;
                    }
                }
                
                // walk the left end up to the start of the next non-trivial snarl
                ++window_start;
                while (window_start < chain.size() && !nontrivial_left_boundary[snarls.structure_boundaries(chain[window_start]).first]) {
                    ++window_start;
                }
                if (window_start < chain.size()) {
                    left_path_positions = query_path_positions(snarls.structure_boundaries(chain[window_start]).first);
                }
                // it's possible that this will be beyond the position of the current end of the window, in which case we bump it over
                if (window_end < window_start) {
                    window_end = window_start;
                    if (window_end < chain.size()) {
                        right_path_positions = query_path_positions(snarls.structure_boundaries(chain[window_end]).second);
                    }
                }
                // update the dista
                if (window_end < chain.size()) {
                    dist = min_path_distance(left_path_positions, right_path_positions, multipass_intervals);
                }
            }
            
            uint64_t snarl_id = snarls.structure_containing(feature.first);
            if (snarl_id != -1) {
                auto& parent_list = snarl_cyclic_descendents[snarl_id];
                parent_list.splice(parent_list.end(), chain_cyclic_descendents[feature.first]);
            }
        }
    }
    
    // TODO: implement
    return std::vector<std::pair<uint64_t, uint64_t>>();
}


std::vector<std::pair<uint64_t, uint64_t>> InconsistencyIdentifier::identify_tight_cycles(const SnarlTree& snarls, const StepIndex step_index,
                                                                                          const std::vector<bool>& nontrivial_left_boundary,
                                                                                          const std::vector<bool>& nontrivial_right_boundary) const {
    
    static const bool debug = false;
    if (debug) {
        std::cerr << "identifying tight cycles with size limit " << max_tight_cycle_size << '\n';
    }
    
    // keep track of which chains/snarls are too big to be identified as a tight cycle
    std::vector<bool> chain_blocked(snarls.chain_size(), false);
    std::vector<bool> snarl_blocked(snarls.structure_size(), false);
    
    // partial lists of tight-cycle containing snarls that bubble up toward the root over the course
    // of iteration
    std::vector<std::list<uint64_t>> chain_cyclic_descendents(snarls.chain_size());
    std::vector<std::list<uint64_t>> snarl_cyclic_descendents(snarls.structure_size());
    
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
            
            if (!nontrivial_left_boundary[start_node] || !nontrivial_right_boundary[end_node]) {
                // this must be a trivial snarl, not worth spending compute to evaluate
                continue;
            }
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
