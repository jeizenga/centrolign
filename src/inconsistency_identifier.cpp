#include "centrolign/inconsistency_identifier.hpp"

#include <map>
#include <deque>
#include <functional>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <list>
#include <array>

#include "centrolign/utility.hpp"


namespace centrolign {

std::vector<std::pair<uint64_t, uint64_t>> InconsistencyIdentifier::identify_inconsistent_bonds(const SnarlTree& snarls, const StepIndex step_index,
                                                                                                const std::vector<bool>& nontrivial_left_boundary) const {
    
    static const bool debug = false;
    if (debug) {
        std::cerr << "looking for regions that have inconsistent indel placement across bonds\n";
    }
    
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
    
    // maximum distance over any pass of a specific path
    auto max_path_distance = [&](const std::unordered_map<uint64_t, std::vector<size_t>>& path_positions_left,
                                 const std::unordered_map<uint64_t, std::vector<size_t>>& path_positions_right,
                                 uint64_t path_id) -> size_t {
        size_t dist = 0;
        auto it1 = path_positions_left.find(path_id);
        auto it2 = path_positions_right.find(path_id);
        assert(it1->second.size() == it2->second.size());
        for (size_t i = 0; i < it1->second.size(); ++i) {
            dist = std::max(dist, it2->second[i] - it1->second[i] + 1);
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
                dists.push_back(it->second[i] - pos_record.second[i] + 1);
            }
        }
        return median(dists.begin(), dists.end());
    };
    
    std::vector<std::pair<uint64_t, uint64_t>> inconsistent_bonds;
    
    // we traverse in pre-order, starting with top-level chains
    std::deque<std::pair<uint64_t, bool>> queue;
    for (uint64_t chain_id = 0; chain_id < snarls.chain_size(); ++chain_id) {
        if (snarls.structure_containing(chain_id) == -1) {
            queue.emplace_back(chain_id, true);
        }
    }
    while (!queue.empty()) {
        auto feature = queue.front();
        queue.pop_front();
        if (!feature.second) {
            // snarl
            for (auto chain_id : snarls.chains_inside(feature.first)) {
                queue.emplace_back(chain_id, true);
            }
        }
        else {
            // chain
                        
            // look for a non-trivial snarl in this chain
            const auto& chain = snarls.structures_inside(feature.first);
            std::vector<size_t> nontrivial_snarls;
            for (size_t i = 0; i < chain.size(); ++i) {
                if (nontrivial_left_boundary[snarls.structure_boundaries(chain[i]).first]) {
                    nontrivial_snarls.push_back(i);
                }
            }
            
            if (nontrivial_snarls.empty()) {
                // there were no non-trivial snarls on this chain, we don't need to recurse into the children
                if (debug) {
                    std::cerr << "no non-trivial snarls\n";
                }
                continue;
            }
            
            if (debug) {
                std::cerr << "contains " << nontrivial_snarls.size() << " non-trivial snarls:\n";
                for (auto i : nontrivial_snarls) {
                    std::cerr << '\t' << chain[i] << ": " << snarls.structure_boundaries(chain[i]).first << ", " << snarls.structure_boundaries(chain[i]).second << '\n';
                }
            }
            
            // create a look-up structure for which pass of a path a position is on
            auto multipass_intervals = query_path_positions(snarls.structure_boundaries(chain.front()).first);
            for (auto it = multipass_intervals.begin(); it != multipass_intervals.end(); ) {
                auto prev = it++;
                // remove paths with a single interval
                if (prev->second.size() == 1) {
                    multipass_intervals.erase(prev);
                }
            }
            
            if (debug) {
                std::cerr << "found " << multipass_intervals.size() << " multi-pass intervals\n";
            }
            
            std::vector<bool> nontrivial_snarl_used(nontrivial_snarls.size(), false);
            
            if (!multipass_intervals.empty()) {
                // this chain is traversed multiple times
                
                
                
                
                // shrink the buckets, since we may have deleted many paths and we will iterate over this structure
                multipass_intervals.rehash(multipass_intervals.size());
                // add the right ends of the intervals
                for (const auto& pos_record : query_path_positions(snarls.structure_boundaries(chain.back()).second)) {
                    auto it = multipass_intervals.find(pos_record.first);
                    if (it != multipass_intervals.end()) {
                        for (auto pos : pos_record.second) {
                            it->second.push_back(pos);
                        }
                        // sort for fast lookup
                        std::sort(it->second.begin(), it->second.end());
                    }
                }
                if (debug) {
                    for (const auto& r : multipass_intervals) {
                        std::cerr << "path " << r.first << ":\n";
                        for (int i = 0; i < r.second.size(); i += 2) {
                            std::cerr << '\t' << r.second[i] << '\t' << r.second[i + 1] << '\n';
                        }
                    }
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
                
                // TODO: replace vector of bools with bitset?
                using pass_set_lengths_t = std::unordered_map<uint64_t, std::map<std::vector<bool>, size_t>>;
                
                auto merge_pass_set_lengths = [](pass_set_lengths_t& into, const pass_set_lengths_t& from) {
                    for (const auto& path_from : from) {
                        auto path_it = into.find(path_from.first);
                        if (path_it == into.end()) {
                            // this path has not been seen before, we simply rhe set lengths over
                            into[path_from.first] = path_from.second;
                        }
                        else {
                            // TODO: it would be cool if there were a find_hint so i could read-modify-write instead of
                            // just replacing
                            auto& path_pass_set_length = path_it->second;
                            for (const auto& pass_set_length : path_from.second) {
                                path_pass_set_length[pass_set_length.first] += pass_set_length.second;
                            }
                        }
                    }
                };
                
                std::vector<pass_set_lengths_t> snarl_pass_set_lengths(nontrivial_snarls.size());
                
                for (size_t i = 0; i < nontrivial_snarls.size(); ++i) {
                    auto& pass_set_length = snarl_pass_set_lengths[i];
                    for (uint64_t chain_id : snarls.chains_inside(chain[nontrivial_snarls[i]])) {
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
                            if (debug) {
                                std::cerr << "snarl " << i << ", add chain of length " << length << " to pass coverage combo ";
                                for (auto w : which_passes) {
                                    std::cerr << w;
                                }
                                std::cerr << '\n';
                            }
                            pass_set_length[pass_record.first][std::move(which_passes)] += length;
                        }
                    }
                }
                
                std::vector<std::pair<size_t, pass_set_lengths_t>> window_path_set_lengths(nontrivial_snarls.size());
                for (size_t i = 0; i < window_path_set_lengths.size(); ++i) {
                    window_path_set_lengths[i] = std::make_pair(i, snarl_pass_set_lengths[i]);
                }
                size_t window_steps = 1;
                
                while (!window_path_set_lengths.empty()) {
                    
                    // the queue for the next higher window size
                    decltype(window_path_set_lengths) next_window_path_set_lengths;
                    
                    // note: iterate backwards so that we will always know if the next window to the right is already used before expanding into it
                    for (size_t i = window_path_set_lengths.size() - 1; i < window_path_set_lengths.size(); --i) {
                        
                        auto window = std::move(window_path_set_lengths[i]);
                        auto& pass_set_length = window.second;
                        
                        if (debug) {
                            std::cerr << "checking window containing " << window_steps << " non-trivial snarls from " << chain[nontrivial_snarls[window.first]] << " to " << chain[nontrivial_snarls[window.first + window_steps - 1]] << " encompassing nodes " << snarls.structure_boundaries(chain[nontrivial_snarls[window.first]]).first << " to " << snarls.structure_boundaries(chain[nontrivial_snarls[window.first + window_steps - 1]]).second << '\n';
                        }
                        
                        bool is_suspicious = false;
                        for (auto it = pass_set_length.begin(), end = pass_set_length.end(); it != end && !is_suspicious; ++it) {
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
                                    
                                    if (length_disjoint1 >= min_inconsistency_disjoint_length &&
                                        length_disjoint2 >= min_inconsistency_disjoint_length &&
                                        (length_disjoint1 + length_disjoint2) / 2 + length_nonoverlapping >= min_inconsistency_total_length) {
                                        // this looks like it could have been the result of an indel that got inconsistent placement over a bond
                                        if (debug) {
                                            std::cerr << "window looks suspicious with pass combination " << pass1 << " and " << pass2 << ": dj1 " << length_disjoint1 << ", dj2 " << length_disjoint2 << ", nol " << length_nonoverlapping << '\n';
                                        }
                                        is_suspicious = true;
                                        break;
                                    }
                                }
                            }
                        }
                        
                        if (is_suspicious) {
                            // looks suspicious, add it to the output and mark it used
                            inconsistent_bonds.emplace_back(snarls.structure_boundaries(chain[nontrivial_snarls[window.first]]).first,
                                                            snarls.structure_boundaries(chain[nontrivial_snarls[window.first + window_steps - 1]]).second);
                            for (size_t j = window.first, n = window.first + window_steps; j < n; ++j) {
                                nontrivial_snarl_used[j] = true;
                            }
                        }
                        else if (window.first + window_steps < nontrivial_snarls.size() &&
                                 !nontrivial_snarl_used[window.first + window_steps] &&
                                 nontrivial_snarls[window.first + window_steps] - nontrivial_snarls[window.first] < max_bond_inconsistency_window) {
                            if (debug) {
                                std::cerr << "queueing window of size " << (window_steps + 1) << " non-trivial snarls starting at " << window.first << "-th non-trivial snarl\n";
                            }
                            // expand the window and add it to the next queue
                            merge_pass_set_lengths(window.second, snarl_pass_set_lengths[window.first + window_steps]);
                            next_window_path_set_lengths.emplace_back(std::move(window));
                        }
                        
                    }
                    
                    ++window_steps;
                    window_path_set_lengths = std::move(next_window_path_set_lengths);
                }
            }
            
            // continue searching in the chains on any snarls that we didn't use
            for (size_t i = 0; i < nontrivial_snarls.size(); ++i) {
                if (!nontrivial_snarl_used[i]) {
                    queue.emplace_back(chain[nontrivial_snarls[i]], false);
                }
            }
        }
    }
    
    return inconsistent_bonds;
}


std::vector<std::pair<uint64_t, uint64_t>> InconsistencyIdentifier::identify_tight_cycles(const SnarlTree& snarls, const StepIndex step_index,
                                                                                          const std::vector<bool>& nontrivial_left_boundary) const {
    
    static const bool debug = false;
    if (debug) {
        std::cerr << "identifying tight cycles with size limit " << max_tight_cycle_size << '\n';
    }
    
    // FIXME: how to handle adjacent snarls that are both flagged??
    
    // keep track of which chains/snarls are too big to be identified as a tight cycle
    std::vector<bool> chain_blocked(snarls.chain_size(), false);
    std::vector<bool> snarl_blocked(snarls.structure_size(), false);
    
    // partial lists of tight-cycle containing snarls that bubble up toward the root over the course
    // of iteration
    std::vector<std::list<uint64_t>> chain_cyclic_descendents(snarls.chain_size());
    std::vector<std::list<uint64_t>> snarl_cyclic_descendents(snarls.structure_size());
    
    for (auto feature : snarls.postorder()) {
        
        uint64_t start_node, end_node;
        if (feature.second) {
            // is a chain
            if (debug) {
                std::cerr << "at chain " << feature.first << " between " << snarls.structure_boundaries(snarls.structures_inside(feature.first).front()).first << " and " << snarls.structure_boundaries(snarls.structures_inside(feature.first).back()).second << '\n';
            }
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
            if (debug) {
                std::cerr << "at snarl " << feature.first << " with boundaries " << snarls.structure_boundaries(feature.first).first << ' ' << snarls.structure_boundaries(feature.first).second << '\n';
            }
            if (snarl_blocked[feature.first]) {
                // already too big in a descendent
                chain_blocked[snarls.chain_containing(feature.first)] = true;
                if (debug) {
                    std::cerr << "blocked\n";
                }
                continue;
            }
            std::tie(start_node, end_node) = snarls.structure_boundaries(feature.first);
            
            if (!nontrivial_left_boundary[start_node]) {
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
                chain_blocked[snarls.chain_containing(feature.first)] = true;
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
