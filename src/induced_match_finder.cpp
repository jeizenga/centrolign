#include "centrolign/induced_match_finder.hpp"

#include <cassert>
#include <limits>
#include <cmath>
#include <cstdlib>


namespace centrolign {


InducedMatchFinder::InducedMatchFinder(const BaseGraph& full_graph, const std::vector<match_set_t>& matches,
                                       const std::vector<std::pair<uint64_t, uint64_t>>& components,
                                       const StepIndex& step_index) : parent(&full_graph), component_path_hits(components.size()) {
    
    // label nodes as coming from one of the 2-disconnected components with DFS
    std::vector<size_t> node_to_component(full_graph.node_size(), -1);
    for (size_t i = 0; i < components.size(); ++i) {
        const auto& component = components[i];
        std::vector<uint64_t> stack(1, component.first);
        node_to_component[component.first] = i;
        node_to_component[component.second] = i;
        while (!stack.empty()) {
            auto node_id = stack.back();
            stack.pop_back();
            for (auto next_id : full_graph.next(node_id)) {
                if (node_to_component[next_id] == -1) {
                    node_to_component[next_id] = i;
                    stack.push_back(next_id);
                }
            }
        }
    }
        
    for (size_t i = 0; i < matches.size(); ++i) {
        
        const auto& match_set = matches[i];
        
        // to keep track of whether we've initialize a hit set for this match set on each component
        std::unordered_set<size_t> component_hits_initialized;
        
        for (size_t j = 0; j < match_set.walks1.size(); ++j) {
            
            std::unordered_set<size_t> overlapping_components;
            
            const auto& walk = match_set.walks1[j];
            for (size_t k = 0; k < walk.size(); ++k) {
                auto component = node_to_component[walk[k]] ;
                if (component != -1) {
                    overlapping_components.insert(component);
                }
            }
            if (overlapping_components.empty()) {
                // we don't need to worry about this match set, it doesn't touch any of
                // the components we're inducing matches on
                continue;
            }
            for (size_t comp : overlapping_components) {
                if (component_hits_initialized.insert(comp).second) {
                    // we still need to initialize the component hit set for this component
                    component_path_hits[comp].emplace_back();
                    auto& hit_set = component_path_hits[comp].back();
                    hit_set.length = walk.size();
                    hit_set.deduplicated_count = match_set.walks1.size();
                }
            }
            
            // find paths that contain this walk as a complete subpath
            std::unordered_set<std::pair<uint64_t, size_t>> extensions;
            for (const auto& step : step_index.path_steps(walk.front())) {
                extensions.insert(step);
            }
            for (size_t k = 1; k < walk.size() && !extensions.empty(); ++k) {
                std::unordered_set<std::pair<uint64_t, size_t>> next_extensions;
                
                for (const auto& step : step_index.path_steps(walk[k])) {
                    if (extensions.count(std::make_pair(step.first, step.second - 1))) {
                        next_extensions.insert(step);
                    }
                }
                extensions = std::move(next_extensions);
            }
            
            // add back the original start position for each walk on each path
            for (const auto& final_extension : extensions) {
                for (size_t comp : overlapping_components) {
                    component_path_hits[comp].back().hit_locations[final_extension.first].emplace_back(final_extension.second + 1 - walk.size(), j);
                }
            }
        }
        
        // check if we should keep this hit set
        for (size_t comp : component_hits_initialized) {
            auto& hit_locations = component_path_hits[comp].back().hit_locations;
            if (hit_locations.empty() || (hit_locations.size() == 1 && hit_locations.begin()->second.size() == 1)) {
                // we can never make a match across two subpaths or two different paths, so get rid of this set
                // TODO: this doesn't catch all useless hit sets, but doing that requires knowledge of the tree
                component_path_hits[comp].pop_back();
            }
            else {
                // make the locations be increasing in path offset
                for (auto it = hit_locations.begin(); it != hit_locations.end(); ++it) {
                    std::sort(it->second.begin(), it->second.end());
                }
            }
        }
    }
}

InducedMatchFinderComponentView InducedMatchFinder::component_view(size_t comp) const {
    
    return InducedMatchFinderComponentView(parent, &component_path_hits[comp]);
}

std::tuple<std::string, size_t, size_t> InducedMatchFinderComponentView::parse_subpath_name(const std::string& subpath_name) const {
    size_t sep = subpath_name.find_last_of(':');
    auto it = std::find(subpath_name.begin() + sep + 1, subpath_name.end(), '-');
    return std::tuple<std::string, size_t, size_t>(subpath_name.substr(0, sep),
                                                   parse_int(std::string(subpath_name.begin() + sep + 1, it)),
                                                   parse_int(std::string(it + 1, subpath_name.end())));
    
}


}
