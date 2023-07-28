#include "centrolign/superbubble_distance_oracle.hpp"

#include <unordered_set>

#include "centrolign/utility.hpp"

namespace centrolign {

using namespace std;


vector<pair<uint64_t, bool>> SuperbubbleDistanceOracle::path_to_root(uint64_t node_id) const {
    
    vector<pair<uint64_t, bool>> path;
    
    path.emplace_back(node_to_superbubble[node_id], false);
    while (path.back() != pair<uint64_t, bool>(superbubble_tree.superbubble_size(), false)) {
        if (path.back().second) {
            // this is a chain
            auto bub_id = superbubble_tree.superbubble_containing(path.back().first);
            if (bub_id == -1) {
                // the outside bubble net graph
                bub_id = superbubble_tree.superbubble_size();
            }
            path.emplace_back(bub_id, false);
        }
        else {
            // all non-exterior bubbles are in a chain
            path.emplace_back(superbubble_tree.chain_containing(path.back().first), true);
        }
    }
    
    return path;
}

size_t SuperbubbleDistanceOracle::min_distance(uint64_t node_id1, uint64_t node_id2) const {
    
    auto path1 = path_to_root(node_id1);
    auto path2 = path_to_root(node_id1);
    
    // find the lowest common ancestor in the snarl tree
    unordered_set<pair<uint64_t, bool>> path1_steps;
    for (const auto& step : path1) {
        path1_steps.insert(step);
    }
    size_t idx2 = 0;
    while (!path1_steps.count(path2[idx2])) {
        ++idx2;
    }
    size_t idx1 = 0;
    while (path1[idx1] != path2[idx2]) {
        ++idx1;
    }
    
    // figure out the distance within the shared component
    size_t dist;
    if (path1[idx1].second) {
        // the lowest shared feature is a chain
        auto chain_idx1 = superbubble_link_index[path1[idx1 - 1].first];
        auto chain_idx2 = superbubble_link_index[path2[idx2 - 1].first];
        assert(chain_idx1 != chain_idx2); // or else the common ancestor is a snarl
        if (chain_idx1 > chain_idx2) {
            // node 1 is further along the chain than node2
            return -1;
        }
        const auto& prefix_sum = chain_prefix_sums[path1[idx1].first];
        dist = prefix_sum.first[chain_idx2] - prefix_sum.first[chain_idx1 + 1];
    }
    else {
        // the lowest shared feature is a superbubble
        pair<uint64_t, bool> feature_id1, feature_id2;
        if (idx1 == 0) {
            feature_id1 = make_pair(node_id1, false);
        }
        else {
            feature_id1 = path1[idx1 - 1];
        }
        if (idx2 == 0) {
            feature_id2 = make_pair(node_id2, false);
        }
        else {
            feature_id2 = path2[idx2 - 1];
        }
        const auto& table = net_graph_tables[path1[idx1].first];
        auto it = table.find(make_pair(feature_id1, feature_id2));
        if (it == table.end()) {
            return -1;
        }
        else {
            dist = it->second;
        }
    }
    
    // walk to the right ends in the first path
    for (size_t i = 0; i < idx1; ++i) {
        if (path1[i].second) {
            // chain
            const auto& prefix_sum = chain_prefix_sums[path1[i].first].first;
            size_t link_idx = superbubble_link_index[path1[i - 1].first];
            dist += prefix_sum.back() - prefix_sum[link_idx + 1];
        }
        else {
            // superbubble
            auto bub_id = path1[i].first;
            const auto& table = net_graph_tables[bub_id];
            auto sink_label = make_pair(superbubble_tree.superbubble_boundaries(bub_id).second, false);
            if (i == 0) {
                // the subfeature is the original node
                dist += table.at(make_pair(make_pair(node_id1, false), sink_label));
            }
            else {
                // the subfeature is a chain, but we are walking from its right end,
                // so we have to adjust the recorded length
                dist += (table.at(make_pair(path1[i - 1], sink_label))
                         - chain_prefix_sums[path1[i - 1].first].first.back());
            }
        }
    }
    
    // walk to the left ends in the second path
    for (size_t i = 0; i < idx2; ++i) {
        if (path2[i].second) {
            // chain
            const auto& prefix_sum = chain_prefix_sums[path2[i].first].first;
            size_t link_idx = superbubble_link_index[path2[i - 1].first];
            dist += prefix_sum[link_idx];
        }
        else {
            // superbubble
            auto bub_id = path2[i].first;
            const auto& table = net_graph_tables[bub_id];
            auto source_label = make_pair(superbubble_tree.superbubble_boundaries(bub_id).first, false);
            if (i == 0) {
                // the subfeature is the original node
                dist += table.at(make_pair(source_label, make_pair(node_id2, false)));
            }
            else {
                // the subfeature is a chain
                dist += table.at(make_pair(source_label, path2[i - 1]));
            }
        }
    }
    
    return dist;
}

}
