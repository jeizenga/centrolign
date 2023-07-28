#include "centrolign/superbubble_distance_oracle.hpp"

#include <unordered_set>

#include "centrolign/utility.hpp"

namespace centrolign {

using namespace std;

const bool SuperbubbleDistanceOracle::debug = false;

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
    auto path2 = path_to_root(node_id2);
    
    if (debug) {
        std::cerr << "measuring distance from " << node_id1 << " -> " << node_id2 << '\n';
        std::cerr << "paths to root:\n";
        for (auto path : {path1, path2}) {
            for (auto s : path) {
                std::cerr << ' ' << s.first << ',' << s.second;
            }
            std::cerr << '\n';
        }
    }
    
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
    
    if (debug) {
        std::cerr << "LCA occurs in " << path1[idx1].first << "," << path1[idx1].second << " at steps " << idx1 << " and " << idx2 << '\n';
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
            if (debug) {
                std::cerr << "unreachable along chain\n";
            }
            return -1;
        }
        const auto& prefix_sum = chain_prefix_sums[path1[idx1].first];
        dist = prefix_sum[chain_idx2] - prefix_sum[chain_idx1 + 1];
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
            if (debug) {
                std::cerr << "unreachable within net graph\n";
            }
            return -1;
        }
        else if (feature_id1.second) {
            // subfeature is a chain, but we are leaving out its right side
            dist = it->second - chain_prefix_sums[feature_id1.first].back();
        }
        else {
            dist = it->second;
        }
    }
    
    if (debug) {
        std::cerr << "init with distance " << dist << " within shared component\n";
    }
    
    // walk to the right ends in the first path
    for (size_t i = 0; i < idx1; ++i) {
        if (path1[i].second) {
            // chain
            const auto& prefix_sum = chain_prefix_sums[path1[i].first];
            size_t link_idx = superbubble_link_index[path1[i - 1].first];
            dist += prefix_sum.back() - prefix_sum[link_idx + 1];
            if (debug) {
                std::cerr << "increment to " << dist << " from path 1 chain " << path1[i].first << "\n";
            }
        }
        else {
            // superbubble
            auto bub_id = path1[i].first;
            const auto& table = net_graph_tables[bub_id];
            auto sink_label = make_pair(superbubble_tree.superbubble_boundaries(bub_id).second, false);
            if (i == 0) {
                // the subfeature is the original node
                dist += table.at(make_pair(make_pair(node_id1, false), sink_label));
                if (debug) {
                    std::cerr << "increment to " << dist << " from path 1 base bubble " << bub_id << "\n";
                }
            }
            else {
                // the subfeature is a chain, but we are walking from its right end,
                // so we have to adjust the recorded length
                dist += (table.at(make_pair(path1[i - 1], sink_label))
                         - chain_prefix_sums[path1[i - 1].first].back());
                if (debug) {
                    std::cerr << "increment to " << dist << " from path 1 containing bubble " << bub_id << "\n";
                }
            }
        }
    }
    
    // walk to the left ends in the second path
    for (size_t i = 0; i < idx2; ++i) {
        if (path2[i].second) {
            // chain
            const auto& prefix_sum = chain_prefix_sums[path2[i].first];
            size_t link_idx = superbubble_link_index[path2[i - 1].first];
            dist += prefix_sum[link_idx];
            if (debug) {
                std::cerr << "increment to " << dist << " from path 2 chain " << path2[i].first << "\n";
            }
        }
        else {
            // superbubble
            auto bub_id = path2[i].first;
            const auto& table = net_graph_tables[bub_id];
            auto source_label = make_pair(superbubble_tree.superbubble_boundaries(bub_id).first, false);
            if (i == 0) {
                // the subfeature is the original node
                dist += table.at(make_pair(source_label, make_pair(node_id2, false)));
                if (debug) {
                    std::cerr << "increment to " << dist << " from path 2 base bubble " << bub_id << "\n";
                }
            }
            else {
                // the subfeature is a chain
                dist += table.at(make_pair(source_label, path2[i - 1]));
                if (debug) {
                    std::cerr << "increment to " << dist << " from path 2 containing bubble " << bub_id << "\n";
                }
            }
        }
    }
    
    return dist;
}

}
