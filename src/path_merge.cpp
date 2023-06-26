#include "centrolign/path_merge.hpp"

namespace centrolign {

using namespace std;


PathMerge::PathMerge(const BaseGraph& graph) : PathMerge(graph, nullptr) {
    
}

PathMerge::PathMerge(const BaseGraph& graph, const SentinelTableau& tableau) : PathMerge(graph, &tableau) {
    src_id = tableau.src_id;
    snk_id = tableau.snk_id;
}

PathMerge::PathMerge(const BaseGraph& graph, const SentinelTableau* tableau) :
    path_head(graph.node_size(), -1),
    index_on_path(graph.node_size(),
                  std::vector<std::pair<size_t, uint64_t>>(graph.path_size() + int(tableau != nullptr),
                                                           std::pair<size_t, uint64_t>(-1, -1))),
    table(graph.node_size(), std::vector<size_t>(graph.path_size() + int(tableau != nullptr), -1)),
    g(&graph)
{
    
    static const bool debug = false;
    
    // identify steps of paths and seed the DP
    for (uint64_t path_id = 0; path_id < graph.path_size(); ++path_id) {
        size_t index = 0;
        for (auto node_id : graph.path(path_id)) {
            for (auto next_id : graph.next(node_id)) {
                // note: we don't need the max here because we're going in increasing order
                table[next_id][path_id] = index;
            }
            index_on_path[node_id][path_id] = std::make_pair(index, path_head[node_id]);
            path_head[node_id] = path_id;
            ++index;
        }
    }
    
    // use DP to fill it in between cross-path edges
    for (auto node_id : topological_order(graph)) {
        auto& row = table[node_id];
        for (auto prev_id : graph.previous(node_id)) {
            const auto& prev_row = table[prev_id];
            for (uint64_t path_id = 0; path_id < graph.path_size(); ++path_id) {
                // cast to int64_t so that -1's get overwritten
                row[path_id] = std::max<int64_t>(prev_row[path_id], row[path_id]);
            }
        }
    }
    
    // add a "pseudo-path" chain for the sentinels
    if (tableau) {
        index_on_path[tableau->src_id][graph.path_size()].first = 0;
        index_on_path[tableau->snk_id][graph.path_size()].first = 1;
        path_head[tableau->src_id] = graph.path_size();
        path_head[tableau->snk_id] = graph.path_size();
        for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
            if (node_id != tableau->src_id) {
                table[node_id][graph.path_size()] = 0;
            }
        }
    }
    
    if (debug) {
        std::cerr << "path lists:\n";
        for (uint64_t n = 0; n < graph.node_size(); ++n) {
            std::cerr << n << ":\t" << path_head[n] << '\n';
        }
        
        std::cerr << "index on paths:\n";
        for (uint64_t n = 0; n < graph.node_size(); ++n) {
            std::cerr << n << ":";
            for (uint64_t p = 0; p < index_on_path[n].size(); ++p) {
                std::cerr << '\t' << (int64_t) index_on_path[n][p].first << ',' << (int64_t) index_on_path[n][p].second;
            }
            std::cerr << '\n';
        }
        
        std::cerr << "table:\n";
        for (uint64_t n = 0; n < graph.node_size(); ++n) {
            std::cerr << n << ":";
            for (uint64_t p = 0; p < table[n].size(); ++p) {
                std::cerr << '\t' << (int64_t) table[n][p];
            }
            std::cerr << '\n';
        }
    }
}

vector<vector<pair<uint64_t, uint64_t>>> PathMerge::chain_forward_edges() const {
    
    // identify the forward links
    vector<vector<pair<uint64_t, uint64_t>>> forward(table.size());
    for (uint64_t node_id = 0; node_id < table.size(); ++node_id) {
        const auto& row = table[node_id];
        for (uint64_t p = 0; p < g->path_size(); ++p) {
            if (row[p] != -1) {
                forward[g->path(p)[row[p]]].emplace_back(node_id, p);
            }
        }
        // TODO: do i really even need forward links on the pseudo-path?
        // i could get rid of the saved src_id if i didn't make these
        if (src_id != -1 && node_id != src_id) {
            forward[src_id].emplace_back(node_id, g->path_size());
        }
    }
    
    return forward;
}

}
