#include "centrolign/test_util.hpp"

#include <set>
#include <map>

namespace centrolign {

using namespace std;

// DFS
bool is_reachable(const BaseGraph& graph, uint64_t id_from, uint64_t id_to) {
    
    vector<bool> traversed(graph.node_size(), false);
    vector<uint64_t> stack(1, id_from);
    traversed[id_from] = true;
    while (!stack.empty()) {
        auto here = stack.back();
        stack.pop_back();
        for (auto nid : graph.next(here)) {
            if (nid == id_to) {
                return true;
            }
            if (!traversed[nid]) {
                traversed[nid] = true;
                stack.push_back(nid);
            }
        }
    }
    return false;
}

vector<vector<uint64_t>> all_paths_internal(const BaseGraph& graph,
                                            uint64_t id_from, uint64_t id_to) {
    
    vector<vector<uint64_t>> paths;
    
    vector<tuple<uint64_t, vector<uint64_t>, size_t>> stack;
    stack.emplace_back(id_from, graph.next(id_from), 0);
    
    while (!stack.empty()) {
        auto& top = stack.back();
        if (get<0>(top) == id_to ||
            (id_to == -1 && graph.next_size(get<0>(top)) == 0)) {
            paths.emplace_back();
            for (auto& r : stack) {
                paths.back().push_back(get<0>(r));
            }
            stack.pop_back();
            continue;
        }
        if (get<2>(top) == get<1>(top).size()) {
            stack.pop_back();
            continue;
        }
        auto next = get<1>(top)[get<2>(top)++];
        stack.emplace_back(next, graph.next(next), 0);
    }
    
    return paths;
}

vector<vector<uint64_t>> all_paths(const BaseGraph& graph, uint64_t id_from) {
    return all_paths_internal(graph, id_from, -1);
}

vector<vector<uint64_t>> all_paths(const BaseGraph& graph,
                                   uint64_t id_from, uint64_t id_to) {
    return all_paths_internal(graph, id_from, id_to);
}

bool graphs_are_equivalent(const BaseGraph& graph1, const BaseGraph& graph2) {
    std::vector<uint64_t> trans1, trans2;
    for (uint64_t n = 0; n < graph1.node_size(); ++n) {
        trans1.push_back(n);
    }
    for (uint64_t n = 0; n < graph2.node_size(); ++n) {
        trans2.push_back(n);
    }
    return translated_graphs_are_identical(graph1, graph2, trans1, trans2);
}

bool translated_graphs_are_identical(const BaseGraph& subgraph1, const BaseGraph& subgraph2,
                                     const std::vector<uint64_t>& back_translation1,
                                     const std::vector<uint64_t>& back_translation2) {
    // same number of nodes
    if (subgraph1.node_size() != subgraph2.node_size()) {
        return false;
    }
    
    size_t num_edges1 = 0, num_edges2 = 0;
    unordered_map<uint64_t, uint64_t> fwd_translation1, fwd_translation2;
    for (uint64_t n = 0; n < subgraph1.node_size(); ++n) {
        fwd_translation1[back_translation1[n]] = n;
        fwd_translation2[back_translation2[n]] = n;
        num_edges1 += subgraph1.next_size(n);
        num_edges2 += subgraph2.next_size(n);
    }
    
    if (num_edges1 != num_edges2) {
        return false;
    }
    
    // exact same nodes
    for (auto r : fwd_translation1) {
        if (!fwd_translation2.count(r.first)) {
            return false;
        }
        if (subgraph1.label(r.second) != subgraph2.label(fwd_translation2[r.first])) {
            return false;
        }
    }
    
    set<pair<uint64_t, uint64_t>> edges1, edges2;
    for (uint64_t n = 0; n < subgraph1.node_size(); ++n) {
        for (auto i : subgraph1.next(n)) {
            edges1.emplace(back_translation1[n], back_translation1[i]);
        }
        for (auto i : subgraph2.next(n)) {
            edges2.emplace(back_translation2[n], back_translation2[i]);
        }
    }
    if (edges1 != edges2) {
        return false;
    }
    
    return true;
}

bool translations_possibly_consistent(const BaseGraph& subgraph1, const BaseGraph& subgraph2,
                                      const std::vector<uint64_t>& back_translation1,
                                      const std::vector<uint64_t>& back_translation2) {

    if (back_translation1.size() != back_translation2.size()) {
        return false;
    }
    
    map<uint64_t, set<uint64_t>> fwd_trans1, fwd_trans2;
    
    for (uint64_t i = 0; i < back_translation1.size(); ++i) {
        fwd_trans1[back_translation1[i]].insert(i);
        fwd_trans2[back_translation2[i]].insert(i);
    }
    
    auto dump_fwd = [&]() {
        int i = 1;
        for (auto fwd_trans : {fwd_trans1, fwd_trans2}) {
            cerr << "forward translation " << i++ << "\n";
            for (const auto& r : fwd_trans) {
                cerr << r.first << ':';
                for (auto v : r.second) {
                    cerr << ' ' << v;
                }
                cerr << '\n';
            }
        }
    };
    
    if (fwd_trans1.size() != fwd_trans2.size()) {
        cerr << "forward translations are not the same size\n";
        dump_fwd();
        return false;
    }
    
    for (const auto& r : fwd_trans1) {
        if (!fwd_trans2.count(r.first)) {
            cerr << "missing forward translation " << r.first << "\n";
            dump_fwd();
            return false;
        }
        const auto& s1 = r.second;
        const auto& s2 = fwd_trans2[r.first];
        
        if (s1.size() != s2.size()) {
            cerr << "node forward translations are not the same size on " << r.first << "\n";
            dump_fwd();
            return false;
        }
        
        for (auto v : s1) {
            if (subgraph1.label(v) != subgraph1.label(*s1.begin())) {
                cerr << "labels do not match internally for forward trans of " << r.first << "\n";
                dump_fwd();
                return false;
            }
        }
        for (auto v : s2) {
            if (subgraph2.label(v) != subgraph1.label(*s1.begin())) {
                cerr << "labels do not match for forward trans of " << r.first << "\n";
                dump_fwd();
                return false;
            }
        }
    }
    
    multiset<pair<uint64_t, uint64_t>> back_edges1, back_edges2;
    
    for (uint64_t n = 0; n < subgraph1.node_size(); ++n) {
        for (auto m : subgraph1.next(n)) {
            back_edges1.emplace(back_translation1[n], back_translation1[m]);
        }
    }
    for (uint64_t n = 0; n < subgraph2.node_size(); ++n) {
        for (auto m : subgraph2.next(n)) {
            back_edges2.emplace(back_translation2[n], back_translation2[m]);
        }
    }
    
    if (back_edges1 != back_edges2) {
        cerr << "translated edges do not match\n";
        return false;
    }
    
    return true;
}

bool possibly_isomorphic(const BaseGraph& graph1,
                         const BaseGraph& graph2) {
    
    if (graph1.node_size() != graph2.node_size()) {
        cerr << "graphs do not have same number of nodes: " << graph1.node_size() << " and " << graph2.node_size() << "\n";
        return false;
    }
    
    set<char> labels;
    for (uint64_t node_id = 0; node_id < graph1.node_size(); ++node_id) {
        labels.insert(graph1.label(node_id));
        labels.insert(graph2.label(node_id));
    }
    
    // divvy up all the analyses by label to increase specificity
    for (auto c : labels) {
        
        // degrees
        std::vector<size_t> in_degrees1, in_degrees2, out_degrees1, out_degrees2;
        // unique neighbor sets
        std::set<std::set<uint64_t>> in_nbrs1, in_nbrs2, out_nbrs1, out_nbrs2;
        // characters of neighbors
        std::vector<std::multiset<char>> in_chars1, in_chars2, out_chars1, out_chars2;
        
        for (uint64_t node_id = 0; node_id < graph1.node_size(); ++node_id) {
            if (graph1.label(node_id) == c) {
                // we match the label on graph 1
                std::set<uint64_t> p1, n1;
                std::multiset<char> pc1, nc1;
                in_degrees1.push_back(graph1.previous_size(node_id));
                out_degrees1.push_back(graph1.next_size(node_id));
                for (auto n : graph1.previous(node_id)) {
                    p1.insert(n);
                    pc1.insert(graph1.label(n));
                }
                for (auto n : graph1.next(node_id)) {
                    n1.insert(n);
                    nc1.insert(graph1.label(n));
                }
                in_nbrs1.insert(p1);
                out_nbrs1.insert(n1);
                in_chars1.push_back(pc1);
                out_chars1.push_back(nc1);
            }
            if (graph2.label(node_id) == c) {
                // we match the label on graph 2
                std::set<uint64_t> p2, n2;
                std::multiset<char> pc2, nc2;
                in_degrees2.push_back(graph2.previous_size(node_id));
                out_degrees2.push_back(graph2.next_size(node_id));
                for (auto n : graph2.previous(node_id)) {
                    p2.insert(n);
                    pc2.insert(graph2.label(n));
                }
                for (auto n : graph2.next(node_id)) {
                    n2.insert(n);
                    nc2.insert(graph2.label(n));
                }
                in_nbrs2.insert(p2);
                out_nbrs2.insert(n2);
                in_chars2.push_back(pc2);
                out_chars2.push_back(nc2);
            }
        }
        
        sort(in_degrees1.begin(), in_degrees1.end());
        sort(in_degrees2.begin(), in_degrees2.end());
        sort(out_degrees1.begin(), out_degrees1.end());
        sort(out_degrees2.begin(), out_degrees2.end());
        sort(in_chars1.begin(), in_chars1.end());
        sort(out_chars1.begin(), out_chars1.end());
        sort(in_chars2.begin(), in_chars2.end());
        sort(out_chars2.begin(), out_chars2.end());
        
        if (in_degrees1 != in_degrees2 ) {
            // in degree distributions are not identical
            cerr << "in degree distributions are not identical for char " << c << '\n';
            for (auto d : in_degrees1) {
                cerr << ' ' << d;
            }
            cerr << '\n';
            for (auto d : in_degrees2) {
                cerr << ' ' << d;
            }
            cerr << '\n';
            return false;
        }
        
        if (out_degrees1 != out_degrees2) {
            // out degree distributions are not identical
            cerr << "out degree distributions are not identical for char " << c << '\n';
            for (auto d : out_degrees1) {
                cerr << ' ' << d;
            }
            cerr << '\n';
            for (auto d : out_degrees2) {
                cerr << ' ' << d;
            }
            cerr << '\n';
            return false;
        }
        
        if (in_chars1 != in_chars2 || out_chars1 != out_chars2) {
            // the sets of label neighborhoods are not identical
            cerr << "label neighborhoods are not identical for char " << c << '\n';
            return false;
        }
        
        if (in_nbrs1.size() != in_nbrs2.size() || out_nbrs1.size() != out_nbrs2.size()) {
            // number of unique neighborhoods is not identical
            cerr << "unique neighbor sets are not identical for char " << c << '\n';
            return false;
        }
        
    }
    
    auto get_neighborhood = [](const BaseGraph& graph, uint64_t node_id) {
        unordered_map<uint64_t, uint64_t> fwd_trans;
        BaseGraph nbd;
        uint64_t center = nbd.add_node(graph.label(node_id));
        fwd_trans[node_id] = center;
        for (auto n : graph.previous(node_id)) {
            auto prev_id = nbd.add_node(graph.label(n));
            fwd_trans[n] = prev_id;
            nbd.add_edge(prev_id, center);
        }
        for (auto n : graph.next(node_id)) {
            auto next_id = nbd.add_node(graph.label(n));
            fwd_trans[n] = next_id;
            nbd.add_edge(center, next_id);
        }
        for (auto r : fwd_trans) {
            if (r.first == node_id) {
                continue;
            }
            for (auto n : graph.previous(r.first)) {
                if (fwd_trans.count(n) && n != node_id) {
                    nbd.add_edge(fwd_trans.at(n), r.second);
                }
            }
            for (auto n : graph.next(r.first)) {
                if (fwd_trans.count(n) && n != node_id) {
                    nbd.add_edge(r.second, fwd_trans.at(n));
                }
            }
        }
        return nbd;
    };
    
    if (graph1.path_size() != graph2.path_size()) {
        cerr << "differing path sizes: " << graph1.path_size()  << " " << graph2.path_size() << '\n';
        return false;
    }
    for (uint64_t path_id1 = 0; path_id1 < graph1.path_size(); ++path_id1) {
        auto path_id2 = graph2.path_id(graph1.path_name(path_id1));
        if (path_id2 == -1) {
            cerr << "path " << graph1.path_name(path_id1) << " is not found in graph 2\n";
            return false;
        }
        auto path1 = graph1.path(path_id1);
        auto path2 = graph2.path(path_id2);
        if (path1.size() != path2.size()) {
            cerr << "path " << graph1.path_name(path_id1) << " is not the same size\n";
            return false;
        }
        for (size_t i = 0; i < path1.size(); ++i) {
            BaseGraph nbd1 = get_neighborhood(graph1, path1[i]);
            BaseGraph nbd2 = get_neighborhood(graph2, path2[i]);
            if (!possibly_isomorphic(nbd1, nbd2)) {
                cerr << "not isomorphic in neighborhood of step " << i << " on path " << graph1.path_name(path_id1) << '\n';
                return false;
            }
        }
    }
    
    return true;
}

void check_alignment(const Alignment& got, const Alignment& expected) {
    if (got != expected) {
        cerr << "got:\n";
        for (auto& p : got) {
            cerr << '\t' << p.node_id1 << ", " << p.node_id2 << '\n';
        }
        cerr << "expected:\n";
        for (auto& p : expected) {
            cerr << '\t' << p.node_id1 << ", " << p.node_id2 << '\n';
        }
        exit(1);
    }
}

}
