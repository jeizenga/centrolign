 #include "centrolign/test_util.hpp"

#include <set>
#include <map>

namespace centrolign {

using namespace std;

bool is_valid_path(const BaseGraph& graph, std::vector<uint64_t>& path) {
    
    for (size_t i = 1; i < path.size(); ++i) {
        bool found = false;
        for (auto n : graph.next(path[i - 1])) {
            if (n == path[i]) {
                found = true;
            }
        }
        if (!found) {
            return false;
        }
    }
    return true;
}

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

bool paths_match(const BaseGraph& graph, const BaseGraph& graph2) {
    if (graph.path_size() != graph2.path_size()) {
        cerr << "graphs do not have same number of paths.\n";
        return false;
    }
    for (uint64_t path_id = 0; path_id < graph.path_size(); ++path_id) {
        
        auto other_path_id = graph2.path_id(graph.path_name(path_id));
        
        if (other_path_id == -1) {
            cerr << "path " << graph.path_name(path_id) << " is missing in second graph.\n";
            return false;
        }
        
        auto seq1 = path_to_string(graph, graph.path(path_id));
        auto seq2 = path_to_string(graph2, graph2.path(other_path_id));
        
        if (seq1 != seq2) {
            cerr << "did not find matching path sequence for path " << graph.path_name(path_id) << " in second graph.\n";
            return false;
        }
    }
    return true;
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

std::string cpp_representation(const BaseGraph& graph, const std::string& name) {
    std::vector<pair<int, int>> edges;
    for (int n = 0; n < graph.node_size(); ++n) {
        for (auto m : graph.next(n)) {
            edges.emplace_back(n, m);
        }
    }
    
    stringstream strm;
    strm << "BaseGraph " << name << ";\n";
    strm << "for (auto c : std::string(\"";
    for (int n = 0; n < graph.node_size(); ++n) {
        char c = graph.label(n);
        if (c <= 5) {
            c = decode_base(c);
        }
        strm << c;
    }
    strm << "\")) {\n";
    strm << "    " << name << ".add_node(c);\n";
    strm << "}\n\n";
    strm << "std::vector<std::pair<int, int>> " << name << "_edges{\n";
    for (size_t i = 0; i < edges.size(); ++i) {
        strm << "    {" << edges[i].first << ", " << edges[i].second << "}";
        if (i + 1 != edges.size()) {
            strm << ",";
        }
        strm << "\n";
    }
    strm << "};\n\n";
    strm << "std::vector<std::vector<int>> " << name << "_paths{\n";
    for (size_t p = 0; p < graph.path_size(); ++p) {
        strm << "    {";
        auto nodes = graph.path(p);
        for (size_t i = 0; i < nodes.size(); ++i) {
            strm << nodes[i];
            if (i + 1 != nodes.size()) {
                strm << ", ";
            }
        }
        strm << "}";
        if (p + 1 != graph.path_size()) {
            strm << ",";
        }
        strm << "\n";
    }
    strm << "};\n\n";
    
    strm << "for (auto e : " << name << "_edges) {\n";
    strm << "   " << name << ".add_edge(e.first, e.second);\n";
    strm << "}\n\n";
    
    strm << "for (size_t i = 0; i < " << name << "_paths.size(); ++i) {\n";
    strm << "    auto p = " << name << ".add_path(std::to_string(i));\n";
    strm << "    for (auto n : " << name << "_paths[i]) {\n";
    strm << "        " << name << ".extend_path(p, n);\n";
    strm << "    }\n";
    strm << "}\n";
    
    return strm.str();
}

std::string pretty_alignment(const Alignment& aln, const std::string& seq1, const std::string& seq2) {
    
    std::stringstream strm;
    
    const size_t line_length = 80;
    size_t idx1 = 0, idx2 = 0;
    size_t gap_size1 = 0, gap_size2 = 0;
    for (size_t i = 0; i < aln.size(); i += line_length) {
        
        strm << idx1 << '\t';
        size_t n = std::min<size_t>(line_length, aln.size() - i);
        vector<size_t> gap_endings;
        for (size_t j = 0; j < n; ++j) {
            if (aln[i + j].node_id1 == AlignedPair::gap) {
                strm << '-';
                ++gap_size1;
            }
            else {
                char c = seq1[aln[i + j].node_id1];
                if (c <= 5) {
                    c = decode_base(c);
                }
                strm << c;
                assert(idx1 == aln[i + j].node_id1);
                if (gap_size1 != 0) {
                    gap_endings.push_back(gap_size1);
                }
                gap_size1 = 0;
                ++idx1;
            }
        }
        if (i + line_length >= aln.size() && gap_size1 != 0) {
            gap_endings.push_back(gap_size1);
        }
        for (size_t j = 0; j < gap_endings.size(); ++j) {
            if (j) {
                strm << ',';
            }
            strm << ' ' << gap_endings[j];
        }
        gap_endings.clear();
        strm << '\n';
        strm << '\t';
        for (size_t j = 0; j < n; ++j) {
            if (aln[i + j].node_id1 == AlignedPair::gap || aln[i + j].node_id2 == AlignedPair::gap ||
                seq1[aln[i + j].node_id1] == seq2[aln[i + j].node_id2]) {
                strm << ' ';
            }
            else {
                strm << 'X';
            }
        }
        strm << '\n';
        strm << idx2 << '\t';
        for (size_t j = 0; j < n; ++j) {
            if (aln[i + j].node_id2 == AlignedPair::gap) {
                strm << '-';
                ++gap_size2;
            }
            else {
                char c = seq2[aln[i + j].node_id2];
                if (c <= 5) {
                    c = decode_base(c);
                }
                strm << c;
                assert(idx2 == aln[i + j].node_id2);
                if (gap_size2 != 0) {
                    gap_endings.push_back(gap_size2);
                }
                gap_size2 = 0;
                ++idx2;
            }
        }
        for (size_t j = 0; j < gap_endings.size(); ++j) {
            if (j) {
                strm << ',';
            }
            else {
                strm << ' ' << gap_endings[j];
            }
        }
        strm << '\n';
        if (i + line_length < aln.size()) {
            strm << '\n';
        }
    }
    
    return strm.str();
}

std::string pretty_alignment(const Alignment& aln, const BaseGraph graph1, const BaseGraph& graph2) {
    
    std::string str1, str2;
    Alignment str_aln;
    for (const auto& ap : aln) {
        str_aln.emplace_back();
        if (ap.node_id1 == AlignedPair::gap) {
            str_aln.back().node_id1 = AlignedPair::gap;
        }
        else {
            str_aln.back().node_id1 = str1.size();
            str1.push_back(graph1.label(ap.node_id1));
        }
        if (ap.node_id2 == AlignedPair::gap) {
            str_aln.back().node_id2 = AlignedPair::gap;
        }
        else {
            str_aln.back().node_id2 = str2.size();
            str2.push_back(graph2.label(ap.node_id2));
        }
    }
    
    return pretty_alignment(str_aln, str1, str2);
}

}
