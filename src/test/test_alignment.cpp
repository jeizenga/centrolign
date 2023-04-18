#include <cstdio>
#include <cstdlib>
#include <cstdlib>
#include <unordered_set>
#include <vector>
#include <list>
#include <string>
#include <random>

#include "centrolign/alignment.hpp"
#include "centrolign/utility.hpp"
#include "centrolign/graph.hpp"
#include "centrolign/random_graph.hpp"

using namespace std;
using namespace centrolign;

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

void reverse_alignment(Alignment& aln) {
    for (auto& p : aln) {
        std::swap(p.node_id1, p.node_id2);
    }
}

int rescore(const Alignment& aln, const BaseGraph& graph1, const BaseGraph& graph2,
            const AlignmentParameters<1>& params) {
    int score = 0;
    for (size_t i = 0; i < aln.size(); ++i) {
        if (aln[i].node_id1 != AlignedPair::gap && aln[i].node_id2 != AlignedPair::gap) {
            if (graph1.label(aln[i].node_id1) == graph2.label(aln[i].node_id2)) {
                score += params.match;
            }
            else {
                score -= params.mismatch;
            }
        }
    }
    
    for (size_t i = 0; i < aln.size();) {
        if (aln[i].node_id1 == AlignedPair::gap) {
            size_t j = i + 1;
            while (j < aln.size() && aln[j].node_id1 == AlignedPair::gap) {
                ++j;
            }
            score -= params.gap_open[0] + (j - i) * params.gap_extend[0];
            i = j;
        }
        else {
            ++i;
        }
    }
    
    for (size_t i = 0; i < aln.size();) {
        if (aln[i].node_id2 == AlignedPair::gap) {
            size_t j = i + 1;
            while (j < aln.size() && aln[j].node_id2 == AlignedPair::gap) {
                ++j;
            }
            score -= params.gap_open[0] + (j - i) * params.gap_extend[0];
            i = j;
        }
        else {
            ++i;
        }
    }
    
    return score;
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

vector<vector<uint64_t>> all_paths(const BaseGraph& graph,
                                   uint64_t id_from, uint64_t id_to) {
    
    vector<vector<uint64_t>> paths;
    
    vector<tuple<uint64_t, vector<uint64_t>, size_t>> stack;
    stack.emplace_back(id_from, graph.next(id_from), 0);
    
    while (!stack.empty()) {
        auto& top = stack.back();
        if (get<0>(top) == id_to) {
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

string path_to_string(const BaseGraph& graph, const vector<uint64_t>& path) {
    string seq;
    for (auto i : path) {
        seq.push_back(graph.label(i));
    }
    return seq;
}

void verify_po_poa(const BaseGraph& graph1, const BaseGraph& graph2,
                   const std::vector<uint64_t>& sources1,
                   const std::vector<uint64_t>& sources2,
                   const std::vector<uint64_t>& sinks1,
                   const std::vector<uint64_t>& sinks2,
                   const AlignmentParameters<1>& params) {
    
    bool reachable1 = false;
    for (auto r : sources1) {
        for (auto n : sinks1) {
            if (is_reachable(graph1, r, n)) {
                reachable1 = true;
                break;
            }
        }
    }
    bool reachable2 = false;
    for (auto r : sources2) {
        for (auto n : sinks2) {
            if (is_reachable(graph2, r, n)) {
                reachable2 = true;
                break;
            }
        }
    }
    
    // for now i'm not interested in the unreachable case
    if (!reachable1 || !reachable2) {
        return;
    }
    
    auto po_poa_alignment = po_poa(graph1, graph2, sources1, sources2,
                                   sinks1, sinks2, params);
    
    
    vector<vector<uint64_t>> paths1, paths2;
    for (auto r : sources1) {
        for (auto n : sinks1) {
            for (auto& p : all_paths(graph1, r, n)) {
                paths1.emplace_back(move(p));
            }
        }
    }
    for (auto r : sources2) {
        for (auto n : sinks2) {
            for (auto& p : all_paths(graph2, r, n)) {
                paths2.emplace_back(move(p));
            }
        }
    }
    
    int best_score = -100000000;
    Alignment best_aln;
    
    for (auto& p1 : paths1) {
        
        auto seq1 = path_to_string(graph1, p1);
        for (auto& p2 : paths2) {
            auto seq2 = path_to_string(graph2, p2);
            auto aln = align_nw(seq1, seq2, params);

            for (auto& ap : aln) {
                if (ap.node_id1 != AlignedPair::gap) {
                    ap.node_id1 = p1[ap.node_id1];
                }
                if (ap.node_id2 != AlignedPair::gap) {
                    ap.node_id2 = p2[ap.node_id2];
                }
            }
            
            int s = rescore(aln, graph1, graph2, params);
            if (s > best_score) {
                best_score = s;
                best_aln = aln;
            }
        }
    }
    
    int popoa_score = rescore(po_poa_alignment, graph1, graph2, params);
    if (best_score != popoa_score) {
        cerr << "failed to find optimal alignment with PO-POA: score " << popoa_score << " vs best score " << best_score << "\n";
        cerr << "graph1:\n";
        print_graph(graph1, cerr);
        cerr << "graph2:\n";
        print_graph(graph2, cerr);
        std::cerr << "sources on 1:\n";
        for (auto s : sources1) {
            std::cerr << ' ' << s;
        }
        std::cerr << '\n';
        std::cerr << "sources on 2:\n";
        for (auto s : sources2) {
            std::cerr << ' ' << s;
        }
        std::cerr << '\n';
        std::cerr << "sinks on 1:\n";
        for (auto s : sinks1) {
            std::cerr << ' ' << s;
        }
        std::cerr << '\n';
        std::cerr << "sinks on 2:\n";
        for (auto s : sinks2) {
            std::cerr << ' ' << s;
        }
        std::cerr << '\n';
        check_alignment(po_poa_alignment, best_aln);
    }
}

int main(int argc, char* argv[]) {
     
    uint64_t gap = AlignedPair::gap;
    
    BaseGraph graph1;
    for (char c : "ACGTGCA") {
        if (c != '\0') {
            graph1.add_node(c);
        }
    }
    
    BaseGraph graph2;
    for (char c : "AGTTTGA") {
        if (c != '\0') {
            graph2.add_node(c);
        }
    }
    
    // both graphs will have the same topology (two bubbles)
    for (auto gp : {&graph1, &graph2}) {
        auto& g = *gp;
        g.add_edge(0, 1);
        g.add_edge(0, 2);
        g.add_edge(1, 3);
        g.add_edge(2, 3);
        g.add_edge(3, 4);
        g.add_edge(3, 5);
        g.add_edge(4, 6);
        g.add_edge(5, 6);
    }
    
    AlignmentParameters<1> params;
    params.match = 1;
    params.mismatch = 1;
    params.gap_extend[0] = 1;
    params.gap_open[0] = 1;
    
//    {
//        Alignment alignment = po_poa(graph1, graph2, {0}, {0}, {6}, {6}, params);
//        Alignment expected;
//        expected.emplace_back(0, 0);
//        expected.emplace_back(2, 1);
//        expected.emplace_back(3, 3);
//        expected.emplace_back(4, 5);
//        expected.emplace_back(6, 6);
//
//        check_alignment(alignment, expected);
//    }
//
//    graph1.add_node('T');
//    graph1.add_edge(7, 0);
//    graph2.add_node('T');
//    graph2.add_edge(6, 7);
//
//    {
//        Alignment alignment = po_poa(graph1, graph2, {7}, {0}, {6}, {7}, params);
//        Alignment expected;
//        expected.emplace_back(7, gap);
//        expected.emplace_back(0, 0);
//        expected.emplace_back(2, 1);
//        expected.emplace_back(3, 3);
//        expected.emplace_back(4, 5);
//        expected.emplace_back(6, 6);
//        expected.emplace_back(gap, 7);
//
//        check_alignment(alignment, expected);
//    }
//
//    // flipped
//    {
//        Alignment alignment = po_poa(graph2, graph1, {0}, {7}, {7}, {6}, params);
//        Alignment expected;
//        expected.emplace_back(gap, 7);
//        expected.emplace_back(0, 0);
//        expected.emplace_back(1, 2);
//        expected.emplace_back(3, 3);
//        expected.emplace_back(5, 4);
//        expected.emplace_back(6, 6);
//        expected.emplace_back(7, gap);
//
//        check_alignment(alignment, expected);
//    }
//
//    // test standard alignment
//    {
//        string seq1 = "ATGGCCTGCCGGA";
//        string seq2 = "TATGGTCTGAACCGG";
//
//        // -ATGGCCTG--CCGGA
//        // TATGGTCTGAACCGG-
//
//        auto alignment = align_nw(seq1, seq2, params);
//        Alignment expected;
//        expected.emplace_back(gap, 0);
//        expected.emplace_back(0, 1);
//        expected.emplace_back(1, 2);
//        expected.emplace_back(2, 3);
//        expected.emplace_back(3, 4);
//        expected.emplace_back(4, 5);
//        expected.emplace_back(5, 6);
//        expected.emplace_back(6, 7);
//        expected.emplace_back(7, 8);
//        expected.emplace_back(gap, 9);
//        expected.emplace_back(gap, 10);
//        expected.emplace_back(8, 11);
//        expected.emplace_back(9, 12);
//        expected.emplace_back(10, 13);
//        expected.emplace_back(11, 14);
//        expected.emplace_back(12, gap);
//
//        check_alignment(alignment, expected);
//    }
    
    // cases that came up in randomized testing
    {
        BaseGraph graph1, graph2;
        
        string seq1 = "GAGCGCA";
        string seq2 = "AATAACC";
        vector<pair<int, int>> edges1{
            {0, 5},
            {1, 6},
            {1, 3},
            {1, 2},
            {1, 4},
            {2, 3},
            {2, 5},
            {3, 6},
            {3, 5},
            {4, 6},
            {5, 6}
        };
        vector<pair<int, int>> edges2{
            {0, 3},
            {0, 4},
            {0, 6},
            {1, 3},
            {1, 4},
            {2, 6},
            {3, 5},
            {3, 6},
            {4, 6},
            {5, 6}
        };
        for (int i = 0; i < seq1.size(); ++i) {
            graph1.add_node(seq1[i]);
        }
        for (int i = 0; i < seq2.size(); ++i) {
            graph2.add_node(seq2[i]);
        }
        for (auto e : edges1) {
            graph1.add_edge(e.first, e.second);
        }
        for (auto e : edges2) {
            graph2.add_edge(e.first, e.second);
        }
        
        vector<uint64_t> sources1{3, 6};
        vector<uint64_t> sources2{0};
        vector<uint64_t> sinks1{3, 6};
        vector<uint64_t> sinks2{5};
        
        verify_po_poa(graph1, graph2, sources1, sources2, sinks1, sinks2, params);
    }
    
    // randomized tests
    
    random_device rd;
    default_random_engine gen(rd());
    
    uniform_int_distribution<int> nodes_distr(5, 10);
    uniform_int_distribution<int> edges_distr(8, 18);
    uniform_int_distribution<int> source_sink_distr(1, 2);
    
    size_t num_trials = 500;
    for (size_t i = 0; i < num_trials; ++i) {
        BaseGraph graph1 = random_graph(nodes_distr(gen), edges_distr(gen), gen);
        BaseGraph graph2 = random_graph(nodes_distr(gen), edges_distr(gen), gen);

        uniform_int_distribution<int> graph1_distr(0, graph1.node_size() - 1);
        uniform_int_distribution<int> graph2_distr(0, graph2.node_size() - 1);

        std::vector<uint64_t> sources1, sources2, sinks1, sinks2;

        for (auto r : {
            make_pair(&sources1, &graph1_distr),
            make_pair(&sinks1, &graph1_distr),
            make_pair(&sources2, &graph2_distr),
            make_pair(&sinks2, &graph2_distr)}) {
            int size = source_sink_distr(gen);
            for (int j = 0; j < size; ++j) {
                r.first->push_back((*r.second)(gen));
            }
            sort(r.first->begin(), r.first->end());
            r.first->resize(unique(r.first->begin(), r.first->end()) - r.first->begin());
        }

        verify_po_poa(graph1, graph2, sources1, sources2, sinks1, sinks2, params);
    }
    
    
    cerr << "passed all tests!" << endl;
}
