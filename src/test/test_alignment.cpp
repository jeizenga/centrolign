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
#include "centrolign/test_util.hpp"

using namespace std;
using namespace centrolign;

void reverse_alignment(Alignment& aln) {
    for (auto& p : aln) {
        std::swap(p.node_id1, p.node_id2);
    }
}


int rescore(const Alignment& aln, const BaseGraph& graph1, const BaseGraph& graph2,
            const AlignmentParameters<1>& params, bool wfa_style) {
    
    auto local_params = params;
    
    if (wfa_style) {
        local_params.match = 0;
        local_params.mismatch = 2 * (params.match + params.mismatch);
        local_params.gap_open[0] = 2 * params.gap_open[0];
        local_params.gap_extend[0] = 2 * params.gap_extend[0] + params.match;
    }
    
    int score = 0;
    for (size_t i = 0; i < aln.size(); ++i) {
        if (aln[i].node_id1 != AlignedPair::gap && aln[i].node_id2 != AlignedPair::gap) {
            if (graph1.label(aln[i].node_id1) == graph2.label(aln[i].node_id2)) {
                score += local_params.match;
            }
            else {
                score -= local_params.mismatch;
            }
        }
    }
    
    for (size_t i = 0; i < aln.size();) {
        if (aln[i].node_id1 == AlignedPair::gap) {
            size_t j = i + 1;
            while (j < aln.size() && aln[j].node_id1 == AlignedPair::gap) {
                ++j;
            }
            score -= local_params.gap_open[0] + (j - i) * local_params.gap_extend[0];
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
            score -= local_params.gap_open[0] + (j - i) * local_params.gap_extend[0];
            i = j;
        }
        else {
            ++i;
        }
    }
    
    return score;
}

bool alignment_is_valid(const Alignment& aln, const BaseGraph& graph1, const BaseGraph& graph2,
                        const std::vector<uint64_t>& sources1,
                        const std::vector<uint64_t>& sources2,
                        const std::vector<uint64_t>& sinks1,
                        const std::vector<uint64_t>& sinks2) {
    
    for (int i = 0; i < aln.size(); ++i) {
        if (aln[i].node_id1 != AlignedPair::gap) {
            if (find(sources1.begin(), sources1.end(), aln[i].node_id1) == sources1.end()) {
                cerr << "alignment starts with non-source on graph 1\n";
                return false;
            }
            break;
        }
    }
    for (int i = 0; i < aln.size(); ++i) {
        if (aln[i].node_id2 != AlignedPair::gap) {
            if (find(sources2.begin(), sources2.end(), aln[i].node_id2) == sources2.end()) {
                cerr << "alignment starts with non-source on graph 2\n";
                return false;
            }
            break;
        }
    }
    for (int i = aln.size() - 1; i >= 0; --i) {
        if (aln[i].node_id1 != AlignedPair::gap) {
            if (find(sinks1.begin(), sinks1.end(), aln[i].node_id1) == sinks1.end()) {
                cerr << "alignment starts with non-source on graph 1\n";
                return false;
            }
            break;
        }
    }
    for (int i = aln.size() - 1; i >= 0; --i) {
        if (aln[i].node_id2 != AlignedPair::gap) {
            if (find(sinks2.begin(), sinks2.end(), aln[i].node_id2) == sinks2.end()) {
                cerr << "alignment starts with non-source on graph 2\n";
                return false;
            }
            break;
        }
    }
    uint64_t prev_id1 = -1, prev_id2 = -1;
    for (auto aln_pair : aln) {
        if (aln_pair.node_id1 != AlignedPair::gap) {
            if (prev_id1 != -1) {
                auto next = graph1.next(prev_id1);
                if (find(next.begin(), next.end(), aln_pair.node_id1) == next.end()) {
                    cerr << "alignments takes non existent graph 1 edge " << prev_id1 << " " << aln_pair.node_id1 << '\n';
                    return false;
                }
            }
            prev_id1 = aln_pair.node_id1;
        }
        if (aln_pair.node_id2 != AlignedPair::gap) {
            if (prev_id2 != -1) {
                auto next = graph2.next(prev_id2);
                if (find(next.begin(), next.end(), aln_pair.node_id2) == next.end()) {
                    cerr << "alignments takes non existent graph 2 edge " << prev_id2 << " " << aln_pair.node_id2 << '\n';
                    return false;
                }
            }
            prev_id2 = aln_pair.node_id2;
        }
    }
    return true;
}

void verify_wfa_po_poa(const BaseGraph& graph1, const BaseGraph& graph2,
                       const std::vector<uint64_t>& sources1,
                       const std::vector<uint64_t>& sources2,
                       const std::vector<uint64_t>& sinks1,
                       const std::vector<uint64_t>& sinks2,
                       const AlignmentParameters<1>& params) {
    // TODO: copypasta
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
    
    auto dump_data = [&]() {
        cerr << "graph1:\n";
        cerr << cpp_representation(graph1, "graph1");
        cerr << "graph2:\n";
        cerr << cpp_representation(graph2, "graph2");
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
    };
    
    auto po_poa_alignment = po_poa(graph1, graph2, sources1, sources2,
                                   sinks1, sinks2, params);
    if (!alignment_is_valid(po_poa_alignment, graph1, graph2, sources1, sources2,
                            sinks1, sinks2)) {
        dump_data();
        exit(1);
    }
    auto wfa_po_poa_alignment = wfa_po_poa(graph1, graph2, sources1, sources2,
                                           sinks1, sinks2, params);
    if (!alignment_is_valid(wfa_po_poa_alignment, graph1, graph2, sources1, sources2,
                            sinks1, sinks2)) {
        cerr << "invalid wfa alignment\n";
        dump_data();
        exit(1);
    }
    // these don't have guarantees, but we want to make sure they runs
    auto pwfa_po_poa_alignment = pwfa_po_poa(graph1, graph2, sources1, sources2,
                                             sinks1, sinks2, params, 5);
    if (!alignment_is_valid(pwfa_po_poa_alignment, graph1, graph2, sources1, sources2,
                            sinks1, sinks2)) {
        cerr << "invalid pwfa alignment\n";
        dump_data();
        exit(1);
    }
    auto del_wfa_alignment = deletion_wfa_po_poa(graph1, graph2, sources1, sources2,
                                                 sinks1, sinks2, params);
    if (!alignment_is_valid(del_wfa_alignment, graph1, graph2, sources1, sources2,
                            sinks1, sinks2)) {
        cerr << "invalid deletion wfa alignment\n";
        dump_data();
        exit(1);
    }
    
    // TODO: it's possible for WFA to prefer solutions with lower WFA scores
    // but higher standard scores (because the lengths of the aligned paths are not
    // necessarily constant)
    int score_po_poa = rescore(po_poa_alignment, graph1, graph2, params, true);
    int score_wfa_po_poa = rescore(wfa_po_poa_alignment, graph1, graph2, params, true);
    if (score_po_poa > score_wfa_po_poa) {
        cerr << "failed to find optimal WFA-scoring alignment, got " << score_wfa_po_poa << ", but could have gotten " << score_po_poa << "\n";
        dump_data();
        check_alignment(wfa_po_poa_alignment, po_poa_alignment);
    }
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
            
            int s = rescore(aln, graph1, graph2, params, false);
            if (s > best_score) {
                best_score = s;
                best_aln = aln;
            }
        }
    }
    
    int popoa_score = rescore(po_poa_alignment, graph1, graph2, params, false);
    if (best_score != popoa_score) {
        cerr << "failed to find optimal alignment with PO-POA: score " << popoa_score << " vs best score " << best_score << "\n";
        std::cerr << cpp_representation(graph1, "graph1") << '\n';
        std::cerr << cpp_representation(graph2, "graph2") << '\n';
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

void test_induced_alignment(const BaseGraph& graph, uint64_t path_id1, uint64_t path_id2,
                            const Alignment& correct) {
    
    auto aln = induced_pairwise_alignment(graph, path_id1, path_id2);
    check_alignment(aln, correct);
    
    auto flipped_correct = correct;
    for (auto& aln_pair : flipped_correct) {
        std::swap(aln_pair.node_id1, aln_pair.node_id2);
    }
    auto flipped_aln = induced_pairwise_alignment(graph, path_id2, path_id1);
    // a pain because of the ordering of gaps
//    check_alignment(flipped_aln, flipped_correct);
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

    {
        Alignment alignment = po_poa(graph1, graph2, {0}, {0}, {6}, {6}, params);
        Alignment alignment_wfa = wfa_po_poa(graph1, graph2, {0}, {0}, {6}, {6}, params);
        Alignment alignment_pwfa = pwfa_po_poa(graph1, graph2, {0}, {0}, {6}, {6}, params, 4);
        Alignment expected;
        expected.emplace_back(0, 0);
        expected.emplace_back(2, 1);
        expected.emplace_back(3, 3);
        expected.emplace_back(4, 5);
        expected.emplace_back(6, 6);

        check_alignment(alignment, expected);
        check_alignment(alignment_wfa, expected);
        check_alignment(alignment_pwfa, expected);
    }

    graph1.add_node('T');
    graph1.add_edge(7, 0);
    graph2.add_node('T');
    graph2.add_edge(6, 7);

    {
        Alignment alignment = po_poa(graph1, graph2, {7}, {0}, {6}, {7}, params);
        Alignment alignment_wfa = wfa_po_poa(graph1, graph2, {7}, {0}, {6}, {7}, params);
        Alignment alignment_pwfa = pwfa_po_poa(graph1, graph2, {7}, {0}, {6}, {7}, params, 4);
        Alignment expected;
        expected.emplace_back(7, gap);
        expected.emplace_back(0, 0);
        expected.emplace_back(2, 1);
        expected.emplace_back(3, 3);
        expected.emplace_back(4, 5);
        expected.emplace_back(6, 6);
        expected.emplace_back(gap, 7);

        check_alignment(alignment, expected);
        check_alignment(alignment_wfa, expected);
        check_alignment(alignment_pwfa, expected);
    }

    // flipped
    {
        Alignment alignment = po_poa(graph2, graph1, {0}, {7}, {7}, {6}, params);
        Alignment alignment_wfa = wfa_po_poa(graph2, graph1, {0}, {7}, {7}, {6}, params);
        Alignment alignment_pwfa = pwfa_po_poa(graph2, graph1, {0}, {7}, {7}, {6}, params, 4);
        Alignment expected;
        expected.emplace_back(gap, 7);
        expected.emplace_back(0, 0);
        expected.emplace_back(1, 2);
        expected.emplace_back(3, 3);
        expected.emplace_back(5, 4);
        expected.emplace_back(6, 6);
        expected.emplace_back(7, gap);

        check_alignment(alignment, expected);
        check_alignment(alignment_wfa, expected);
        check_alignment(alignment_pwfa, expected);
    }

    // test standard alignment
    {
        string seq1 = "ATGGCCTGCCGGA";
        string seq2 = "TATGGTCTGAACCGG";

        // -ATGGCCTG--CCGGA
        // TATGGTCTGAACCGG-

        auto alignment = align_nw(seq1, seq2, params);
        Alignment expected;
        expected.emplace_back(gap, 0);
        expected.emplace_back(0, 1);
        expected.emplace_back(1, 2);
        expected.emplace_back(2, 3);
        expected.emplace_back(3, 4);
        expected.emplace_back(4, 5);
        expected.emplace_back(5, 6);
        expected.emplace_back(6, 7);
        expected.emplace_back(7, 8);
        expected.emplace_back(gap, 9);
        expected.emplace_back(gap, 10);
        expected.emplace_back(8, 11);
        expected.emplace_back(9, 12);
        expected.emplace_back(10, 13);
        expected.emplace_back(11, 14);
        expected.emplace_back(12, gap);

        check_alignment(alignment, expected);
    }

    // deletion wfa
    {
        BaseGraph short_graph;
        for (auto c : std::string("ACTGCA")) {
            short_graph.add_node(c);
        }

        std::vector<std::pair<int, int>> short_graph_edges{
            {0, 1},
            {0, 2},
            {1, 3},
            {2, 3},
            {3, 4},
            {3, 5},
            {4, 5}
        };

        for (auto e : short_graph_edges) {
            short_graph.add_edge(e.first, e.second);
        }

        BaseGraph long_graph;
        for (auto c : std::string("ACG") + string(50, 'T') + string("CA")) {
            long_graph.add_node(c);
        }

        for (size_t i = 0; i + 1 < long_graph.node_size(); ++i) {
            long_graph.add_edge(i, i + 1);
        }


        vector<uint64_t> sources_short{0};
        vector<uint64_t> sources_long{0};
        vector<uint64_t> sinks_short{5};
        vector<uint64_t> sinks_long{long_graph.node_size() - 1};

        AlignmentParameters<1> params;
        params.match = 1;
        params.mismatch = 1;
        params.gap_open[0] = 1;
        params.gap_extend[0] = 1;

        auto aln = deletion_wfa_po_poa(short_graph, long_graph, sources_short, sources_long,
                                       sinks_short, sinks_long, params);
        if (!alignment_is_valid(aln, short_graph, long_graph, sources_short, sources_long,
                                sinks_short, sinks_long)) {
            std::cerr << "failed deletion wfa test\n";
            exit(1);
        }

    }

    // cases that came up in randomized testing
    {
        BaseGraph graph1;
        for (auto c : std::string("ACACACACACACGGC")) {
            graph1.add_node(c);
        }

        std::vector<std::pair<int, int>> graph1_edges{
            {0, 1},
            {0, 12},
            {0, 2},
            {1, 2},
            {2, 3},
            {2, 4},
            {2, 13},
            {2, 5},
            {3, 4},
            {3, 13},
            {3, 5},
            {4, 5},
            {5, 6},
            {5, 6},
            {6, 7},
            {6, 14},
            {7, 8},
            {8, 9},
            {8, 11},
            {8, 11},
            {10, 1},
            {10, 12},
            {10, 2},
            {10, 6},
            {12, 2},
            {13, 5},
            {14, 8}
        };

        for (auto e : graph1_edges) {
            graph1.add_edge(e.first, e.second);
        }

        BaseGraph graph2;
        for (auto c : std::string("ACCACCACCACCACCGAAC")) {
            graph2.add_node(c);
        }

        std::vector<std::pair<int, int>> graph2_edges{
            {0, 1},
            {1, 2},
            {1, 18},
            {2, 3},
            {3, 4},
            {4, 5},
            {5, 6},
            {6, 7},
            {6, 10},
            {7, 8},
            {7, 12},
            {8, 9},
            {9, 10},
            {9, 15},
            {10, 11},
            {11, 12},
            {12, 13},
            {13, 14},
            {13, 17},
            {15, 11},
            {16, 1},
            {18, 3}
        };

        for (auto e : graph2_edges) {
            graph2.add_edge(e.first, e.second);
        }


        vector<uint64_t> sources1{3, 12};
        vector<uint64_t> sources2{0};
        vector<uint64_t> sinks1{13};
        vector<uint64_t> sinks2{13, 16};

        verify_po_poa(graph1, graph2, sources1, sources2, sinks1, sinks2, params);
        verify_wfa_po_poa(graph1, graph2, sources1, sources2, sinks1, sinks2, params);
    }

    {
        BaseGraph graph1;
        for (auto c : std::string("AACAACAACAACAACTTCGTC")) {
            graph1.add_node(c);
        }

        std::vector<std::pair<int, int>> graph1_edges{
            {0, 1},
            {0, 19},
            {1, 2},
            {2, 3},
            {2, 4},
            {3, 4},
            {4, 5},
            {5, 6},
            {5, 18},
            {6, 7},
            {7, 8},
            {7, 20},
            {8, 9},
            {8, 15},
            {9, 10},
            {9, 16},
            {10, 11},
            {11, 12},
            {12, 13},
            {13, 14},
            {13, 14},
            {15, 10},
            {15, 16},
            {16, 11},
            {16, 14},
            {17, 1},
            {17, 19},
            {18, 7},
            {19, 2},
            {19, 4},
            {20, 9},
            {20, 15}
        };
        for (auto e : graph1_edges) {
            graph1.add_edge(e.first, e.second);
        }

        BaseGraph graph2;
        for (auto c : std::string("AAAAGAAAAAAAATACTGG")) {
            graph2.add_node(c);
        }

        std::vector<std::pair<int, int>> graph2_edges{
            {0, 1},
            {1, 2},
            {1, 17},
            {2, 3},
            {2, 4},
            {3, 4},
            {4, 5},
            {4, 18},
            {4, 6},
            {5, 6},
            {5, 16},
            {6, 7},
            {7, 8},
            {8, 9},
            {9, 10},
            {9, 16},
            {9, 12},
            {10, 11},
            {10, 15},
            {10, 12},
            {11, 12},
            {12, 13},
            {13, 14},
            {15, 12},
            {16, 11},
            {16, 15},
            {16, 12},
            {17, 3},
            {17, 4},
            {18, 6}
        };

        for (auto e : graph2_edges) {
            graph2.add_edge(e.first, e.second);
        }

        vector<uint64_t> sources1{3};
        vector<uint64_t> sources2{4, 9};
        vector<uint64_t> sinks1{10, 19};
        vector<uint64_t> sinks2{12};

        verify_po_poa(graph1, graph2, sources1, sources2, sinks1, sinks2, params);
        verify_wfa_po_poa(graph1, graph2, sources1, sources2, sinks1, sinks2, params);
    }
    
    {
        BaseGraph graph1;
        for (auto c : std::string("AAAAAAAAAGCTT")) {
            graph1.add_node(c);
        }
        
        std::vector<std::pair<int, int>> graph1_edges{
            {0, 1},
            {0, 2},
            {1, 2},
            {2, 3},
            {3, 4},
            {3, 5},
            {3, 10},
            {3, 12},
            {3, 6},
            {3, 7},
            {3, 6},
            {4, 5},
            {4, 10},
            {4, 12},
            {4, 6},
            {5, 6},
            {6, 7},
            {7, 8},
            {7, 11},
            {8, 9},
            {10, 6},
            {11, 9},
            {12, 6}
        };
        
        for (auto e : graph1_edges) {
            graph1.add_edge(e.first, e.second);
        }

        BaseGraph graph2;
        for (auto c : std::string("ATCAAACAAACAACCAAACACCAAG")) {
            graph2.add_node(c);
        }
        
        std::vector<std::pair<int, int>> graph2_edges{
            {0, 1},
            {1, 2},
            {1, 5},
            {2, 3},
            {2, 20},
            {3, 4},
            {3, 5},
            {4, 5},
            {5, 6},
            {5, 21},
            {5, 7},
            {6, 7},
            {7, 8},
            {8, 9},
            {8, 23},
            {8, 24},
            {9, 10},
            {10, 11},
            {11, 12},
            {12, 13},
            {13, 14},
            {14, 15},
            {15, 16},
            {16, 17},
            {16, 22},
            {16, 19},
            {17, 18},
            {17, 19},
            {18, 19},
            {20, 4},
            {21, 7},
            {22, 18},
            {22, 19},
            {23, 10},
            {24, 10}
        };
        
        for (auto e : graph2_edges) {
            graph2.add_edge(e.first, e.second);
        }

        vector<uint64_t> sources1{1};
        vector<uint64_t> sources2{0, 6};
        vector<uint64_t> sinks1{11};
        vector<uint64_t> sinks2{7};
        
        verify_po_poa(graph1, graph2, sources1, sources2, sinks1, sinks2, params);
        verify_wfa_po_poa(graph1, graph2, sources1, sources2, sinks1, sinks2, params);
    }

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
        verify_wfa_po_poa(graph1, graph2, sources1, sources2, sinks1, sinks2, params);
    }

    {
        BaseGraph graph1;
        for (auto c : std::string("CCACCCCCCCACT")) {
            graph1.add_node(c);
        }

        std::vector<std::pair<int, int>> graph1_edges{
            {0, 1},
            {0, 4},
            {1, 2},
            {2, 3},
            {3, 4},
            {3, 6},
            {3, 11},
            {4, 5},
            {5, 6},
            {5, 10},
            {6, 7},
            {6, 12},
            {7, 8},
            {8, 9},
            {8, 11},
            {10, 7},
            {10, 12},
            {12, 8}
        };

        for (auto e : graph1_edges) {
            graph1.add_edge(e.first, e.second);
        }

        BaseGraph graph2;
        for (auto c : std::string("CCCCCCTTCCCCCCA")) {
            graph2.add_node(c);
        }

        std::vector<std::pair<int, int>> graph2_edges{
            {0, 1},
            {0, 10},
            {0, 2},
            {1, 2},
            {2, 3},
            {3, 4},
            {3, 11},
            {4, 5},
            {4, 7},
            {4, 12},
            {5, 6},
            {6, 7},
            {7, 8},
            {7, 12},
            {7, 14},
            {8, 9},
            {8, 13},
            {10, 2},
            {11, 5},
            {12, 9},
            {12, 13},
            {14, 9},
            {14, 13}
        };

        for (auto e : graph2_edges) {
            graph2.add_edge(e.first, e.second);
        }

        vector<uint64_t> sources1{1, 2};
        vector<uint64_t> sources2{1, 7};
        vector<uint64_t> sinks1{1, 3};
        vector<uint64_t> sinks2{5, 14};

        verify_po_poa(graph1, graph2, sources1, sources2, sinks1, sinks2, params);
        verify_wfa_po_poa(graph1, graph2, sources1, sources2, sinks1, sinks2, params);
    }

    // randomized tests

    random_device rd;
    default_random_engine gen(rd());

    uniform_int_distribution<int> nodes_distr(5, 10);
    uniform_int_distribution<int> edges_distr(8, 18);
    uniform_int_distribution<int> source_sink_distr(1, 2);
    uniform_int_distribution<int> challenge_nodes_distr(10, 25);

    size_t num_trials = 500;
    for (size_t i = 0; i < num_trials; ++i) {
        BaseGraph graph1, graph2;
        if (i % 2 == 0) {
            graph1 = random_graph(nodes_distr(gen), edges_distr(gen), gen);
            graph2 = random_graph(nodes_distr(gen), edges_distr(gen), gen);
        }
        else {
            graph1 = random_challenge_graph(challenge_nodes_distr(gen), gen);
            graph2 = random_challenge_graph(challenge_nodes_distr(gen), gen);
        }

        uniform_int_distribution<int> graph1_distr(0, graph1.node_size() - 1);
        uniform_int_distribution<int> graph2_distr(0, graph2.node_size() - 1);

        std::vector<uint64_t> sources1, sources2, sinks1, sinks2;

        for (auto r : {
            make_pair(&sources1, &graph1_distr),
            make_pair(&sinks1, &graph1_distr),
            make_pair(&sources2, &graph2_distr),
            make_pair(&sinks2, &graph2_distr)
        }) {
            int size = source_sink_distr(gen);
            for (int j = 0; j < size; ++j) {
                r.first->push_back((*r.second)(gen));
            }
            sort(r.first->begin(), r.first->end());
            r.first->resize(unique(r.first->begin(), r.first->end()) - r.first->begin());
        }

        verify_po_poa(graph1, graph2, sources1, sources2, sinks1, sinks2, params);
        verify_wfa_po_poa(graph1, graph2, sources1, sources2, sinks1, sinks2, params);
    }

    // tests for induced pairwise alignments

    {
        BaseGraph graph;
        for (auto c : string("AAACAA")) {
            graph.add_node(c);
        }
        vector<pair<int, int>> edges{
            {0, 1},
            {0, 5},
            {1, 2},
            {1, 3},
            {2, 4},
            {3, 4},
            {4, 5}
        };
        for (auto e : edges) {
            graph.add_edge(e.first, e.second);
        }

        vector<vector<int>> paths{
            {0, 1, 2, 4, 5},
            {0, 1, 3, 4, 5},
            {0, 5}
        };

        int p = 0;
        for (auto& path : paths) {
            auto path_id = graph.add_path(to_string(p++));
            for (auto n : path) {
                graph.extend_path(path_id, n);
            }
        }

        Alignment aln01{
            {0, 0},
            {1, 1},
            {2, 2},
            {3, 3},
            {4, 4}
        };

        Alignment aln02_12{
            {0, 0},
            {1, AlignedPair::gap},
            {2, AlignedPair::gap},
            {3, AlignedPair::gap},
            {4, 1}
        };

        test_induced_alignment(graph, 0, 1, aln01);
        test_induced_alignment(graph, 0, 2, aln02_12);
        test_induced_alignment(graph, 1, 2, aln02_12);
    }

    {
        BaseGraph graph;
        for (auto c : string("AACCTAA")) {
            graph.add_node(c);
        }
        vector<pair<int, int>> edges{
            {0, 1},
            {1, 2},
            {1, 4},
            {2, 3},
            {3, 5},
            {4, 5},
            {5, 6}
        };
        for (auto e : edges) {
            graph.add_edge(e.first, e.second);
        }

        vector<vector<int>> paths{
            {0, 1, 2, 3, 5, 6},
            {1, 4, 5}
        };

        int p = 0;
        for (auto& path : paths) {
            auto path_id = graph.add_path(to_string(p++));
            for (auto n : path) {
                graph.extend_path(path_id, n);
            }
        }

        Alignment aln01{
            {0, AlignedPair::gap},
            {1, 0},
            {2, AlignedPair::gap},
            {3, AlignedPair::gap},
            {AlignedPair::gap, 1},
            {4, 2},
            {5, AlignedPair::gap}
        };

        test_induced_alignment(graph, 0, 1, aln01);
    }
    
    cerr << "passed all tests!" << endl;
}
