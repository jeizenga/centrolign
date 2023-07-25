#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <random>

#include "centrolign/anchorer.hpp"
#include "centrolign/chain_merge.hpp"
#include "centrolign/path_merge.hpp"
#include "centrolign/modify_graph.hpp"
#include "centrolign/graph.hpp"
#include "centrolign/utility.hpp"
#include "centrolign/test_util.hpp"
#include "centrolign/match_finder.hpp"

using namespace std;
using namespace centrolign;


class TestAnchorer : public Anchorer {
public:
    using Anchorer::AnchorGraph;
    using Anchorer::exhaustive_chain_dp;
    using Anchorer::sparse_chain_dp;
    using Anchorer::anchor_weight;
    using Anchorer::sparse_affine_chain_dp;
    using Anchorer::gap_open;
    using Anchorer::gap_extend;
    using Anchorer::length_scale;
    using Anchorer::edge_weight;
    using Anchorer::post_switch_distances;
};


// DFS to walk out all fixed-length paths
unordered_map<string, vector<vector<uint64_t>>> k_mer_walks(const BaseGraph& graph, size_t k) {
    
    assert(k != 0);
    
    unordered_map<string, vector<vector<uint64_t>>> to_return;
    
    for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
        
        // records of (node id, next edge)
        vector<pair<uint64_t, size_t>> stack;
        stack.emplace_back(node_id, 0);
        
        while (!stack.empty()) {
            
            if (stack.size() == k) {
                
                string seq;
                vector<uint64_t> walk;
                for (auto& rec : stack) {
                    seq.push_back(graph.label(rec.first));
                    walk.push_back(rec.first);
                }
                to_return[seq].push_back(walk);
                
                stack.pop_back();
            }
            else if (stack.back().second == graph.next_size(stack.back().first)) {
                stack.pop_back();
            }
            else {
                uint64_t next = graph.next(stack.back().first)[stack.back().second++];
                stack.emplace_back(next, 0);
            }
        }
        
    }
    
    return to_return;
}

std::vector<match_set_t> generate_anchor_set(const BaseGraph& graph1,
                                             const BaseGraph& graph2,
                                             size_t k) {
    
    auto kmers1 = k_mer_walks(graph1, k);
    auto kmers2 = k_mer_walks(graph2, k);
    
    std::vector<match_set_t> anchors;
    for (pair<const string, vector<vector<uint64_t>>>& entry : kmers1) {
        if (kmers2.count(entry.first)) {
            anchors.emplace_back();
            anchors.back().walks1 = move(entry.second);
            anchors.back().walks2 = move(kmers2[entry.first]);
        }
    }
    
    return anchors;
}

void print_anchor_set(const vector<match_set_t>& anchors, size_t i) {
    cerr << i << "\twalks on graph1:\n";
    for (const auto& walk : anchors[i].walks1) {
        cerr << "\t\t";
        for (size_t j = 0; j < walk.size(); ++j) {
            if (j) {
                cerr << ',';
            }
            cerr << walk[j];
        }
        cerr << '\n';
    }
    cerr << "\twalks on graph2:\n";
    for (const auto& walk : anchors[i].walks2) {
        cerr << "\t\t";
        for (size_t j = 0; j < walk.size(); ++j) {
            if (j) {
                cerr << ',';
            }
            cerr << walk[j];
        }
        cerr << '\n';
    }
}

template<class XMerge>
void print_chain(const TestAnchorer& anchorer, const vector<anchor_t>& chain, bool affine,
                 const XMerge& xmerge1, const XMerge& xmerge2,
                 vector<vector<size_t>>& switch_dists1, vector<vector<size_t>>& switch_dists2) {
    for (size_t i = 0; i < chain.size(); ++i) {
        auto& link = chain[i];
        cerr << i << " (count1 " << link.count1 << ", count2 " << link.count2 << ", score " << anchorer.anchor_weight(link.count1, link.count2, link.walk1.size()) << "):\n";
        for (auto walk : {link.walk1, link.walk2}) {
            cerr << '\t';
            for (size_t j = 0; j < walk.size(); ++j) {
                if (j) {
                    cerr << ',';
                }
                cerr << walk[j];
            }
            cerr << '\n';
        }
        if (affine && i + 1 != chain.size()) {
            double wt = anchorer.edge_weight(chain[i].walk1.back(), chain[i + 1].walk1.front(),
                                             chain[i].walk2.back(), chain[i + 1].walk2.front(),
                                             xmerge1, xmerge2, switch_dists1, switch_dists2);
            cerr << "(edge weight " << wt << ")\n";
        }
    }
}

void test_sparse_dynamic_programming(const BaseGraph& graph1,
                                     const BaseGraph& graph2,
                                     const std::vector<match_set_t>& anchors,
                                     bool affine) {
    
    PathMerge chain_merge1(graph1);
    PathMerge chain_merge2(graph2);
    
    TestAnchorer anchorer;
    anchorer.gap_open[0] = 0.25;
    anchorer.gap_extend[0] = 0.25;
    anchorer.gap_open[1] = .4;
    anchorer.gap_extend[1] = 0.15;
    
    std::vector<anchor_t> exhaustive_chain, sparse_chain;
    {
        auto anchors_copy = anchors;
        if (affine) {
            sparse_chain = anchorer.sparse_affine_chain_dp(anchors_copy,
                                                           graph1,
                                                           graph2,
                                                           chain_merge1,
                                                           chain_merge2,
                                                           anchorer.gap_open,
                                                           anchorer.gap_extend, anchors_copy.size(), true);
        }
        else {
            sparse_chain = anchorer.sparse_chain_dp(anchors_copy,
                                                    graph1,
                                                    chain_merge1,
                                                    chain_merge2, anchors_copy.size(), true);
        }
    }
    {
        auto anchors_copy = anchors;
        exhaustive_chain = anchorer.exhaustive_chain_dp(anchors_copy,
                                                        graph1, graph2,
                                                        chain_merge1,
                                                        chain_merge2,
                                                        affine, anchors_copy.size());
    }
    
    
    double exhaustive_score = 0.0, sparse_score = 0.0;
    for (auto& link : exhaustive_chain) {
        exhaustive_score += anchorer.anchor_weight(link.count1, link.count2, link.walk1.size());
    }
    for (auto& link : sparse_chain) {
        sparse_score += anchorer.anchor_weight(link.count1, link.count2, link.walk1.size());
    }
    
    std::vector<vector<size_t>> switch_dists1, switch_dists2;
    if (affine) {
        switch_dists1 = anchorer.post_switch_distances(graph1, chain_merge1);
        switch_dists2 = anchorer.post_switch_distances(graph2, chain_merge2);
        
        for (size_t i = 1; i < exhaustive_chain.size(); ++i) {
            exhaustive_score += anchorer.edge_weight(exhaustive_chain[i-1].walk1.back(), exhaustive_chain[i].walk1.front(),
                                                     exhaustive_chain[i-1].walk2.back(), exhaustive_chain[i].walk2.front(),
                                                     chain_merge1, chain_merge2, switch_dists1, switch_dists2);
        }
        for (size_t i = 1; i < sparse_chain.size(); ++i) {
            sparse_score += anchorer.edge_weight(sparse_chain[i-1].walk1.back(), sparse_chain[i].walk1.front(),
                                                 sparse_chain[i-1].walk2.back(), sparse_chain[i].walk2.front(),
                                                 chain_merge1, chain_merge2, switch_dists1, switch_dists2);
        }
    }
    
    if (abs(exhaustive_score - sparse_score) > 1e-6) {
        cerr << "did not find equivalent chains with sparse and exhaustive DP, affine? " << affine << "\n";
        cerr << "anchor sets:\n";
        for (size_t i = 0; i < anchors.size(); ++i) {
            print_anchor_set(anchors, i);
        }
        cerr << "graphs:\n";
        cerr << cpp_representation(graph1, "graph1") << '\n';
        cerr << cpp_representation(graph2, "graph2") << '\n';
        cerr << "exhaustive chain (score " << exhaustive_score << "):\n";
        print_chain(anchorer, exhaustive_chain, affine,
                    chain_merge1, chain_merge2,
                    switch_dists1, switch_dists2);
        cerr << "sparse chain (score " << sparse_score << "):\n";
        print_chain(anchorer, sparse_chain, affine,
                    chain_merge1, chain_merge2,
                    switch_dists1, switch_dists2);
        exit(1);
    }
}

unordered_map<string, size_t> substring_counts(const string& str) {
    unordered_map<string, size_t> counts;
    for (size_t i = 0; i < str.size(); ++i) {
        for (size_t j = i + 1; j <= str.size(); ++j) {
            ++counts[str.substr(i, j - i)];
        }
    }
    return counts;
}

vector<tuple<string, size_t, size_t>> minimal_rare_matches(const string& str1,
                                                           const string& str2,
                                                           size_t max_count) {
    
    
    auto counts1 = substring_counts(str1);
    auto counts2 = substring_counts(str2);
    
    unordered_map<pair<size_t, size_t>, vector<string>> matches_by_count;
    
    for (auto rec : counts1) {
        if (counts2.count(rec.first)) {
            matches_by_count[make_pair(rec.second, counts2[rec.first])].push_back(rec.first);
        }
    }
    
    vector<tuple<string, size_t, size_t>> mrms;
    
    for (auto rec : matches_by_count) {
        if (rec.first.first > max_count || rec.first.second > max_count) {
            continue;
        }
//        cerr << "considering matches with counts " << rec.first.first << " " << rec.first.second << '\n';
        for (size_t i = 0; i < rec.second.size(); ++i) {
            bool contains = false;
//            cerr << "checking sequence " << rec.second[i] << '\n';
            for (size_t j = 0; j < rec.second.size(); ++j) {
                if (i == j) {
                    continue;
                }
                if (rec.second[i].find(rec.second[j]) != string::npos) {
                    contains = true;
//                    cerr << rec.second[i] << " contains equal count " << rec.second[j] << '\n';
                    break;
                }
            }
            if (!contains) {
//                cerr << rec.second[i] << " has no containing string" << endl;
                mrms.emplace_back(rec.second[i], rec.first.first, rec.first.second);
            }
        }
    }
    
    return mrms;
}

string walk_to_sequence(const BaseGraph& graph, vector<uint64_t>& walk) {
    string seq;
    for (auto n : walk) {
        seq.push_back(decode_base(graph.label(n)));
    }
    return seq;
}

void test_minimal_rare_matches(const string& seq1, const string& seq2, size_t max_count) {
    
    BaseGraph graph1 = make_base_graph("seq1", seq1);
    BaseGraph graph2 = make_base_graph("seq2", seq2);
    
    add_sentinels(graph1, '!', '$');
    add_sentinels(graph2, '#', '%');
    
    PathMerge chain_merge1(graph1);
    PathMerge chain_merge2(graph2);
        
    MatchFinder match_finder;
    match_finder.max_count = max_count;
    
    vector<tuple<string, size_t, size_t>> matches;
    for (auto match : match_finder.find_matches(graph1, graph2)) {
        matches.emplace_back(walk_to_sequence(graph1, match.walks1.front()),
                             match.walks1.size(), match.walks2.size());
    }
    
    vector<tuple<string, size_t, size_t>> expected = minimal_rare_matches(seq1, seq2, max_count);
    sort(expected.begin(), expected.end());
    sort(matches.begin(), matches.end());
    if (matches != expected) {
        cerr << "minimal rare matches failed on sequences " << seq1 << " and " << seq2 << " with max count " << max_count << '\n';
        for (auto m : {expected, matches}) {
            if (m == matches) {
                cerr << "obtained:\n";
            }
            else {
                cerr << "expected:\n";
            }
            for (auto mrm : m) {
                cerr << get<0>(mrm) << ' ' << get<1>(mrm) << ' ' << get<2>(mrm) << '\n';
            }
        }
        exit(1);
    }
}


int main(int argc, char* argv[]) {
     
    random_device rd;
    default_random_engine gen(rd());
    
    // heaviest weight path
    {
        TestAnchorer::AnchorGraph graph;
        uint64_t n0 = graph.add_node(0, 0, 0, 1.0);
        uint64_t n1 = graph.add_node(0, 0, 0, 0.5);
        uint64_t n2 = graph.add_node(0, 0, 0, 0.0);
        uint64_t n3 = graph.add_node(0, 0, 0, 0.5);
        uint64_t n4 = graph.add_node(0, 0, 0, 1.0);
        uint64_t n5 = graph.add_node(0, 0, 0, 1.0);
        uint64_t n6 = graph.add_node(0, 0, 0, 0.5);
        uint64_t n7 = graph.add_node(0, 0, 0, 1.0);

        graph.add_edge(n0, n2);
        graph.add_edge(n1, n2);
        graph.add_edge(n2, n3);
        graph.add_edge(n2, n4);
        graph.add_edge(n3, n5);
        graph.add_edge(n4, n5);
        graph.add_edge(n5, n6);
        graph.add_edge(n5, n7);

        vector<uint64_t> expected{n0, n2, n4, n5, n7};
        assert(graph.heaviest_weight_path() == expected);
    }
    // heaviest path with edge weights
    {
        TestAnchorer::AnchorGraph graph;
        uint64_t n0 = graph.add_node(0, 0, 0, 1.0);
        uint64_t n1 = graph.add_node(0, 0, 0, 1.5);
        uint64_t n2 = graph.add_node(0, 0, 0, 1.0);

        graph.add_edge(n0, n1, -1.0);
        graph.add_edge(n0, n2, 0.0);
        graph.add_edge(n1, n2, -1.0);

        vector<uint64_t> expected{n0, n2};
        assert(graph.heaviest_weight_path() == expected);
    }


    // minimal rare matches tests
    {
        string seq1 = "AAGAG";
        string seq2 = "AAGAG";
        test_minimal_rare_matches(seq1, seq2, 1);
    }
    {
        string seq1 = "GCGACACGACNNG";
        string seq2 = "ATTCGACGCGACA";
        test_minimal_rare_matches(seq1, seq2, 3);
    }
    for (int len : {5, 10, 20, 40}) {

        for (int rep = 0; rep < 5; ++rep) {

            string seq1 = random_sequence(len, gen);
            string seq2 = random_sequence(len, gen);
            string seq3 = mutate_sequence(seq1, 0.1, 0.1, gen);

            for (int max_count : {1, 2, 3, 4, 5}) {

                test_minimal_rare_matches(seq1, seq2, max_count);
                test_minimal_rare_matches(seq1, seq3, max_count);
            }
        }
    }

    {
        BaseGraph graph1;
        for (auto c : std::string("CCAAACCAACCGAACCCAACCCAACCCACT")) {
            graph1.add_node(c);
        }
        
        std::vector<std::pair<int, int>> graph1_edges{
            {0, 1},
            {0, 27},
            {1, 2},
            {2, 3},
            {3, 4},
            {4, 5},
            {5, 6},
            {6, 7},
            {7, 8},
            {7, 11},
            {8, 9},
            {8, 25},
            {8, 28},
            {8, 14},
            {9, 10},
            {9, 11},
            {9, 26},
            {10, 11},
            {10, 26},
            {11, 12},
            {12, 13},
            {13, 14},
            {14, 15},
            {14, 19},
            {15, 16},
            {15, 17},
            {16, 17},
            {17, 18},
            {18, 19},
            {19, 20},
            {20, 21},
            {20, 29},
            {21, 22},
            {22, 23},
            {22, 24},
            {23, 24},
            {25, 10},
            {25, 11},
            {25, 26},
            {26, 12},
            {27, 2},
            {27, 4},
            {28, 10},
            {28, 11},
            {28, 26},
            {29, 22}
        };
        
        std::vector<std::vector<int>> graph1_paths{
            {0, 27, 2, 3, 4, 5, 6, 7, 8, 9, 10, 26, 12, 13, 14, 15, 17, 18, 19, 20, 29, 22, 23, 24},
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 28, 11, 12, 13, 14, 19, 20, 21, 22, 23, 24},
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 25, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 24}
        };
        
        for (auto e : graph1_edges) {
            graph1.add_edge(e.first, e.second);
        }
        
        for (size_t i = 0; i < graph1_paths.size(); ++i) {
            auto p = graph1.add_path(std::to_string(i));
            for (auto n : graph1_paths[i]) {
                graph1.extend_path(p, n);
            }
        }
        
        BaseGraph graph2;
        for (auto c : std::string("AAACAAAACAAAACAAAGCAAAACAATACC")) {
            graph2.add_node(c);
        }
        
        std::vector<std::pair<int, int>> graph2_edges{
            {0, 1},
            {0, 28},
            {1, 2},
            {2, 3},
            {3, 4},
            {4, 5},
            {5, 6},
            {6, 7},
            {7, 8},
            {7, 29},
            {8, 9},
            {9, 10},
            {10, 11},
            {11, 12},
            {11, 15},
            {12, 13},
            {13, 14},
            {14, 15},
            {14, 25},
            {15, 16},
            {16, 17},
            {16, 18},
            {16, 26},
            {17, 18},
            {17, 26},
            {18, 19},
            {18, 27},
            {19, 20},
            {20, 21},
            {21, 22},
            {22, 23},
            {23, 24},
            {25, 16},
            {25, 21},
            {26, 19},
            {26, 27},
            {27, 20},
            {28, 2},
            {29, 9}
        };
        
        std::vector<std::vector<int>> graph2_paths{
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 15, 16, 17, 26, 27, 20, 21, 22, 23, 24},
            {0, 28, 2, 3, 4, 5, 6, 7, 29, 9, 10, 11, 12, 13, 14, 25, 16, 18, 19, 20, 21, 22, 23, 24}
        };
        
        for (auto e : graph2_edges) {
            graph2.add_edge(e.first, e.second);
        }
        
        for (size_t i = 0; i < graph2_paths.size(); ++i) {
            auto p = graph2.add_path(std::to_string(i));
            for (auto n : graph2_paths[i]) {
                graph2.extend_path(p, n);
            }
        }
        
        auto anchors = generate_anchor_set(graph1, graph2, 2);
        test_sparse_dynamic_programming(graph1, graph2, anchors, true);
    }
    
    {
        BaseGraph graph1;
        for (auto c : std::string("CATTGTCCTAAGGAAT")) {
            graph1.add_node(c);
        }

        std::vector<std::pair<int, int>> graph1_edges{
            {0, 15},
            {0, 6},
            {0, 14},
            {0, 11},
            {1, 7},
            {1, 10},
            {1, 6},
            {1, 8},
            {1, 9},
            {2, 3},
            {2, 14},
            {2, 15},
            {2, 11},
            {3, 13},
            {3, 9},
            {4, 11},
            {4, 6},
            {5, 11},
            {5, 6},
            {5, 15},
            {5, 10},
            {6, 15},
            {6, 13},
            {7, 15},
            {7, 11},
            {8, 10},
            {8, 12},
            {9, 15},
            {9, 12},
            {9, 11},
            {10, 13},
            {10, 14},
            {11, 13},
            {13, 14},
            {14, 15}
        };

        std::vector<std::vector<int>> graph1_paths{
            {0, 11, 13, 14, 15},
            {1, 7, 11, 13, 14, 15},
            {2, 3, 9, 12},
            {1, 8, 10, 14, 15},
            {4, 6, 15},
            {5, 10, 14, 15}
        };

        for (auto e : graph1_edges) {
            graph1.add_edge(e.first, e.second);
        }

        for (size_t i = 0; i < graph1_paths.size(); ++i) {
            auto p = graph1.add_path(std::to_string(i));
            for (auto n : graph1_paths[i]) {
                graph1.extend_path(p, n);
            }
        }

        BaseGraph graph2;
        for (auto c : std::string("TGGCTGTGGCACGGAG")) {
            graph2.add_node(c);
        }

        std::vector<std::pair<int, int>> graph2_edges{
            {0, 1},
            {1, 11},
            {1, 13},
            {1, 2},
            {2, 14},
            {2, 10},
            {2, 13},
            {2, 5},
            {3, 12},
            {4, 6},
            {4, 12},
            {5, 15},
            {5, 7},
            {5, 8},
            {5, 10},
            {6, 9},
            {6, 8},
            {6, 15},
            {6, 12},
            {7, 9},
            {7, 12},
            {8, 15},
            {8, 9},
            {8, 14},
            {9, 10},
            {9, 13},
            {9, 12},
            {9, 11},
            {10, 14},
            {10, 11},
            {11, 14},
            {11, 13},
            {11, 15},
            {12, 13},
            {13, 15}
        };

        std::vector<std::vector<int>> graph2_paths{
            {0, 1, 11, 15},
            {4, 12, 13, 15},
            {4, 6, 9, 10, 14},
            {3, 12, 13, 15},
            {0, 1, 2, 5, 8, 14},
            {0, 1, 2, 5, 7, 9, 10, 11, 14}
        };

        for (auto e : graph2_edges) {
            graph2.add_edge(e.first, e.second);
        }

        for (size_t i = 0; i < graph2_paths.size(); ++i) {
            auto p = graph2.add_path(std::to_string(i));
            for (auto n : graph2_paths[i]) {
                graph2.extend_path(p, n);
            }
        }

        auto anchors = generate_anchor_set(graph1, graph2, 2);
        test_sparse_dynamic_programming(graph1, graph2, anchors, true);
    }

    {
        BaseGraph graph1;
        for (auto c : std::string("AAGAAGCACCTCAACT")) {
            graph1.add_node(c);
        }

        std::vector<std::pair<int, int>> graph1_edges{
            {0, 6},
            {0, 11},
            {0, 8},
            {0, 3},
            {1, 2},
            {1, 5},
            {1, 9},
            {1, 13},
            {1, 11},
            {2, 4},
            {2, 11},
            {2, 8},
            {3, 6},
            {3, 8},
            {3, 5},
            {4, 13},
            {4, 9},
            {5, 8},
            {5, 15},
            {5, 10},
            {5, 6},
            {6, 10},
            {6, 15},
            {6, 14},
            {6, 13},
            {6, 11},
            {7, 15},
            {7, 13},
            {8, 10},
            {8, 12},
            {9, 15},
            {10, 12},
            {11, 14},
            {11, 13},
            {13, 14}
        };

        std::vector<std::vector<int>> graph1_paths{
            {1, 5, 10, 12},
            {1, 2, 4, 9, 15},
            {0, 3, 8, 12},
            {0, 6, 11, 13, 14},
            {7, 13, 14}
        };

        for (auto e : graph1_edges) {
            graph1.add_edge(e.first, e.second);
        }

        for (size_t i = 0; i < graph1_paths.size(); ++i) {
            auto p = graph1.add_path(std::to_string(i));
            for (auto n : graph1_paths[i]) {
                graph1.extend_path(p, n);
            }
        }

        BaseGraph graph2;
        for (auto c : std::string("TAGGAGCGCTCGTGGA")) {
            graph2.add_node(c);
        }

        std::vector<std::pair<int, int>> graph2_edges{
            {0, 4},
            {0, 10},
            {0, 14},
            {0, 7},
            {0, 12},
            {0, 15},
            {1, 5},
            {1, 4},
            {1, 10},
            {1, 14},
            {1, 13},
            {1, 7},
            {3, 7},
            {3, 6},
            {4, 9},
            {4, 14},
            {4, 11},
            {4, 15},
            {4, 6},
            {5, 6},
            {5, 10},
            {5, 14},
            {5, 8},
            {5, 12},
            {6, 10},
            {6, 14},
            {7, 11},
            {7, 15},
            {9, 13},
            {9, 15},
            {9, 12},
            {10, 13},
            {10, 14},
            {11, 12},
            {12, 15}
        };

        std::vector<std::vector<int>> graph2_paths{
            {1, 5, 14},
            {0, 4, 6, 10, 13},
            {3, 7, 15},
            {2},
            {1, 4, 11, 12, 15},
            {0, 4, 9, 12, 15},
            {1, 5, 8}
        };

        for (auto e : graph2_edges) {
            graph2.add_edge(e.first, e.second);
        }

        for (size_t i = 0; i < graph2_paths.size(); ++i) {
            auto p = graph2.add_path(std::to_string(i));
            for (auto n : graph2_paths[i]) {
                graph2.extend_path(p, n);
            }
        }

        auto anchors = generate_anchor_set(graph1, graph2, 2);
        test_sparse_dynamic_programming(graph1, graph2, anchors, true);
    }

    {
        BaseGraph graph1;
        for (auto c : std::string("ATTCGGGCAC")) {
            graph1.add_node(c);
        }

        std::vector<std::pair<int, int>> graph1_edges{
            {0, 6},
            {0, 1},
            {0, 3},
            {0, 8},
            {0, 4},
            {1, 5},
            {1, 2},
            {1, 7},
            {1, 8},
            {1, 4},
            {2, 9},
            {3, 6},
            {3, 9},
            {4, 9},
            {5, 6},
            {6, 7},
            {6, 8},
            {8, 9}
        };

        std::vector<std::vector<int>> graph1_paths{
            {0, 6, 7},
            {0, 1, 8, 9},
            {0, 3, 9},
            {0, 1, 5, 6, 8, 9},
            {0, 1, 2, 9},
            {0, 1, 4, 9}
        };

        for (auto e : graph1_edges) {
            graph1.add_edge(e.first, e.second);
        }

        for (size_t i = 0; i < graph1_paths.size(); ++i) {
            auto p = graph1.add_path(std::to_string(i));
            for (auto n : graph1_paths[i]) {
                graph1.extend_path(p, n);
            }
        }

        BaseGraph graph2;
        for (auto c : std::string("AGCTCAAACA")) {
            graph2.add_node(c);
        }

        std::vector<std::pair<int, int>> graph2_edges{
            {0, 1},
            {0, 2},
            {1, 4},
            {1, 2},
            {1, 5},
            {2, 3},
            {2, 5},
            {3, 7},
            {3, 9},
            {4, 5},
            {4, 8},
            {5, 7},
            {5, 8},
            {5, 9},
            {6, 8},
            {6, 9},
            {7, 9},
            {8, 9}
        };

        std::vector<std::vector<int>> graph2_paths{
            {0, 2, 3, 7, 9},
            {0, 1, 5, 8, 9},
            {6, 9},
            {0, 1, 4, 5, 8, 9}
        };

        for (auto e : graph2_edges) {
            graph2.add_edge(e.first, e.second);
        }

        for (size_t i = 0; i < graph2_paths.size(); ++i) {
            auto p = graph2.add_path(std::to_string(i));
            for (auto n : graph2_paths[i]) {
                graph2.extend_path(p, n);
            }
        }

        auto anchors = generate_anchor_set(graph1, graph2, 2);
        test_sparse_dynamic_programming(graph1, graph2, anchors, true);
    }

    {
        BaseGraph graph1;
        for (auto c : std::string("GGCAGTTAGCCATCGC")) {
            graph1.add_node(c);
        }

        std::vector<std::pair<int, int>> graph1_edges{
            {0, 2},
            {0, 8},
            {0, 6},
            {0, 4},
            {1, 10},
            {1, 9},
            {1, 13},
            {1, 7},
            {1, 4},
            {2, 8},
            {2, 7},
            {2, 6},
            {2, 13},
            {3, 8},
            {3, 4},
            {3, 11},
            {3, 9},
            {4, 6},
            {4, 9},
            {4, 11},
            {5, 9},
            {5, 7},
            {6, 14},
            {6, 10},
            {7, 8},
            {7, 10},
            {8, 15},
            {8, 13},
            {8, 14},
            {9, 10},
            {9, 11},
            {10, 11},
            {12, 15},
            {12, 13},
            {13, 15}
        };

        std::vector<std::vector<int>> graph1_paths{
            {1, 7, 10, 11},
            {0, 4, 9, 11},
            {0, 2, 13, 15},
            {3, 8, 14},
            {1, 4, 6, 10, 11},
            {12, 13, 15},
            {5, 9, 10, 11}
        };

        for (auto e : graph1_edges) {
            graph1.add_edge(e.first, e.second);
        }

        for (size_t i = 0; i < graph1_paths.size(); ++i) {
            auto p = graph1.add_path(std::to_string(i));
            for (auto n : graph1_paths[i]) {
                graph1.extend_path(p, n);
            }
        }
        BaseGraph graph2;
        for (auto c : std::string("ACTACATCGGCCGCCG")) {
            graph2.add_node(c);
        }

        std::vector<std::pair<int, int>> graph2_edges{
            {0, 13},
            {0, 15},
            {0, 11},
            {0, 3},
            {1, 8},
            {1, 12},
            {1, 3},
            {1, 7},
            {1, 4},
            {1, 13},
            {2, 3},
            {2, 6},
            {2, 10},
            {3, 15},
            {3, 8},
            {3, 12},
            {4, 13},
            {4, 14},
            {4, 8},
            {5, 6},
            {5, 12},
            {5, 15},
            {5, 11},
            {6, 12},
            {6, 8},
            {6, 7},
            {6, 14},
            {7, 15},
            {8, 14},
            {9, 10},
            {10, 12},
            {11, 14},
            {11, 12},
            {12, 13},
            {14, 15}
        };

        std::vector<std::vector<int>> graph2_paths{
            {5, 12, 13},
            {2, 6, 7, 15},
            {0, 3, 8, 14, 15},
            {9, 10, 12, 13},
            {1, 4, 14, 15},
            {0, 11, 14, 15}
        };

        for (auto e : graph2_edges) {
            graph2.add_edge(e.first, e.second);
        }

        for (size_t i = 0; i < graph2_paths.size(); ++i) {
            auto p = graph2.add_path(std::to_string(i));
            for (auto n : graph2_paths[i]) {
                graph2.extend_path(p, n);
            }
        }

        auto anchors = generate_anchor_set(graph1, graph2, 2);
        test_sparse_dynamic_programming(graph1, graph2, anchors, true);
    }

    // sparse chaining tests
    {
        // baby 2-bubble graphs

        BaseGraph graph1;
        uint64_t n10 = graph1.add_node('A');
        uint64_t n11 = graph1.add_node('A');
        uint64_t n12 = graph1.add_node('G');
        uint64_t n13 = graph1.add_node('T');
        uint64_t n14 = graph1.add_node('G');
        uint64_t n15 = graph1.add_node('C');
        uint64_t n16 = graph1.add_node('C');

        graph1.add_edge(n10, n11);
        graph1.add_edge(n10, n12);
        graph1.add_edge(n11, n13);
        graph1.add_edge(n12, n13);
        graph1.add_edge(n13, n14);
        graph1.add_edge(n13, n15);
        graph1.add_edge(n14, n16);
        graph1.add_edge(n15, n16);

        uint64_t p10 = graph1.add_path("1");
        graph1.extend_path(p10, n10);
        graph1.extend_path(p10, n11);
        graph1.extend_path(p10, n13);
        graph1.extend_path(p10, n14);
        graph1.extend_path(p10, n16);
        uint64_t p11 = graph1.add_path("2");
        graph1.extend_path(p11, n10);
        graph1.extend_path(p11, n12);
        graph1.extend_path(p11, n13);
        graph1.extend_path(p11, n15);
        graph1.extend_path(p11, n16);

        BaseGraph graph2;
        uint64_t n20 = graph2.add_node('G');
        uint64_t n21 = graph2.add_node('A');
        uint64_t n22 = graph2.add_node('G');
        uint64_t n23 = graph2.add_node('A');
        uint64_t n24 = graph2.add_node('T');
        uint64_t n25 = graph2.add_node('C');
        uint64_t n26 = graph2.add_node('C');

        graph2.add_edge(n20, n21);
        graph2.add_edge(n20, n22);
        graph2.add_edge(n21, n23);
        graph2.add_edge(n22, n23);
        graph2.add_edge(n23, n24);
        graph2.add_edge(n23, n25);
        graph2.add_edge(n24, n26);
        graph2.add_edge(n25, n26);

        uint64_t p20 = graph2.add_path("3");
        graph2.extend_path(p20, n20);
        graph2.extend_path(p20, n21);
        graph2.extend_path(p20, n23);
        graph2.extend_path(p20, n25);
        graph2.extend_path(p20, n26);
        uint64_t p21 = graph2.add_path("4");
        graph2.extend_path(p21, n20);
        graph2.extend_path(p21, n22);
        graph2.extend_path(p21, n23);
        graph2.extend_path(p21, n24);
        graph2.extend_path(p21, n26);

        {
            vector<match_set_t> anchors(2);
            anchors[0].walks1.emplace_back(vector<uint64_t>{n10, n11});
            anchors[0].walks2.emplace_back(vector<uint64_t>{n21, n23});
            anchors[1].walks1.emplace_back(vector<uint64_t>{n15, n16});
            anchors[1].walks2.emplace_back(vector<uint64_t>{n25, n26});

            for (auto affine : {true, false}) {
                test_sparse_dynamic_programming(graph1, graph2, anchors, affine);
            }
        }

        // these aren't real matches anymore, but the algorithm doesn't pay
        // attention to the sequence, so whatever
        {
            // put them out of order on graph2
            vector<match_set_t> anchors(2);
            anchors[0].walks1.emplace_back(vector<uint64_t>{n10, n11});
            anchors[0].walks2.emplace_back(vector<uint64_t>{n25, n26});
            anchors[1].walks1.emplace_back(vector<uint64_t>{n15, n16});
            anchors[1].walks2.emplace_back(vector<uint64_t>{n21, n23});

            for (auto affine : {true, false}) {
                test_sparse_dynamic_programming(graph1, graph2, anchors, affine);
            }
        }
        {
            // make it have to choose a better score
            vector<match_set_t> anchors(3);
            anchors[0].walks1.emplace_back(vector<uint64_t>{n10, n11});
            anchors[0].walks2.emplace_back(vector<uint64_t>{n20, n22});
            anchors[1].walks1.emplace_back(vector<uint64_t>{n14});
            anchors[1].walks2.emplace_back(vector<uint64_t>{n24});
            anchors[2].walks1.emplace_back(vector<uint64_t>{n15});
            anchors[2].walks2.emplace_back(vector<uint64_t>{n25});
            anchors[2].walks2.emplace_back(vector<uint64_t>{n26});

            for (auto affine : {true, false}) {
                test_sparse_dynamic_programming(graph1, graph2, anchors, affine);
            }
        }
    }

    // sparse affine dynamic programming on simple linear graphs
    {
        BaseGraph graph1 = make_base_graph("1", "AACCGTT");
        BaseGraph graph2 = make_base_graph("2", "AAGCCTT");

        vector<match_set_t> anchors(3);
        anchors[0].walks1.emplace_back(vector<uint64_t>{0, 1});
        anchors[0].walks2.emplace_back(vector<uint64_t>{0, 1});
        anchors[1].walks1.emplace_back(vector<uint64_t>{2, 3});
        anchors[1].walks2.emplace_back(vector<uint64_t>{3, 4});
        anchors[2].walks1.emplace_back(vector<uint64_t>{5, 6});
        anchors[2].walks2.emplace_back(vector<uint64_t>{5, 6});

        TestAnchorer anchorer;
        anchorer.length_scale = false; // keep the match weights simple
        anchorer.gap_open[0] = 1.0;
        anchorer.gap_extend[0] = 1.0;
        anchorer.gap_open[1] = 3.0;
        anchorer.gap_extend[1] = 0.5;

        PathMerge chain_merge1(graph1);
        PathMerge chain_merge2(graph2);

        auto chain = anchorer.sparse_affine_chain_dp(anchors, graph1, graph2, chain_merge1, chain_merge2,
                                                     anchorer.gap_open, anchorer.gap_extend, anchors.size(), true);
        bool correct = (chain.size() == 2);
        correct &= (chain[0].walk1 == vector<uint64_t>{0, 1});
        correct &= (chain[0].walk2 == vector<uint64_t>{0, 1});
        correct &= (chain[1].walk1 == vector<uint64_t>{5, 6});
        correct &= (chain[1].walk2 == vector<uint64_t>{5, 6});
        assert(correct);
    }

    vector<pair<int, int>> graph_sizes{{7, 10}, {10, 18}, {16, 25}, {16, 35}, {20, 80}};
    for (auto size : graph_sizes) {
        for (int rep = 0; rep < 5; ++rep) {
            BaseGraph graph1 = random_graph(size.first, size.second, gen);
            BaseGraph graph2 = random_graph(size.first, size.second, gen);
            add_random_path_cover(graph1, gen);
            add_random_path_cover(graph2, gen);
//            cerr << "graph1:\n";
//            cerr << cpp_representation(graph1, "graph1");
//            cerr << "graph2:\n";
//            cerr << cpp_representation(graph2, "graph2");
            for (int k : {2, 3}) {
                auto anchors = generate_anchor_set(graph1, graph2, k);
                for (auto affine : {true, false}) {
                    test_sparse_dynamic_programming(graph1, graph2, anchors, affine);
                }
            }
        }
    }

    vector<int> challenge_graph_sizes{10, 20, 30};
    for (auto size : challenge_graph_sizes) {
        for (int rep = 0; rep < 5; ++rep) {
            BaseGraph graph1 = random_challenge_graph(size, gen);
            BaseGraph graph2 = random_challenge_graph(size, gen);
            add_random_path_cover(graph1, gen);
            add_random_path_cover(graph2, gen);
            for (int k : {2, 3, 4}) {
                auto anchors = generate_anchor_set(graph1, graph2, k);
                for (auto affine : {true, false}) {
                    test_sparse_dynamic_programming(graph1, graph2, anchors, affine);
                }
            }
        }
    }
    
    cerr << "passed all tests!" << endl;
}
