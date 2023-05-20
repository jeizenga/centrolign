#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <random>

#include "centrolign/anchorer.hpp"
#include "centrolign/chain_merge.hpp"
#include "centrolign/modify_graph.hpp"
#include "centrolign/graph.hpp"
#include "centrolign/utility.hpp"
#include "centrolign/test_util.hpp"

using namespace std;
using namespace centrolign;

class TestAnchorer : public Anchorer {
public:
    using Anchorer::anchor_set_t;
    using Anchorer::AnchorGraph;
    using Anchorer::find_matches;
    using Anchorer::exhaustive_chain_dp;
    using Anchorer::sparse_chain_dp;
    using Anchorer::anchor_weight;
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

std::vector<TestAnchorer::anchor_set_t> generate_anchor_set(const BaseGraph& graph1,
                                                            const BaseGraph& graph2,
                                                            size_t k) {
    
    auto kmers1 = k_mer_walks(graph1, k);
    auto kmers2 = k_mer_walks(graph2, k);
    
    std::vector<TestAnchorer::anchor_set_t> anchors;
    for (pair<const string, vector<vector<uint64_t>>>& entry : kmers1) {
        if (kmers2.count(entry.first)) {
            anchors.emplace_back();
            anchors.back().walks1 = move(entry.second);
            anchors.back().walks2 = move(kmers2[entry.first]);
        }
    }
    
    return anchors;
}

void print_anchor_set(const vector<TestAnchorer::anchor_set_t>& anchors, size_t i) {
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

void print_chain(const vector<anchor_t>& chain) {
    for (size_t i = 0; i < chain.size(); ++i) {
        auto& link = chain[i];
        cerr << i << " (count1 " << link.count1 << ", count2 " << link.count2 << "):\n";
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
    }
}

void test_sparse_dynamic_programming(const BaseGraph& graph1,
                                     const BaseGraph& graph2,
                                     const std::vector<TestAnchorer::anchor_set_t>& anchors) {
    
    ChainMerge chain_merge1(graph1);
    ChainMerge chain_merge2(graph2);
    
    TestAnchorer anchorer;
    
    std::vector<anchor_t> exhaustive_chain, sparse_chain;
    {
        auto anchors_copy = anchors;
        sparse_chain = anchorer.sparse_chain_dp(anchors_copy,
                                                graph1,
                                                chain_merge1,
                                                chain_merge2);
    }
    {
        auto anchors_copy = anchors;
        exhaustive_chain = anchorer.exhaustive_chain_dp(anchors_copy,
                                                        chain_merge1,
                                                        chain_merge2);
    }
    
    double exhaustive_score = 0.0, sparse_score = 0.0;
    for (auto& link : exhaustive_chain) {
        exhaustive_score += anchorer.anchor_weight(link.count1, link.count2, link.walk1.size());
    }
    for (auto& link : sparse_chain) {
        sparse_score += anchorer.anchor_weight(link.count1, link.count2, link.walk1.size());
    }
    
    if (abs(exhaustive_score - sparse_score) > 1e-6) {
        cerr << "did not find equivalent chains with sparse and exhaustive DP\n";
        cerr << "graph 1:\n";
        print_graph(graph1, cerr);
        cerr << "graph 2:\n";
        print_graph(graph2, cerr);
        cerr << "anchor sets:\n";
        for (size_t i = 0; i < anchors.size(); ++i) {
            print_anchor_set(anchors, i);
        }
        cerr << "exhaustive chain:\n";
        print_chain(exhaustive_chain);
        cerr << "sparse chain:\n";
        print_chain(sparse_chain);
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
    
    ChainMerge chain_merge1(graph1);
    ChainMerge chain_merge2(graph2);
    
    TestAnchorer anchorer;
    anchorer.max_count = max_count;
    
    vector<tuple<string, size_t, size_t>> matches;
    for (auto match : anchorer.find_matches(graph1, graph2)) {
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
            vector<TestAnchorer::anchor_set_t> anchors(2);
            anchors[0].walks1.emplace_back(vector<uint64_t>{n10, n11});
            anchors[0].walks2.emplace_back(vector<uint64_t>{n21, n23});
            anchors[1].walks1.emplace_back(vector<uint64_t>{n15, n16});
            anchors[1].walks2.emplace_back(vector<uint64_t>{n25, n26});
            
            test_sparse_dynamic_programming(graph1, graph2, anchors);
        }
        
        // these aren't real matches anymore, but the algorithm doesn't pay
        // attention to the sequence, so whatever
        {
            // put them out of order on graph2
            vector<TestAnchorer::anchor_set_t> anchors(2);
            anchors[0].walks1.emplace_back(vector<uint64_t>{n10, n11});
            anchors[0].walks2.emplace_back(vector<uint64_t>{n25, n26});
            anchors[1].walks1.emplace_back(vector<uint64_t>{n15, n16});
            anchors[1].walks2.emplace_back(vector<uint64_t>{n21, n23});
            
            test_sparse_dynamic_programming(graph1, graph2, anchors);
        }
        {
            // make it have to choose a better score
            vector<TestAnchorer::anchor_set_t> anchors(3);
            anchors[0].walks1.emplace_back(vector<uint64_t>{n10, n11});
            anchors[0].walks2.emplace_back(vector<uint64_t>{n20, n22});
            anchors[1].walks1.emplace_back(vector<uint64_t>{n14});
            anchors[1].walks2.emplace_back(vector<uint64_t>{n24});
            anchors[2].walks1.emplace_back(vector<uint64_t>{n15});
            anchors[2].walks2.emplace_back(vector<uint64_t>{n25});
            anchors[2].walks2.emplace_back(vector<uint64_t>{n26});
            
            test_sparse_dynamic_programming(graph1, graph2, anchors);
        }
    }
    
    vector<pair<int, int>> graph_sizes{{7, 10}, {16, 25}, {16, 35}, {20, 80}};
    for (auto size : graph_sizes) {
        for (int rep = 0; rep < 5; ++rep) {
            BaseGraph graph1 = random_graph(size.first, size.second, gen);
            BaseGraph graph2 = random_graph(size.first, size.second, gen);
            add_random_path_cover(graph1, gen);
            add_random_path_cover(graph2, gen);
            for (int k : {2, 3}) {
                auto anchors = generate_anchor_set(graph1, graph2, k);
                test_sparse_dynamic_programming(graph1, graph2, anchors);
            }
        }
    }
    
    cerr << "passed all tests!" << endl;
}
