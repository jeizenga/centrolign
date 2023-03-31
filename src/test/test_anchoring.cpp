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
#include "centrolign/random_graph.hpp"

using namespace std;
using namespace centrolign;

class TestAnchorer : public Anchorer {
public:
    using Anchorer::anchor_set_t;
    using Anchorer::AnchorGraph;
    using Anchorer::find_matches;
};

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
        seq.push_back(decode_base(graph.base(n)));
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
    
    // full anchoring test
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
    
    
    cerr << "passed all tests!" << endl;
}
