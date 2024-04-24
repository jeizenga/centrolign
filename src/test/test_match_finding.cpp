#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <unordered_set>
#include <vector>
#include <list>
#include <string>
#include <random>
#include <sstream>
#include <algorithm>

#include "centrolign/score_function.hpp"
#include "centrolign/path_merge.hpp"
#include "centrolign/match_finder.hpp"
#include "centrolign/modify_graph.hpp"
#include "centrolign/graph.hpp"
#include "centrolign/utility.hpp"
#include "centrolign/test_util.hpp"

using namespace std;
using namespace centrolign;

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
        size_t total_count = rec.first.first * rec.first.second;
        if (total_count > max_count || total_count == 0) {
            continue;
        }
        for (size_t i = 0; i < rec.second.size(); ++i) {
            bool contains = false;
            for (size_t j = 0; j < rec.second.size(); ++j) {
                if (i == j) {
                    continue;
                }
                if (rec.second[i].find(rec.second[j]) != string::npos) {
                    contains = true;
                    break;
                }
            }
            if (!contains) {
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
    
    auto tableau1 = add_sentinels(graph1, '!', '$');
    auto tableau2 = add_sentinels(graph2, '#', '%');
    
    PathMerge<> chain_merge1(graph1);
    PathMerge<> chain_merge2(graph2);
    
    ScoreFunction score_function;
    score_function.anchor_score_function = ScoreFunction::InverseCount;
    MatchFinder match_finder(score_function);
    match_finder.max_count = max_count;
    
    vector<tuple<string, size_t, size_t>> expected = minimal_rare_matches(seq1, seq2, max_count);
    
    sort(expected.begin(), expected.end());
    for (bool use_path_esa : {true, false}) {
        
        match_finder.path_matches = use_path_esa;
        
        for (bool use_color_set_size : {true, false}) {
            
            match_finder.use_color_set_size = use_color_set_size;
            
            
            vector<tuple<string, size_t, size_t>> matches;
            for (auto match : match_finder.find_matches(graph1, graph2, tableau1, tableau2)) {
                matches.emplace_back(walk_to_sequence(graph1, match.walks1.front()),
                                     match.walks1.size(), match.walks2.size());
            }
            
            sort(matches.begin(), matches.end());
            if (matches != expected) {
                cerr << "minimal rare matches failed on sequences " << seq1 << " and " << seq2 << " with max count " << max_count << ", path esa? " << use_path_esa << '\n';
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
    }
}

int main(int argc, char* argv[]) {
    
    random_device rd;
    default_random_engine gen(rd());
    
    {
        string seq1 = "AAGAG";
        string seq2 = "AAGAG";
        test_minimal_rare_matches(seq1, seq2, 1);
    }
    {
        string seq1 = "GCGACACGACNNG";
        string seq2 = "ATTCGACGCGACA";
        test_minimal_rare_matches(seq1, seq2, 5);
    }
    {
        string seq1 = "TCAAGATCTC";
        string seq2 = "TGACATTCTC";
        test_minimal_rare_matches(seq1, seq2, 1);
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
