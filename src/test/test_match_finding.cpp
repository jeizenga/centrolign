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
    
    ScoreFunction score_function;
    score_function.anchor_score_function = ScoreFunction::InverseCount;
    
    vector<tuple<string, size_t, size_t>> expected = minimal_rare_matches(seq1, seq2, max_count);
    
    sort(expected.begin(), expected.end());
    for (bool use_path_esa : {true, false}) {
                
        for (bool use_color_set_size : {true, false}) {
            
            std::vector<match_set_t> raw_matches;
            if (use_path_esa) {
                PathMatchFinder match_finder(score_function);
                match_finder.max_count = max_count;
                match_finder.use_color_set_size = use_color_set_size;
                raw_matches = match_finder.find_matches(graph1, graph2, tableau1, tableau2);
            }
            else {
                
                GESAMatchFinder match_finder(score_function);
                match_finder.max_count = max_count;
                match_finder.use_color_set_size = use_color_set_size;
                raw_matches = match_finder.find_matches(graph1, graph2, tableau1, tableau2);
            }
            
            vector<tuple<string, size_t, size_t>> matches;
            for (auto match : raw_matches) {
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

void test_count_index_equivalence(const BaseGraph& graph1, const BaseGraph& graph2,
                                  const SentinelTableau& tableau1, const SentinelTableau& tableau2,
                                  size_t max_count) {
    
    ScoreFunction score_function;
    score_function.anchor_score_function = ScoreFunction::InverseCount;
    PathMatchFinder match_finder(score_function);
    match_finder.max_count = max_count;
    
    match_finder.use_color_set_size = true;
    auto matches_css = match_finder.find_matches(graph1, graph2, tableau1, tableau2);
    match_finder.use_color_set_size = false;
    auto matches_ruq = match_finder.find_matches(graph1, graph2, tableau1, tableau2);
    
    for (auto match_ptr : {&matches_css, &matches_ruq}) {
        auto& matches = *match_ptr;
        for (auto& match_set : matches) {
            sort(match_set.walks1.begin(), match_set.walks1.end());
            sort(match_set.walks2.begin(), match_set.walks2.end());
        }
        sort(matches.begin(), matches.end(),
             [](const match_set_t& a, const match_set_t& b) {
            return a.walks1 < b.walks1 || (a.walks1 == b.walks1 && a.walks2 < b.walks2);
        });
    }
    
    if (matches_css.size() != matches_ruq.size()) {
        std::cerr << "CSS and RUQ find different numbers of matches\n";
        exit(1);
    }
    for (size_t i = 0; i < matches_ruq.size(); ++i) {
        if (matches_ruq[i].walks1 != matches_css[i].walks1 ||
            matches_ruq[i].walks2 != matches_css[i].walks2) {
            std::cerr << "CSS and RUQ find different match sets\n";
            exit(1);
        }
    }
}

int main(int argc, char* argv[]) {
    
    random_device rd;
    default_random_engine gen(rd());
    
    // cases from automated testing
    {
        BaseGraph graph1;
        for (auto c : std::string("CACCCGCTTT")) {
            graph1.add_node(c);
        }

        std::vector<std::pair<int, int>> graph1_edges{
            {0, 1},
            {0, 7},
            {0, 1},
            {1, 2},
            {1, 2},
            {2, 3},
            {2, 6},
            {2, 8},
            {2, 9},
            {2, 4},
            {3, 4},
            {5, 1},
            {5, 7},
            {6, 4},
            {7, 2},
            {8, 4},
            {8, 4},
            {9, 4}
        };

        std::vector<std::vector<int>> graph1_paths{
            {0, 1, 2, 4},
            {5, 7, 2, 6, 4},
            {0, 7, 2, 3, 4},
            {0, 1, 2, 8, 4},
            {0, 7, 2, 9, 4}
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
        for (auto c : std::string("CGCCCAGACC")) {
            graph2.add_node(c);
        }

        std::vector<std::pair<int, int>> graph2_edges{
            {0, 1},
            {0, 7},
            {0, 8},
            {0, 2},
            {0, 6},
            {0, 3},
            {0, 5},
            {0, 2},
            {1, 2},
            {1, 6},
            {1, 3},
            {1, 5},
            {2, 3},
            {2, 5},
            {3, 4},
            {5, 4},
            {6, 3},
            {6, 5},
            {7, 2},
            {7, 6},
            {7, 3},
            {7, 5},
            {8, 2},
            {8, 6},
            {8, 3},
            {8, 5},
            {8, 4},
            {9, 1},
            {9, 7},
            {9, 8},
            {9, 2},
            {9, 6},
            {9, 3},
            {9, 5},
            {9, 4}
        };

        std::vector<std::vector<int>> graph2_paths{
            {0, 2, 3, 4},
            {9, 7, 5, 4},
            {9, 1, 6, 5, 4},
            {0, 8, 5, 4}
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

        auto tableau1 = add_sentinels(graph1, '!', '$');
        auto tableau2 = add_sentinels(graph2, '#', '%');

        test_count_index_equivalence(graph1, graph2, tableau1, tableau2, 1);
    }
    
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
    
    
    
    vector<int> graph_sizes{10, 20, 30, 50};
    for (auto size : graph_sizes) {
        for (int rep = 0; rep < 5; ++rep) {
            
            BaseGraph graph1 = random_challenge_graph(size, gen);
            BaseGraph graph2 = random_challenge_graph(size, gen);
            add_random_path_cover(graph1, gen);
            add_random_path_cover(graph2, gen);
            
            auto tableau1 = add_sentinels(graph1, '!', '$');
            auto tableau2 = add_sentinels(graph2, '#', '%');
            
            for (int max_count : {1, 2, 3}) {
                test_count_index_equivalence(graph1, graph2, tableau1, tableau2, max_count);
            }
        }
    }
    
    cerr << "passed all tests!" << endl;
}
