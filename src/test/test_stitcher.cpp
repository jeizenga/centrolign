#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <utility>
#include <string>

#include "centrolign/stitcher.hpp"
#include "centrolign/graph.hpp"
#include "centrolign/modify_graph.hpp"
#include "centrolign/chain_merge.hpp"
#include "centrolign/anchorer.hpp"
#include "centrolign/match_finder.hpp"
#include "centrolign/test_util.hpp"

using namespace std;
using namespace centrolign;

class TestStitcher : public Stitcher {
public:
    TestStitcher() : Stitcher() {}
    using Stitcher::identify_despecification_partition;
    
};

tuple<int, double, bool> check_despecification_partition(const vector<pair<size_t, size_t>>& partition,
                                                         const std::vector<anchor_t>& anchors,
                                                         int64_t min_indel_fuzz_length, double indel_fuzz_score_proportion) {
    
    double deleted_score = 0;
    bool valid = true;
    for (auto& interval : partition) {
        if (interval.first == 0 || interval.second == anchors.size()) {
            valid = false;
            //std::cerr << "fail on first/last\n";
            break;
        }
        
        int num_sv_indels = 0;
        for (size_t i = interval.first; i <= interval.second; ++i) {
            if (abs(anchors[i].gap_before) >= min_indel_fuzz_length) {
                ++num_sv_indels;
            }
        }
        
        if (num_sv_indels != 1) {
            //std::cerr << "fail on num indels " << num_sv_indels << "\n";
            valid = false;
            break;
        }
        
        double interval_score = 0.0;
        for (size_t i = interval.first; i < interval.second; ++i) {
            interval_score += anchors[i].score;
        }
        
        if (interval_score >= indel_fuzz_score_proportion * (anchors[interval.first - 1].score + anchors[interval.second].score)) {
            //std::cerr << "fail on score lim\n";
            valid = false;
            break;
        }
        
        deleted_score += interval_score;
    }
    
    return tuple<int, double, bool>(partition.size(), deleted_score, valid);
}

vector<pair<size_t, size_t>> brute_force_despecification_partition(const std::vector<anchor_t>& anchors,
                                                                   int64_t min_indel_fuzz_length, double indel_fuzz_score_proportion) {
    
    assert(anchors.size() < 64);
    
    vector<pair<size_t, size_t>> best_partition;
    int64_t best_num_despecified = -1;
    double best_deleted_score = -100000000;
    
    for (uint64_t set = 0; set < (1 << anchors.size()); ++set) {
        
        vector<pair<size_t, size_t>> partition;
        for (size_t i = 0; i < anchors.size(); ++i) {
            if (set & (1 << i)) {
                if (i == 0 || !(set & (1 << (i - 1)))) {
                    partition.emplace_back(i, i + 1);
                }
                else {
                    partition.back().second = i + 1;
                }
            }
        }
        
        int num_despecified;
        double deleted_score;
        bool valid;
        tie(num_despecified, deleted_score, valid) = check_despecification_partition(partition, anchors, min_indel_fuzz_length, indel_fuzz_score_proportion);
        
        if (valid && (num_despecified > best_num_despecified || (num_despecified == best_num_despecified &&
                                                                 deleted_score < best_deleted_score))) {
            best_num_despecified = num_despecified;
            best_deleted_score = deleted_score;
            best_partition = std::move(partition);
        }
    }
    
    return best_partition;
}


void test_despecification_partition(const std::vector<anchor_t>& anchors,
                                    int64_t min_indel_fuzz_length, double indel_fuzz_score_proportion) {
    
    TestStitcher stitcher;
    stitcher.min_indel_fuzz_length = min_indel_fuzz_length;
    stitcher.indel_fuzz_score_proportion = indel_fuzz_score_proportion;
    
    auto got = stitcher.identify_despecification_partition(anchors);
    auto expected = brute_force_despecification_partition(anchors, min_indel_fuzz_length, indel_fuzz_score_proportion);
    
    int num_despecified_got, num_despecified_expected;
    double deleted_score_got, deleted_score_expected;
    bool valid_got, valid_expected;
    
    tie(num_despecified_got, deleted_score_got, valid_got) = check_despecification_partition(got, anchors, min_indel_fuzz_length, indel_fuzz_score_proportion);
    tie(num_despecified_expected, deleted_score_expected, valid_expected) = check_despecification_partition(expected, anchors, min_indel_fuzz_length, indel_fuzz_score_proportion);
    
    
    double eps = 0.0000000001;
    
    if (!valid_expected) {
        cerr << "test failure: invalid test\n";
    }
    else if (make_pair(num_despecified_got, -deleted_score_got) > make_pair(num_despecified_expected, -deleted_score_expected + eps)) {
        cerr << "test failure: nonoptimal expectation " << num_despecified_got << ", " << deleted_score_got << ", " << num_despecified_expected << ", " << deleted_score_expected << "\n";
    }
    else if (!valid_got) {
        cerr << "test failure: invalid partition\n";
    }
    else if (make_pair(num_despecified_got, -deleted_score_got + eps) < make_pair(num_despecified_expected, -deleted_score_expected)) {
        cerr << "test failure: suboptimal partition\n";
    }
    if (!valid_expected || !valid_got || make_pair(num_despecified_got, -deleted_score_got) > make_pair(num_despecified_expected, -deleted_score_expected + eps) ||
        make_pair(num_despecified_got, -deleted_score_got + eps) < make_pair(num_despecified_expected, -deleted_score_expected)) {
        std::cerr << "min length: " << min_indel_fuzz_length << '\n';
        std::cerr << "proportion: " << indel_fuzz_score_proportion << '\n';
        std::cerr << "anchors:\n";
        for (auto a : anchors) {
            std::cerr << a.score << ", " << a.gap_before << '\n';
        }
        std::cerr << "opt partition, num " << num_despecified_expected << ", deleted " << deleted_score_expected << '\n';
        for (auto p : expected) {
            std::cerr << '\t' << p.first << ':' << p.second << '\n';
        }
        std::cerr << "obtained partition, num " << num_despecified_got << ", deleted " << deleted_score_got << '\n';
        for (auto p : got) {
            std::cerr << '\t' << p.first << ':' << p.second << '\n';
        }
        exit(1);
    }
}


int main(int argc, char* argv[]) {
    
    random_device rd;
    default_random_engine gen(rd());
    
    {
        int64_t min_indel_fuzz_length = 88;
        double indel_fuzz_score_proportion = 0.767724;
        
        std::vector<std::pair<double, int>> anchor_data{
            {68.5195, 9},
            {74.3832, 63},
            {16.0234, 46},
            {80.74, 88},
            {98.8417, 63},
            {40.4285, 19},
            {8.71209, 68}
        };
        
        std::vector<anchor_t> anchors(anchor_data.size());
        for (size_t i = 0; i < anchors.size(); ++i) {
            anchors[i].score = anchor_data[i].first;
            anchors[i].gap_before = anchor_data[i].second;
        }
        
        test_despecification_partition(anchors, min_indel_fuzz_length, indel_fuzz_score_proportion);
    }
    
    {
        int64_t min_indel_fuzz_length = 52;
        double indel_fuzz_score_proportion = 0.767724;
        
        std::vector<std::pair<double, int>> anchor_data{
            {43.4928, 33},
            {12.3185, 71},
            {81.563, 12},
            {93.5945, 13},
            {71.1301, 27},
            {89.5764, 67},
            {30.6563, 37}
        };
        
        std::vector<anchor_t> anchors(anchor_data.size());
        for (size_t i = 0; i < anchors.size(); ++i) {
            anchors[i].score = anchor_data[i].first;
            anchors[i].gap_before = anchor_data[i].second;
        }
        
        test_despecification_partition(anchors, min_indel_fuzz_length, indel_fuzz_score_proportion);
    }
    
    {
        int64_t min_indel_fuzz_length = 34;
        double indel_fuzz_score_proportion = 0.481206;
        
        std::vector<std::pair<double, int>> anchor_data{
            {68.5026, 23},
            {62.0347, 72},
            {4.88778, 63},
            {53.3706, 9},
            {95.8233, 30},
            {6.86524, 100}
        };
        
        std::vector<anchor_t> anchors(anchor_data.size());
        for (size_t i = 0; i < anchors.size(); ++i) {
            anchors[i].score = anchor_data[i].first;
            anchors[i].gap_before = anchor_data[i].second;
        }
        
        test_despecification_partition(anchors, min_indel_fuzz_length, indel_fuzz_score_proportion);
    }
    
    {
        int64_t min_indel_fuzz_length = 15;
        double indel_fuzz_score_proportion = 0.447142;
        
        std::vector<std::pair<double, int>> anchor_data{
            {39.3462, 88},
            {50.8966, 65},
            {21.4857, 9},
            {66.9799, 1},
            {93.1467, 9}
        };
        
        std::vector<anchor_t> anchors(anchor_data.size());
        for (size_t i = 0; i < anchors.size(); ++i) {
            anchors[i].score = anchor_data[i].first;
            anchors[i].gap_before = anchor_data[i].second;
        }
        
        test_despecification_partition(anchors, min_indel_fuzz_length, indel_fuzz_score_proportion);
    }
    
    {
        int64_t min_indel_fuzz_length = 3;
        double indel_fuzz_score_proportion = 0.473987;
        
        std::vector<std::pair<double, int>> anchor_data{
            {61.0955, 3},
            {18.2809, 60},
            {92.3824, 97},
            {56.9985, 55},
            {49.8453, 98}
        };
        
        std::vector<anchor_t> anchors(anchor_data.size());
        for (size_t i = 0; i < anchors.size(); ++i) {
            anchors[i].score = anchor_data[i].first;
            anchors[i].gap_before = anchor_data[i].second;
        }
        
        test_despecification_partition(anchors, min_indel_fuzz_length, indel_fuzz_score_proportion);
    }
    
    {
        int64_t min_indel_fuzz_length = 50;
        double indel_fuzz_score_proportion = 0.05;
        
        std::vector<anchor_t> anchors(3);
        anchors[0].score = 1.0;
        anchors[0].gap_after = 100;
        anchors[1].gap_before = 100;
        anchors[1].score = 0.01;
        anchors[1].gap_after = 0;
        anchors[2].gap_before = 0;
        anchors[2].score = 1.0;
        
        test_despecification_partition(anchors, min_indel_fuzz_length, indel_fuzz_score_proportion);
    }
    
    {
        uniform_int_distribution<int> num_anchors_distr(1, 15);
        uniform_int_distribution<int> indel_len_distr(0, 100);
        uniform_real_distribution<double> score_distr(0.0, 100.0);
        
        uniform_int_distribution<int> indel_min_length_distr(1, 100);
        uniform_real_distribution<double> prop_distr(0.0, 1.0);
        
        
        size_t num_reps = 100;
        for (size_t r = 0; r < num_reps; ++r) {
            size_t n = num_anchors_distr(gen);
            vector<anchor_t> anchors(n);
            for (size_t i = 0; i < n; ++i) {
                anchors[i].score = score_distr(gen);
                anchors[i].gap_before = indel_len_distr(gen);
            }
            int64_t min_indel_fuzz_length = indel_min_length_distr(gen);
            double indel_fuzz_score_proportion = prop_distr(gen);
            
            test_despecification_partition(anchors, min_indel_fuzz_length, indel_fuzz_score_proportion);
        }
    }
    
    Stitcher stitcher;
    
    BaseGraph graph1, graph2;
    for (auto c : string("ACCAGTCGTTGA")) {
        graph1.add_node(c);
    }
    for (auto c : string("GATCGTGAACTATGC")) {
        graph2.add_node(c);
    }
    vector<pair<int, int>> edges1{
        {0, 1},
        {0, 2},
        {1, 3},
        {2, 3},
        {3, 4},
        {4, 5},
        {5, 6},
        {5, 7},
        {6, 8},
        {7, 8},
        {8, 9},
        {8, 10},
        {9, 10},
        {10, 11}
    };
    vector<pair<int, int>> edges2{
        {0, 1},
        {1, 2},
        {1, 3},
        {2, 3},
        {3, 4},
        {4, 5},
        {4, 6},
        {5, 7},
        {6, 7},
        {7, 8},
        {7, 9},
        {8, 10},
        {9, 10},
        {10, 11},
        {10, 12},
        {11, 13},
        {12, 13},
        {13, 14}
    };
    for (auto e : edges1) {
        graph1.add_edge(e.first, e.second);
    }
    for (auto e : edges2) {
        graph2.add_edge(e.first, e.second);
    }
    vector<vector<int>> paths1{
        {0, 1, 3, 4, 5, 6, 8,    10, 11},
        {0, 2, 3, 4, 5, 7, 8, 9, 10, 11}
    };
    vector<vector<int>> paths2{
        {0, 1,    3, 4, 5, 7, 8, 10, 11, 13, 14},
        {0, 1, 2, 3, 4, 6, 7, 9, 10, 12, 13, 14}
    };
    int n = 0;
    for (auto p : paths1) {
        auto i = graph1.add_path(to_string(n++));
        for (auto j : p) {
            graph1.extend_path(i, j);
        }
    }
    for (auto p : paths2) {
        auto i = graph2.add_path(to_string(n++));
        for (auto j : p) {
            graph2.extend_path(i, j);
        }
    }
    
    auto tableau1 = add_sentinels(graph1, '^', '$');
    auto tableau2 = add_sentinels(graph2, '^', '$');
    
    ChainMerge chain_merge1(graph1, tableau1);
    ChainMerge chain_merge2(graph2, tableau2);
    
    vector<vector<pair<int, int>>> anchors{
        {{0, 1}, {1, 3}},
        {{4, 4}, {5, 5}},
        {{9, 12}, {10, 13}}
    };
    
    vector<anchor_t> anchor_chain;
    for (auto a : anchors) {
        anchor_chain.emplace_back();
        anchor_chain.back().count1 = 1;
        anchor_chain.back().count2 = 1;
        for (auto b : a) {
            anchor_chain.back().walk1.push_back(b.first);
            anchor_chain.back().walk2.push_back(b.second);
        }
    }
    
    Alignment expected;
    expected.emplace_back(AlignedPair::gap, 0);
    expected.emplace_back(0, 1);
    expected.emplace_back(1, 3);
    expected.emplace_back(3, AlignedPair::gap);
    expected.emplace_back(4, 4);
    expected.emplace_back(5, 5);
    expected.emplace_back(AlignedPair::gap, 7);
    expected.emplace_back(6, 9);
    expected.emplace_back(8, 10);
    expected.emplace_back(9, 12);
    expected.emplace_back(10, 13);
    expected.emplace_back(11, 14);
    
    vector<vector<anchor_t>> segments(1, anchor_chain);
    
    auto got = stitcher.stitch(segments,
                               graph1, graph2,
                               tableau1, tableau2,
                               chain_merge1, chain_merge2);
    
    check_alignment(got, expected);
    
    
    cerr << "passed all tests!" << endl;
    return 0;
}
