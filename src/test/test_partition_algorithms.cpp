#include <cstdio>
#include <cstdlib>
#include <cstdlib>
#include <cassert>

#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <utility>
#include <sstream>

#include "centrolign/test_util.hpp"

#include "centrolign/average_constrained_partition.hpp"

using namespace std;
using namespace centrolign;

pair<int, bool> check_partition(const vector<pair<size_t, size_t>>& partition,
                                vector<pair<int, int>>& data, int min_avg) {
    int total_score = 0;
    bool valid = true;
    for (auto& interval : partition) {
        int interval_score = 0;
        int interval_weight = 0;
        for (size_t i = interval.first; i < interval.second; ++i) {
            interval_score += data[i].first;
            interval_weight += data[i].second;
        }
        
        if (interval_score < min_avg * interval_weight) {
            valid = false;
            break;
        }
        
        total_score += interval_score;
    }
    
    return make_pair(total_score, valid);
}

vector<pair<size_t, size_t>> brute_force_average_constrained_partition(vector<pair<int, int>>& data, int min_avg) {
    
    assert(data.size() < 64);
    
    vector<pair<size_t, size_t>> best_partition;
    int best_score = -100000000;
    
    for (uint64_t set = 0; set < (1 << data.size()); ++set) {
        
        vector<pair<size_t, size_t>> partition;
        for (size_t i = 0; i < data.size(); ++i) {
            if (set & (1 << i)) {
                if (i == 0 || !(set & (1 << (i - 1)))) {
                    partition.emplace_back(i, i + 1);
                }
                else {
                    partition.back().second = i + 1;
                }
            }
        }
        
        int total_score;
        bool valid;
        tie(total_score, valid) = check_partition(partition, data, min_avg);
        
        if (valid && total_score > best_score) {
            best_score = total_score;
            best_partition = move(partition);
        }
    }
    
    return best_partition;
}

void do_test(vector<pair<int, int>>& data, int min_avg) {
    
    auto got = average_constrained_partition(data, min_avg);
    auto expected = brute_force_average_constrained_partition(data, min_avg);
    
    int score_got, score_expected;
    bool valid_got, valid_expected;
    
    tie(score_got, valid_got) = check_partition(got, data, min_avg);
    tie(score_expected, valid_expected) = check_partition(expected, data, min_avg);
    
    // check for bugs in the test
    assert(valid_expected);
    assert(score_got <= score_expected);
    
    if (!valid_got) {
        cerr << "test failure: invalid partition\n";
    }
    if (score_got != score_expected) {
        cerr << "test failure: suboptimal partition\n";
    }
    if (!valid_got || score_got != score_expected) {
        std::cerr << "min average: " << min_avg << '\n';
        std::cerr << "data:\n";
        for (auto d : data) {
            std::cerr << "\t{" << d.first << "," << d.second << "},\n";
        }
        std::cerr << "opt partition, score " << score_expected << '\n';
        for (auto p : expected) {
            std::cerr << '\t' << p.first << ':' << p.second << '\n';
        }
        std::cerr << "obtained partition, score " << score_got << '\n';
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
        vector<pair<int, int>> data{
            {-1, 1},
            {1, 1},
            {1, 1},
            {-3, 1},
            {2, 1},
            {1, 1},
            {0, 1}
        };
        
        do_test(data, 1);
        do_test(data, 2);
        do_test(data, 3);
    }
    
    uniform_int_distribution<int> data_size_distr(5, 13);
    uniform_int_distribution<int> value_distr(-2, 6);
    uniform_int_distribution<int> weight_distr(1, 2);
    uniform_int_distribution<int> min_avg_distr(0, 2);
    
    int num_reps = 50;
    for (int rep = 0; rep < num_reps; ++rep) {
        
        vector<pair<int, int>> data(data_size_distr(gen));
        for (int i = 0; i < data.size(); ++i) {
            data[i] = make_pair(value_distr(gen), weight_distr(gen));
        }
        int min_avg = min_avg_distr(gen);
        do_test(data, min_avg);
    }
    
    cerr << "passed all tests!" << endl;
    return 0;
}
