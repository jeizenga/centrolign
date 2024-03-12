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

#include "centrolign/partitioner.hpp"

using namespace std;
using namespace centrolign;

class TestPartitioner : public Partitioner {
public:
    using Partitioner::Partitioner;
    using Partitioner::average_constrained_partition;
    using Partitioner::maximum_weight_partition;
    using Partitioner::window_average_constrained_partition;
};

pair<double, bool> check_average_constrained_partition(const vector<pair<size_t, size_t>>& partition,
                                                       vector<pair<double, double>>& data, double min_avg, double penalty) {
    double total_score = 0;
    bool valid = true;
    for (auto& interval : partition) {
        double interval_score = 0;
        double interval_weight = 0;
        for (size_t i = interval.first; i < interval.second; ++i) {
            interval_score += data[i].first;
            interval_weight += data[i].second;
        }
        
        if (interval_score < min_avg * interval_weight) {
            valid = false;
            break;
        }
        
        total_score += interval_score - penalty;
    }
    
    return make_pair(total_score, valid);
}

pair<double, bool> check_window_average_constrained_partition(const vector<pair<size_t, size_t>>& partition,
                                                              const vector<pair<double, double>>& data,
                                                              double window, double min_avg, double penalty) {
    
    
    double total_score = 0;
    bool valid = true;
    for (auto& interval : partition) {
        double interval_score = 0;
        double interval_weight = 0;
        for (size_t i = interval.first; i < interval.second; ++i) {
            interval_score += data[i].first;
            interval_weight += data[i].second;
        }
        
        if (interval_weight <= window && interval_score < min_avg * interval_weight) {
            valid = false;
            break;
        }
        
        if (interval_weight > window) {
            // the left adjusted windows
            for (size_t i = interval.first; i < interval.second; ++i) {
                size_t j = i;
                double window_score = 0.0;
                double window_weight = 0.0;
                while (j < interval.second && window_weight < window) {
                    window_score += data[j].first;
                    window_weight += data[j].second;
                    ++j;
                }
                
                if (window_weight < window) {
                    break;
                }
                
                double overhang = (window  - (window_weight - data[j - 1].second)) / data[j - 1].second;
                double weighted_score = window_score - data[j - 1].first + overhang * data[j - 1].first;
                if (weighted_score / window < min_avg) {
                    valid = false;
                    break;
                }
            }
            
            if (!valid) {
                break;
            }
            
            // right-adjusted windows
            for (int64_t i = interval.second - 1; i >= (int64_t) interval.first; --i) {
                int64_t j = i;
                double window_score = 0.0;
                double window_weight = 0.0;
                while (j >= (int64_t) interval.first && window_weight < window) {
                    window_score += data[j].first;
                    window_weight += data[j].second;
                    --j;
                }
                
                if (window_weight < window) {
                    break;
                }
                
                double overhang = (window  - (window_weight - data[j + 1].second)) / data[j + 1].second;
                double weighted_score = window_score - data[j + 1].first + overhang * data[j + 1].first;
                if (weighted_score / window < min_avg) {
                    valid = false;
                    break;
                }
            }
            
            if (!valid) {
                break;
            }
        }
        
        total_score += interval_score - penalty;
    }
    
    return make_pair(total_score, valid);
}

vector<pair<size_t, size_t>> brute_force_average_constrained_partition(vector<pair<double, double>>& data, double min_avg, double penalty) {
    
    assert(data.size() < 64);
    
    vector<pair<size_t, size_t>> best_partition;
    double best_score = -100000000;
    
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
        
        double total_score;
        bool valid;
        tie(total_score, valid) = check_average_constrained_partition(partition, data, min_avg, penalty);
        
        if (valid && total_score > best_score) {
            best_score = total_score;
            best_partition = move(partition);
        }
    }
    
    return best_partition;
}


vector<pair<size_t, size_t>> brute_force_window_average_constrained_partition(const vector<pair<double, double>>& data,
                                                                              double window, double min_avg, double penalty) {
    
    assert(data.size() < 64);
    
    vector<pair<size_t, size_t>> best_partition;
    double best_score = -100000000;
    
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
        
        double total_score;
        bool valid;
        tie(total_score, valid) = check_window_average_constrained_partition(partition, data, window, min_avg, penalty);
        
//        cerr << "valid " << valid << ", score " << total_score << '\n';
//        for (const auto& p : partition) {
//            cerr << '\t' << p.first << ":" << p.second << '\n';
//        }
        
        if (valid && total_score > best_score) {
            best_score = total_score;
            best_partition = move(partition);
        }
    }
    
    return best_partition;
}

int score_max_weight_partition(vector<pair<size_t, size_t>>& partition, vector<double>& data, double penalty) {
    
    double score = 0;
    for (auto p : partition) {
        score -= penalty;
        for (auto i = p.first; i < p.second; ++i) {
            score += data[i];
        }
    }
    
    return score;
}

vector<pair<size_t, size_t>> brute_force_maximum_weight_partition(vector<double>& data, double penalty) {
    
    assert(data.size() < 64);
    
    vector<pair<size_t, size_t>> best_partition;
    double best_score = -100000000;
    
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
        
        double score = score_max_weight_partition(partition, data, penalty);
        
        if (score > best_score) {
            best_score = score;
            best_partition = move(partition);
        }
    }
    
    return best_partition;
}


void test_window_average_constrained_partition(vector<pair<double, double>>& data, double window, double min_avg, double penalty) {
    
    ScoreFunction score_function;
    score_function.score_scale = 1.0;
    TestPartitioner partitioner(score_function);
    partitioner.window_length = window;
    partitioner.minimum_segment_score = penalty;
    partitioner.minimum_segment_average = min_avg;
    
    auto got = partitioner.window_average_constrained_partition(data);
    auto expected = brute_force_window_average_constrained_partition(data, window, min_avg, penalty);
    
    double score_got, score_expected;
    bool valid_got, valid_expected;
    
    tie(score_got, valid_got) = check_window_average_constrained_partition(got, data, window, min_avg, penalty);
    tie(score_expected, valid_expected) = check_window_average_constrained_partition(expected, data, window, min_avg, penalty);
    
    // check for bugs in the test
    assert(valid_expected);
    assert(score_got <= score_expected);
    
    if (!valid_got) {
        cerr << "test failure: invalid partition\n";
    }
    else if (score_got != score_expected) {
        cerr << "test failure: suboptimal partition\n";
    }
    if (!valid_got || score_got != score_expected) {
        std::cerr << "window: " << window << '\n';
        std::cerr << "min average: " << min_avg << '\n';
        std::cerr << "penalty: " << penalty << '\n';
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

void test_average_constrained_partition(vector<pair<double, double>>& data, double min_avg, double penalty) {
    
    ScoreFunction score_function;
    score_function.score_scale = 1.0;
    TestPartitioner partitioner(score_function);
    partitioner.minimum_segment_score = penalty;
    partitioner.minimum_segment_average = min_avg;
    
    auto got = partitioner.average_constrained_partition(data);
    auto expected = brute_force_average_constrained_partition(data, min_avg, penalty);
    
    double score_got, score_expected;
    bool valid_got, valid_expected;
    
    tie(score_got, valid_got) = check_average_constrained_partition(got, data, min_avg, penalty);
    tie(score_expected, valid_expected) = check_average_constrained_partition(expected, data, min_avg, penalty);
    
    // check for bugs in the test
    assert(valid_expected);
    assert(score_got <= score_expected);
    
    if (!valid_got) {
        cerr << "test failure: invalid partition\n";
    }
    else if (score_got != score_expected) {
        cerr << "test failure: suboptimal partition\n";
    }
    if (!valid_got || score_got != score_expected) {
        std::cerr << "min average: " << min_avg << '\n';
        std::cerr << "penalty: " << penalty << '\n';
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

void test_maximum_weight_partition(vector<double>& data, double penalty) {
    
    ScoreFunction score_function;
    score_function.score_scale = 1.0;
    TestPartitioner partitioner(score_function);
    partitioner.minimum_segment_score = penalty;
    partitioner.minimum_segment_average = 0.0;
    
    auto got = partitioner.maximum_weight_partition(data);
    auto expected = brute_force_maximum_weight_partition(data, penalty);
    
    double score_got = score_max_weight_partition(got, data, penalty);
    double score_expected = score_max_weight_partition(expected, data, penalty);
    
    if (score_got != score_expected) {
        std::cerr << "max weight test failed\n";
        std::cerr << "gap penalty: " << penalty << '\n';
        std::cerr << "data:\n";
        for (auto d : data) {
            std::cerr << "\t" << d << ",\n";
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

void do_tests(vector<pair<double, double>>& data, double window, double min_avg, double penalty) {
    
    vector<double> unweighted;
    for (auto d : data) {
        unweighted.push_back(d.first);
    }
    
    test_window_average_constrained_partition(data, window, min_avg, penalty);
    test_average_constrained_partition(data, min_avg, penalty);
    test_maximum_weight_partition(unweighted, penalty);
}

int main(int argc, char* argv[]) {
    
    random_device rd;
    default_random_engine gen(rd());
    
    // simple case to bootstrap
    {
        vector<pair<double, double>> data{
            {-1.0, 1.0},
            {1.0, 1.0},
            {1.0, 1.0},
            {-3.0, 1.0},
            {2.0, 1.0},
            {1.0, 1.0},
            {0.0, 1.0}
        };
        
        do_tests(data, 1, 1, 0);
        do_tests(data, 1, 2, 0);
        do_tests(data, 1, 3, 0);
    }
    
    // failures from randomized testing
    {
        vector<pair<double, double>> data{
            {1,1},
            {0,1},
            {2,3},
            {-1,3},
            {0,1},
            {0,1},
            {-2,2}
        };
        do_tests(data, 2, 1, 2);
    }
    
    {
        vector<pair<double, double>> data{
            {0,2},
            {-1,3},
            {4,2},
            {3,2},
            {1,3},
            {0,2},
            {5,3},
            {1,1},
            {5,2},
            {1,1},
            {3,2},
            {5,2},
            {4,2},
            {-2,3}
        };
        do_tests(data, 3, 2, 1);
    }
    
    {
        vector<pair<double, double>> data{
            {2,1},
            {5,1},
            {0,1},
            {1,1},
            {2,1},
            {6,3},
            {-1,2},
            {2,3},
            {5,3},
            {1,3},
            {1,2},
            {-2,1},
            {5,3}
        };
        do_tests(data, 2, 0, 0);
    }
    
    {
        vector<pair<double, double>> data{
            {0,1},
            {4,1},
            {1,3},
            {0,2},
            {2,3},
            {1,1},
            {0,2},
            {-2,2},
            {5,1},
            {-1,2},
            {4,2}
        };
        do_tests(data, 3, 1, 2);
    }
    
    {
        vector<pair<double, double>> data{
            {4,1},
            {0,1},
            {4,2},
            {6,2},
            {6,3},
            {2,1},
            {4,3},
            {6,2},
            {6,2},
            {0,1},
            {1,1}
        };
        do_tests(data, 3, 1, 1);
    }
    
    {
        vector<pair<double, double>> data{
            {4,1},
            {6,1},
            {5,2},
            {0,2},
            {1,2},
            {2,1},
            {-1,2},
            {5,1},
            {5,3},
            {5,1},
            {0,1},
            {5,2},
            {5,1},
            {1,2}
        };
        do_tests(data, 3, 2, 2);
    }
    
    {
        vector<pair<double, double>> data{
            {4,2},
            {-2,2},
            {-1,1},
            {-1,2},
            {2,3},
            {1,3},
            {-1,3},
            {-2,1},
            {1,1},
            {-2,3},
            {2,3},
            {3,1}
        };
        do_tests(data, 3, 1, 1);
    }
    
    {
        vector<pair<double, double>> data{
            {4,1},
            {6,1},
            {5,2},
            {0,2},
            {1,2},
            {2,1},
            {-1,2},
            {5,1},
            {5,3},
            {5,1},
            {0,1},
            {5,2},
            {5,1},
            {1,2}
        };
        do_tests(data, 3, 2, 2);
    }
    
    {
        vector<pair<double, double>> data{
            {0,1},
            {-2,1},
            {1,1},
            {0,2},
            {5,1},
            {-1,2},
            {3,1},
            {5,2},
            {3,1},
            {3,1},
            {2,2},
            {5,1},
            {5,1},
            {2,1},
            {-1,2}
        };
        do_tests(data, 2, 1, 2);
    }
    
    {
        vector<pair<double, double>> data{
            {2,2},
            {5,2},
            {-1,2},
            {3,1},
            {6,2},
            {-2,1},
            {3,2},
            {3,2}
        };
        do_tests(data, 3, 2, 0);
    }
    
    uniform_int_distribution<int> data_size_distr(5, 15);
    uniform_int_distribution<int> value_distr(-2, 6);
    uniform_int_distribution<int> weight_distr(1, 3);
    uniform_int_distribution<int> window_distr(1, 4);
    uniform_int_distribution<int> min_avg_distr(0, 2);
    uniform_int_distribution<int> penalty_distr(0, 2);
    
    int num_reps = 100;
    for (int rep = 0; rep < num_reps; ++rep) {
        
        vector<pair<double, double>> data(data_size_distr(gen));
        for (int i = 0; i < data.size(); ++i) {
            data[i] = make_pair(value_distr(gen), weight_distr(gen));
        }
        double window = window_distr(gen);
        double min_avg = min_avg_distr(gen);
        double penalty = penalty_distr(gen);
        do_tests(data, window, min_avg, penalty);
    }
    
    cerr << "passed all tests!" << endl;
    return 0;
}
