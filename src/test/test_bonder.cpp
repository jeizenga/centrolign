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

#include "centrolign/bonder.hpp"

using namespace std;
using namespace centrolign;

class TestBonder : public Bonder {
public:
    using Bonder::Bonder;
    using Bonder::longest_partition;
};

std::pair<double, bool> check_longest_partition(const std::vector<std::pair<size_t, size_t>>& partition,
                                                const std::vector<std::tuple<double, double, double>>& shared,
                                                const std::vector<std::tuple<double, double, double>>& intervening,
                                                double min_opt_prop, double min_length) {
    
    double score = 0.0;
    bool valid = true;
    for (auto interval : partition) {
        
        double opt_total = 0.0;
        double sec_total = 0.0;
        double length_total = 0.0;
        for (size_t i = interval.first; i < interval.second; ++i) {
            
            if (i != interval.first) {
                length_total += std::get<0>(intervening[i - 1]);
                opt_total += std::get<1>(intervening[i - 1]);
                sec_total += std::get<2>(intervening[i - 1]);
            }
            
            length_total += std::get<0>(shared[i]);
            opt_total += std::get<1>(shared[i]);
            sec_total += std::get<2>(shared[i]);
        }
        valid = valid && (sec_total > min_opt_prop * opt_total);
        score += length_total - min_length;
    }
    
    return std::make_pair(score, valid);
}

vector<pair<size_t, size_t>> brute_force_longest_partition(const std::vector<std::tuple<double, double, double>>& shared,
                                                           const std::vector<std::tuple<double, double, double>>& intervening,
                                                           double min_opt_prop, double min_length) {
    
    assert(shared.size() < 64);
    
    vector<pair<size_t, size_t>> best_partition;
    double best_score = -100000000;
    
    for (uint64_t set = 0; set < (1 << shared.size()); ++set) {
        
        vector<pair<size_t, size_t>> partition;
        for (size_t i = 0; i < shared.size(); ++i) {
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
        tie(total_score, valid) = check_longest_partition(partition, shared, intervening, min_opt_prop, min_length);
        
        if (valid && total_score > best_score) {
            best_score = total_score;
            best_partition = move(partition);
        }
    }
    
    return best_partition;
}

void test_longest_partition(const std::vector<std::tuple<double, double, double>>& shared,
                            const std::vector<std::tuple<double, double, double>>& intervening,
                            double min_opt_prop, double min_length) {

    TestBonder bonder;
    bonder.min_opt_proportion = min_opt_prop;
    bonder.min_length = min_length;
    
    auto got = bonder.longest_partition(shared, intervening);
    auto expected = brute_force_longest_partition(shared, intervening, min_opt_prop, min_length);
    
    double score_got, score_expected;
    bool valid_got, valid_expected;
    
    tie(score_got, valid_got) = check_longest_partition(got, shared, intervening, min_opt_prop, min_length);
    tie(score_expected, valid_expected) = check_longest_partition(expected, shared, intervening, min_opt_prop, min_length);
    
    if (!valid_got) {
        cerr << "test failure: invalid partition\n";
    }
    else if (score_got != score_expected) {
        cerr << "test failure: suboptimal partition\n";
    }
    else if (!valid_expected || score_got > score_expected) {
        cerr << "test failure: test did not identify optimal\n";
    }
    if (!valid_got || score_got != score_expected || !valid_expected || score_got > score_expected) {
        std::cerr << "min opt prop: " << min_opt_prop << '\n';
        std::cerr << "min length: " << min_length << '\n';
        std::cerr << "shared:\n";
        for (auto d : shared) {
            std::cerr << "\t{" << std::get<0>(d) << "," << std::get<1>(d) << "," << std::get<2>(d) << "},\n";
        }
        std::cerr << "intervening:\n";
        for (auto d : intervening) {
            std::cerr << "\t{" << std::get<0>(d) << "," << std::get<1>(d) << "," << std::get<2>(d) << "},\n";
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
        std::vector<std::tuple<double, double, double>> shared{
            {1.0, 1.0, 0.8},
            {1.0, 1.0, 0.8}
        };
        std::vector<std::tuple<double, double, double>> intervening{
            {1.0, 1.0, 0.0}
        };
        double min_opt_prop = 0.75;
        double min_length = 0.0;
        test_longest_partition(shared, intervening, min_opt_prop, min_length);
    }

    {
        std::vector<std::tuple<double, double, double>> shared{
            {4.13927,0.57112,0.0509141},
            {2.74866,0.853637,0.517178},
            {7.04387,0.240957,0.182393},
            {9.00506,0.702367,0.136439},
            {5.76231,0.493398,0.242135}
        };
        std::vector<std::tuple<double, double, double>> intervening{
            {7.70762,0.670648,0.490583},
            {3.89254,0.308082,0.126494},
            {7.71464,0.085307,0.0389388},
            {5.01921,0.0525329,0.029716}
        };
        double min_opt_prop = 0.484905;
        double min_length = 7.74479;
        test_longest_partition(shared, intervening, min_opt_prop, min_length);
    }
    
    {
        std::vector<std::tuple<double, double, double>> shared{
            {6.23226,0.156539,0.0927762},
            {1.93783,0.137939,0.0723573},
            {7.96105,0.241554,0.0448032},
            {6.20772,0.185558,0.118737},
            {6.46017,0.449736,0.281866}
        };
        std::vector<std::tuple<double, double, double>> intervening{
            {8.08854,0.200853,0.130206},
            {1.37111,0.950428,0.0914957},
            {6.96297,0.0630396,0.0227176},
            {3.31642,0.0347763,0.0188846}
        };
        double min_opt_prop = 0.417389;
        double min_length = 14.8837;
        test_longest_partition(shared, intervening, min_opt_prop, min_length);
    }
    
    {
        std::vector<std::tuple<double, double, double>> shared{
            {9.58639,0.950358,0.355838},
            {7.04428,0.86397,0.369675},
            {1.80685,0.968014,0.847936},
            {1.88018,0.81249,0.124165},
            {6.46014,0.079866,0.0347434},
            {1.22745,0.709431,0.426963},
            {3.15084,0.980257,0.124857},
            {5.35253,0.588144,0.118747},
            {4.7296,0.534231,0.3005},
            {5.82534,0.40199,0.383159}
        };
        std::vector<std::tuple<double, double, double>> intervening{
            {5.26577,0.498366,0.404198},
            {6.79596,0.507775,0.426392},
            {9.18299,0.147504,0.11051},
            {7.94498,0.0238321,0.0125835},
            {6.16545,0.539503,0.522885},
            {5.03586,0.184744,0.00660905},
            {4.08447,0.0327913,0.0232785},
            {2.24852,0.708383,0.505098},
            {7.73638,0.931438,0.638612}
        };
        double min_opt_prop = 0.632286;
        double min_length = 13.3669;
        test_longest_partition(shared, intervening, min_opt_prop, min_length);
    }

    
    uniform_real_distribution<double> min_prop_distr(0.0, 1.0);
    uniform_real_distribution<double> length_distr(1.0, 10.0);
    uniform_real_distribution<double> min_length_distr(0.0, 20.0);
    uniform_real_distribution<double> score_distr(0.0, 1.0);
    uniform_real_distribution<double> sec_multiplier_distr(0.0, 1.1);
    
    vector<int> sizes{5, 10, 15};
    int num_reps = 20;
    for (auto size : sizes) {
        for (int rep = 0; rep < num_reps; ++rep) {
            
            vector<tuple<double, double, double>> shared(size), intervening(size - 1);
            for (size_t i = 0; i < shared.size(); ++i) {
                get<0>(shared[i]) = length_distr(gen);
                get<1>(shared[i]) = score_distr(gen);
                get<2>(shared[i]) = get<1>(shared[i]) * sec_multiplier_distr(gen);
            }
            for (size_t i = 0; i < intervening.size(); ++i) {
                get<0>(intervening[i]) = length_distr(gen);
                get<1>(intervening[i]) = score_distr(gen);
                get<2>(intervening[i]) = get<1>(intervening[i]) * sec_multiplier_distr(gen);
            }
            
            double min_prop = min_prop_distr(gen);
            double min_length = min_length_distr(gen);
            
            test_longest_partition(shared, intervening, min_prop, min_length);
        }
    }
    
    
    cerr << "passed all tests!" << endl;
}
