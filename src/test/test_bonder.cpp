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
    using Bonder::longest_windowed_partition;
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

std::pair<double, bool> check_longest_windowed_partition(const std::vector<std::pair<size_t, size_t>>& partition,
                                                         const std::vector<std::tuple<double, double, double>>& shared,
                                                         const std::vector<std::tuple<double, double, double>>& intervening,
                                                         double min_opt_prop, double min_length, double window_length, bool print_failure = false) {
    double score = 0.0;
    bool valid = true;
    for (auto interval : partition) {
        
        double opt_total = 0.0;
        double sec_total = 0.0;
        double length_total = 0.0;
        for (size_t i = interval.first; i < interval.second; ++i) {
            
            int64_t w = i;
            bool from_shared = true;
            double curr_window_length = 0.0;
            double window_opt_total = 0.0;
            double window_sec_total = 0.0;
            while (!(w == interval.second - 1 && !from_shared)) {
                double to_add = from_shared ? std::get<0>(shared[w]) : std::get<0>(intervening[w]);
                if (curr_window_length + to_add >= window_length) {
                    double overhang = (window_length - curr_window_length) / to_add;
                    const auto& final_seg = from_shared ? shared[w] : intervening[w];
                    window_opt_total += overhang * std::get<1>(final_seg);
                    window_sec_total += overhang * std::get<2>(final_seg);
                    bool valid_here = (window_sec_total >= min_opt_prop * window_opt_total);
                    if (!valid_here && print_failure) {
                        std::cerr << "failed rightward from " << i << " shared\n";
                    }
                    valid = valid && valid_here;
                    break;
                }
                curr_window_length += to_add;
                if (from_shared) {
                    window_opt_total += std::get<1>(shared[w]);
                    window_sec_total += std::get<2>(shared[w]);
                }
                else {
                    window_opt_total += std::get<1>(intervening[w]);
                    window_sec_total += std::get<2>(intervening[w]);
                    ++w;
                }
                from_shared = !from_shared;
            }
            
            w = i;
            from_shared = true;
            curr_window_length = 0.0;
            window_opt_total = 0.0;
            window_sec_total = 0.0;
            while (!(w == int64_t(interval.first) - 1)) {
                double to_add = from_shared ? std::get<0>(shared[w]) : std::get<0>(intervening[w]);
                if (curr_window_length + to_add >= window_length) {
                    double overhang = (window_length - curr_window_length) / to_add;
                    const auto& final_seg = from_shared ? shared[w] : intervening[w];
                    window_opt_total += overhang * std::get<1>(final_seg);
                    window_sec_total += overhang * std::get<2>(final_seg);
                    bool valid_here = (window_sec_total >= min_opt_prop * window_opt_total);
                    if (!valid_here && print_failure) {
                        std::cerr << "failed leftward from " << i << " shared\n";
                    }
                    valid = valid && valid_here;
                    break;
                }
                curr_window_length += to_add;
                if (from_shared) {
                    window_opt_total += std::get<1>(shared[w]);
                    window_sec_total += std::get<2>(shared[w]);
                    --w;
                }
                else {
                    window_opt_total += std::get<1>(intervening[w]);
                    window_sec_total += std::get<2>(intervening[w]);
                }
                from_shared = !from_shared;
            }
            
            if (i + 1 != interval.second) {
                
                int64_t w = i;
                from_shared = false;
                curr_window_length = 0.0;
                window_opt_total = 0.0;
                window_sec_total = 0.0;
                while (!(w == interval.second - 1 && !from_shared)) {
                    double to_add = from_shared ? std::get<0>(shared[w]) : std::get<0>(intervening[w]);
                    if (curr_window_length + to_add >= window_length) {
                        double overhang = (window_length - curr_window_length) / to_add;
                        const auto& final_seg = from_shared ? shared[w] : intervening[w];
                        window_opt_total += overhang * std::get<1>(final_seg);
                        window_sec_total += overhang * std::get<2>(final_seg);
                        bool valid_here = (window_sec_total >= min_opt_prop * window_opt_total);
                        if (!valid_here && print_failure) {
                            std::cerr << "failed rightward from " << i << " intervening\n";
                        }
                        valid = valid && valid_here;
                        break;
                    }
                    curr_window_length += to_add;
                    if (from_shared) {
                        window_opt_total += std::get<1>(shared[w]);
                        window_sec_total += std::get<2>(shared[w]);
                    }
                    else {
                        window_opt_total += std::get<1>(intervening[w]);
                        window_sec_total += std::get<2>(intervening[w]);
                        ++w;
                    }
                    from_shared = !from_shared;
                }
                
                w = i;
                from_shared = false;
                curr_window_length = 0.0;
                window_opt_total = 0.0;
                window_sec_total = 0.0;
                while (!(w == int64_t(interval.first) - 1)) {
                    double to_add = from_shared ? std::get<0>(shared[w]) : std::get<0>(intervening[w]);
                    if (curr_window_length + to_add >= window_length) {
                        double overhang = (window_length - curr_window_length) / to_add;
                        const auto& final_seg = from_shared ? shared[w] : intervening[w];
                        window_opt_total += overhang * std::get<1>(final_seg);
                        window_sec_total += overhang * std::get<2>(final_seg);
                        bool valid_here = (window_sec_total >= min_opt_prop * window_opt_total);
                        if (!valid_here && print_failure) {
                            std::cerr << "failed leftward from " << i << " intervening\n";
                        }
                        valid = valid && valid_here;
                        break;
                    }
                    curr_window_length += to_add;
                    if (from_shared) {
                        window_opt_total += std::get<1>(shared[w]);
                        window_sec_total += std::get<2>(shared[w]);
                        --w;
                    }
                    else {
                        window_opt_total += std::get<1>(intervening[w]);
                        window_sec_total += std::get<2>(intervening[w]);
                    }
                    from_shared = !from_shared;
                }
            }
            
            if (i != interval.first) {
                if (std::get<0>(intervening[i - 1]) >= window_length) {
                    // we don't permit entirely intervening windows
                    valid = false;
                    if (print_failure) {
                        std::cerr << "fail on too long intervening between " << (i - 1) << " and " << i << '\n';
                    }
                }
                
                length_total += std::get<0>(intervening[i - 1]);
                opt_total += std::get<1>(intervening[i - 1]);
                sec_total += std::get<2>(intervening[i - 1]);
            }
            
            length_total += std::get<0>(shared[i]);
            opt_total += std::get<1>(shared[i]);
            sec_total += std::get<2>(shared[i]);
        }
        // this condition captures intervals shorter than the window
        if (length_total <= window_length) {
            bool valid_here = (sec_total > min_opt_prop * opt_total);
            if (print_failure && !valid_here) {
                std::cerr << "fail on full prop from interval " << interval.first << ":" << interval.second << '\n';
            }
            valid = valid && valid_here;
        }
        score += length_total - min_length;
    }
    
    return std::make_pair(score, valid);
}

vector<pair<size_t, size_t>> brute_force_longest_windowed_partition(const std::vector<std::tuple<double, double, double>>& shared,
                                                                    const std::vector<std::tuple<double, double, double>>& intervening,
                                                                    double min_opt_prop, double min_length, double window_length) {
    
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
        tie(total_score, valid) = check_longest_windowed_partition(partition, shared, intervening, min_opt_prop, min_length, window_length);
        
        if (valid && total_score > best_score) {
            best_score = total_score;
            best_partition = move(partition);
        }
    }
    
    return best_partition;
}

void apply_test(const std::vector<std::tuple<double, double, double>>& shared,
                const std::vector<std::tuple<double, double, double>>& intervening,
                const vector<pair<size_t, size_t>>& got, const vector<pair<size_t, size_t>>& expected,
                double score_got, double score_expected, bool valid_got, bool valid_expected,
                double min_opt_prop, double min_length, double window_length) {
    
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
        std::cerr << "window length: " << window_length << '\n';
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
    
    
    apply_test(shared, intervening, got, expected, score_got, score_expected, valid_got, valid_expected,
               min_opt_prop, min_length, -1.0);
    
}

void test_longest_windowed_partition(const std::vector<std::tuple<double, double, double>>& shared,
                                      const std::vector<std::tuple<double, double, double>>& intervening,
                                      double min_opt_prop, double min_length, double window_length) {
    
    TestBonder bonder;
    bonder.min_opt_proportion = min_opt_prop;
    bonder.min_length = min_length;
    bonder.window_length = window_length;
    
    auto got = bonder.longest_windowed_partition(shared, intervening);
    auto expected = brute_force_longest_windowed_partition(shared, intervening, min_opt_prop, min_length, window_length);
    
    double score_got, score_expected;
    bool valid_got, valid_expected;
    
    tie(score_got, valid_got) = check_longest_windowed_partition(got, shared, intervening, min_opt_prop, min_length, window_length, true);
    tie(score_expected, valid_expected) = check_longest_windowed_partition(expected, shared, intervening, min_opt_prop, min_length, window_length);
    
    apply_test(shared, intervening, got, expected, score_got, score_expected, valid_got, valid_expected,
               min_opt_prop, min_length, window_length);
    
}

void do_tests(const std::vector<std::tuple<double, double, double>>& shared,
              const std::vector<std::tuple<double, double, double>>& intervening,
              double min_opt_prop, double min_length, double window_length) {
    
    test_longest_partition(shared, intervening, min_opt_prop, min_length);
    test_longest_windowed_partition(shared, intervening, min_opt_prop, min_length, window_length);
}


int main(int argc, char* argv[]) {
     
    random_device rd;
    default_random_engine gen(rd());
    
    {
        std::vector<std::tuple<double, double, double>> shared{
            {6.77295,0.221092,0.102639},
            {1.47505,0.60579,0.330855},
            {9.76515,0.0435343,0.0386732},
            {6.39488,0.226435,0.00340727},
            {6.94636,0.883661,0.463603}
        };
        std::vector<std::tuple<double, double, double>> intervening{
            {3.8005,0.843752,0.770744},
            {9.36865,0.243723,0.122616},
            {3.76464,0.685706,0.565402},
            {3.05252,0.882991,0.0544103}
        };
        double min_opt_prop = 0.562209;
        double min_length = 14.3296;
        double window_length = 12.3956;
        do_tests(shared, intervening, min_opt_prop, min_length, window_length);
    }
    
    {
        std::vector<std::tuple<double, double, double>> shared{
            {5.8195,0.902362,0.78278},
            {1.98314,0.730391,0.791058},
            {4.75961,0.31335,0.0279827},
            {4.34004,0.299881,0.205333},
            {6.9448,0.561213,0.328167}
        };
        std::vector<std::tuple<double, double, double>> intervening{
            {7.27921,0.217796,0.0542346},
            {8.2105,0.00892272,0.00840213},
            {2.94546,0.783561,0.282201},
            {3.09021,0.132805,0.0407658}
        };
        double min_opt_prop = 0.420019;
        double min_length = 1.02473;
        double window_length = 28.1783;
        do_tests(shared, intervening, min_opt_prop, min_length, window_length);
    }
    
    {
        std::vector<std::tuple<double, double, double>> shared{
            {8.03596,0.40715,0.0232886},
            {7.37587,0.926492,0.224252},
            {1.54726,0.441565,0.0993635},
            {5.87644,0.819689,0.173465},
            {6.12154,0.357603,0.295189}
        };
        std::vector<std::tuple<double, double, double>> intervening{
            {1.75105,0.996891,0.575165},
            {7.81787,0.195799,0.0551634},
            {2.10464,0.877787,0.798768},
            {1.50624,0.618965,0.482695}
        };
        double min_opt_prop = 0.354464;
        double min_length = 7.06864;
        double window_length = 12.3921;
        do_tests(shared, intervening, min_opt_prop, min_length, window_length);
    }
    
    {
        std::vector<std::tuple<double, double, double>> shared{
            {2.11585,0.8836,0.829148},
            {2.70877,0.803905,0.128772},
            {9.14373,0.286264,0.136492},
            {6.82174,0.099902,0.0403279},
            {1.66202,0.420703,0.245367}
        };
        std::vector<std::tuple<double, double, double>> intervening{
            {1.44801,0.532382,0.325309},
            {3.58604,0.0964342,0.00296816},
            {4.09387,0.14986,0.162214},
            {2.05962,0.513442,0.196855}
        };
        double min_opt_prop = 0.521271;
        double min_length = 5.19205;
        double window_length = 17.9201;
        do_tests(shared, intervening, min_opt_prop, min_length, window_length);
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
        double window_length = 1.0;
        do_tests(shared, intervening, min_opt_prop, min_length, window_length);
    }
    
    {
        std::vector<std::tuple<double, double, double>> shared{
            {2.22899,0.0859514,0.0454957},
            {7.53831,0.969004,0.778769},
            {4.53274,0.767475,0.491838},
            {8.62815,0.897002,0.845958},
            {2.01937,0.547733,0.57212}
        };
        std::vector<std::tuple<double, double, double>> intervening{
            {9.30858,0.61885,0.0527381},
            {7.26259,0.913271,0.791344},
            {1.98094,0.443888,0.105537},
            {3.38906,0.9802,0.736625}
        };
        double min_opt_prop = 0.129975;
        double min_length = 6.06805;
        double window_length = 2.73029;
        do_tests(shared, intervening, min_opt_prop, min_length, window_length);
    }
    
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
        double window_length = 1.0;
        do_tests(shared, intervening, min_opt_prop, min_length, window_length);
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
        double window_length = 1.0;
        do_tests(shared, intervening, min_opt_prop, min_length, window_length);
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
        double window_length = 1.0;
        do_tests(shared, intervening, min_opt_prop, min_length, window_length);
    }

    
    uniform_real_distribution<double> min_prop_distr(0.0, 1.0);
    uniform_real_distribution<double> length_distr(1.0, 10.0);
    uniform_real_distribution<double> min_length_distr(0.0, 20.0);
    uniform_real_distribution<double> window_length_distr(1.0, 30.0);
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
            double window_length = window_length_distr(gen);
            
            do_tests(shared, intervening, min_prop, min_length, window_length);
        }
    }
    
    
    cerr << "passed all tests!" << endl;
}
