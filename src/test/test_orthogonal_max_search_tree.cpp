#include <cstdio>
#include <cstdlib>
#include <cstdlib>

#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <cassert>
#include <limits>

#include "centrolign/orthogonal_max_search_tree.hpp"

using namespace std;
using namespace centrolign;

template<class Generator>
vector<tuple<double, double, double>> random_data(size_t size, Generator& gen) {
    
    uniform_real_distribution<double> distr(0.0, 1.0);
    
    vector<tuple<double, double, double>> data;
    for (size_t i = 0; i < size; ++i) {
        data.emplace_back(distr(gen), distr(gen), distr(gen));
    }
    return data;
}

template<class Generator>
vector<tuple<double, double, double>> random_updates(const vector<tuple<double, double, double>>& data,
                                                      size_t num_updates, Generator& gen) {
    
    uniform_real_distribution<double> val_distr(0.0, 1.0);
    uniform_int_distribution<size_t> data_distr(0, data.size() - 1);
    
    vector<tuple<double, double, double>> updates;
    for (size_t i = 0; i < num_updates; ++i) {
        size_t idx = data_distr(gen);
        updates.emplace_back(get<0>(data[idx]), get<1>(data[idx]), val_distr(gen));
    }
    
    return updates;
}

double mininf = std::numeric_limits<double>::lowest();

double range_max(vector<tuple<double, double, double>>& data,
                 double m1, double M1, double m2, double M2) {
    
    double rm = mininf;
    for (auto d : data) {
        double x, y, v;
        tie(x, y, v) = d;
        if (x >= m1 && x < M1 && y >= m2 && y < M2) {
            rm = max(rm, v);
        }
    }
    return rm;
}

void apply_update(OrthogonalMaxSearchTree<double, double, double>& tree,
                  vector<tuple<double, double, double>>& data,
                  tuple<double, double, double>& update) {
    
    bool found = false;
    
    auto it = tree.find(get<0>(update), get<1>(update));
    tree.update(it, get<2>(update));
    for (size_t i = 0; i < data.size(); ++i) {
        auto& d = data[i];
        if (get<0>(d) == get<0>(update) && get<1>(d) == get<1>(update)) {
            get<2>(d) = get<2>(update);
            found = true;
            break;
        }
    }
    
    assert(found);
}

void print_inputs(vector<tuple<double, double, double>>& data,
                  vector<tuple<double, double, double>>& updates) {
    
    cerr << "data:\n";
    for (auto d : data) {
        cerr << std::get<0>(d) << '\t' << std::get<1>(d) << '\t' << std::get<2>(d) << '\n';
    }
    cerr << "updates:\n";
    for (auto d : updates) {
        cerr << std::get<0>(d) << '\t' << std::get<1>(d) << '\t' << std::get<2>(d) << '\n';
    }
}

bool test_queries(OrthogonalMaxSearchTree<double, double, double>& tree,
                  vector<tuple<double, double, double>>& data) {
    
    
    auto its1 = minmax_element(data.begin(), data.end(),
                               [](tuple<double, double, double> a,
                                  tuple<double, double, double> b){
        return std::get<0>(a) < std::get<0>(b);
    });
    
    auto its2 = minmax_element(data.begin(), data.end(),
                               [](tuple<double, double, double> a,
                                  tuple<double, double, double> b){
        return std::get<1>(a) < std::get<1>(b);
    });
    
    double m1 = get<0>(*its1.first);
    double M1 = get<0>(*its1.second);
    double m2 = get<0>(*its1.first);
    double M2 = get<0>(*its2.second);
    
    double gm1 = m1 - 0.1 * (M1 - m1);
    double gm2 = m2 - 0.1 * (M1 - m1);
    double gM1 = M1 + 0.1 * (M1 - m1);
    double gM2 = M2 + 0.1 * (M2 - m2);
    
    std::vector<double> grid1, grid2;
    for (size_t i = 0; i < 11; ++i) {
        grid1.push_back(gm1 + i * (gM1 - gm1) / 10.0);
        grid2.push_back(gm2 + i * (gM2 - gm2) / 10.0);
    }
    
    // also get the values themselves for edge cases
    for (auto& d : data) {
        grid1.push_back(std::get<0>(d));
        grid2.push_back(std::get<1>(d));
    }
    
    
    for (size_t i = 0; i < grid1.size(); ++i) {
        for (size_t j = i + 1; j < grid1.size(); ++j) {
            for (size_t k = 0; k < grid2.size(); ++k) {
                for (size_t l = k + 1; l < grid2.size(); ++l) {
                    
                    double m1 = grid1[i];
                    double M1 = grid1[j];
                    double m2 = grid2[k];
                    double M2 = grid2[l];
                    
                    auto it = tree.range_max(m1, M1, m2, M2);
                    auto direct = range_max(data, m1, M1, m2, M2);
                    
                    if (it == tree.end()) {
                        if (direct != mininf) {
                            cerr << "failed on non-empty query\n";
                            return false;
                        }
                    }
                    else {
                        if (get<2>(*it) != direct) {
                            cerr << "failed on incorrect query\n";
                            return false;
                        }
                    }
                }
            }
        }
    }
    return true;
}

void do_test(vector<tuple<double, double, double>>& data,
             vector<tuple<double, double, double>>& updates) {
    
    auto original = data;
    auto copy = data;
    OrthogonalMaxSearchTree<double, double, double> tree(data);
    
    bool pass = test_queries(tree, copy);
    if (!pass) {
        cerr << "fail on initial queries\n";
        print_inputs(original, updates);
        exit(1);
    }
    
    for (auto update : updates) {
        apply_update(tree, copy, update);
        bool pass = test_queries(tree, copy);
        if (!pass) {
            cerr << "fail on updated queries\n";
            print_inputs(original, updates);
            exit(1);
        }
    }
}



int main(int argc, char* argv[]) {
    
    {
        vector<tuple<double, double, double>> data{
            {0.0188351,    0.771493,    0.256978},
            {0.567013,    0.668294,    0.682313},
            {0.902988,    0.233216,    0.11039},
            {0.257719,    0.77038,    0.313141},
            {0.672538,    0.00697451,    0.411589},
            {0.727974,    0.114095,    0.888307}
        };
        vector<tuple<double, double, double>> updates{
            {0.567013,    0.668294,    0.996294}
        };
        do_test(data, updates);
    }
    
    
    {
        vector<tuple<double, double, double>> data{
            {0, 2, 0},
            {0, 1, 1},
            {2, 2, 2}
        };
        vector<tuple<double, double, double>> updates{
            {0, 2, 2},
            {0, 1, 0},
            {2, 2, 1}
        };
        do_test(data, updates);
    }

    random_device rd;
    default_random_engine gen(rd());
    
    uniform_int_distribution<size_t> size_distr(1, 12);
    uniform_int_distribution<size_t> num_updates_distr(1, 5);
    size_t num_tests = 10;
    for (size_t i = 0; i < num_tests; ++i) {
        size_t size = size_distr(gen);
        size_t num_updates = num_updates_distr(gen);
        auto data = random_data(size, gen);
        auto updates = random_updates(data, num_updates, gen);
        do_test(data, updates);
    }
    
    
    cerr << "passed all tests!" << endl;
}
