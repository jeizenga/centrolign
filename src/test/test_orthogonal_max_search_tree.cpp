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
    
    // real data that led to a failure
    {
        using key_t = tuple<int64_t, size_t, size_t, size_t>;
        vector<tuple<key_t, size_t, double>> data{
            {{5,0,0,1}, 16, -1.79769e+308},
            {{-4,1,0,0}, 14, -1.79769e+308},
            {{2,2,0,0}, 1, -1.79769e+308},
            {{1,2,0,1}, 2, -1.79769e+308},
            {{-2,2,0,2}, 5, -1.79769e+308},
            {{-3,2,0,3}, 6, -1.79769e+308},
            {{-4,2,0,4}, 7, -1.79769e+308},
            {{-7,2,0,5}, 10, -1.79769e+308},
            {{-8,2,0,6}, 11, -1.79769e+308},
            {{-9,2,0,8}, 12, -1.79769e+308},
            {{-9,2,0,9}, 12, -1.79769e+308},
            {{-10,2,0,11}, 13, -1.79769e+308},
            {{-14,2,0,12}, 17, -1.79769e+308},
            {{-15,2,0,13}, 18, -1.79769e+308},
            {{-16,2,0,14}, 19, -1.79769e+308},
            {{-10,2,0,15}, 13, -1.79769e+308},
            {{-15,2,0,16}, 18, -1.79769e+308},
            {{-14,2,0,17}, 17, -1.79769e+308},
            {{3,2,1,0}, 1, -1.79769e+308},
            {{2,2,1,1}, 2, -1.79769e+308},
            {{-1,2,1,2}, 5, -1.79769e+308},
            {{-2,2,1,3}, 6, -1.79769e+308},
            {{-3,2,1,4}, 7, -1.79769e+308},
            {{-6,2,1,5}, 10, -1.79769e+308},
            {{-7,2,1,6}, 11, -1.79769e+308},
            {{-8,2,1,8}, 12, -1.79769e+308},
            {{-8,2,1,9}, 12, -1.79769e+308},
            {{-9,2,1,11}, 13, -1.79769e+308},
            {{-13,2,1,12}, 17, -1.79769e+308},
            {{-14,2,1,13}, 18, -1.79769e+308},
            {{-15,2,1,14}, 19, -1.79769e+308},
            {{-9,2,1,15}, 13, -1.79769e+308},
            {{-14,2,1,16}, 18, -1.79769e+308},
            {{-13,2,1,17}, 17, -1.79769e+308},
            {{7,2,2,0}, 1, -1.79769e+308},
            {{6,2,2,1}, 2, -1.79769e+308},
            {{3,2,2,2}, 5, -1.79769e+308},
            {{2,2,2,3}, 6, -1.79769e+308},
            {{1,2,2,4}, 7, -1.79769e+308},
            {{-2,2,2,5}, 10, -1.79769e+308},
            {{-3,2,2,6}, 11, -1.79769e+308},
            {{-4,2,2,8}, 12, -1.79769e+308},
            {{-4,2,2,9}, 12, -1.79769e+308},
            {{-5,2,2,11}, 13, -1.79769e+308},
            {{-9,2,2,12}, 17, -1.79769e+308},
            {{-10,2,2,13}, 18, -1.79769e+308},
            {{-11,2,2,14}, 19, -1.79769e+308},
            {{-5,2,2,15}, 13, -1.79769e+308},
            {{-10,2,2,16}, 18, -1.79769e+308},
            {{-9,2,2,17}, 17, -1.79769e+308},
            {{11,2,3,0}, 1, -1.79769e+308},
            {{10,2,3,1}, 2, -1.79769e+308},
            {{7,2,3,2}, 5, -1.79769e+308},
            {{6,2,3,3}, 6, -1.79769e+308},
            {{5,2,3,4}, 7, -1.79769e+308},
            {{2,2,3,5}, 10, -1.79769e+308},
            {{1,2,3,6}, 11, -1.79769e+308},
            {{0,2,3,8}, 12, -1.79769e+308},
            {{0,2,3,9}, 12, -1.79769e+308},
            {{-1,2,3,11}, 13, -1.79769e+308},
            {{-5,2,3,12}, 17, -1.79769e+308},
            {{-6,2,3,13}, 18, -1.79769e+308},
            {{-7,2,3,14}, 19, -1.79769e+308},
            {{-1,2,3,15}, 13, -1.79769e+308},
            {{-6,2,3,16}, 18, -1.79769e+308},
            {{-5,2,3,17}, 17, -1.79769e+308},
            {{16,2,4,0}, 1, -1.79769e+308},
            {{15,2,4,1}, 2, -1.79769e+308},
            {{12,2,4,2}, 5, -1.79769e+308},
            {{11,2,4,3}, 6, -1.79769e+308},
            {{10,2,4,4}, 7, -1.79769e+308},
            {{7,2,4,5}, 10, -1.79769e+308},
            {{6,2,4,6}, 11, -1.79769e+308},
            {{5,2,4,8}, 12, -1.79769e+308},
            {{5,2,4,9}, 12, -1.79769e+308},
            {{4,2,4,11}, 13, -1.79769e+308},
            {{0,2,4,12}, 17, -1.79769e+308},
            {{-1,2,4,13}, 18, -1.79769e+308},
            {{-2,2,4,14}, 19, -1.79769e+308},
            {{4,2,4,15}, 13, -1.79769e+308},
            {{-1,2,4,16}, 18, -1.79769e+308},
            {{0,2,4,17}, 17, -1.79769e+308},
            {{1,2,6,0}, 1, -1.79769e+308},
            {{0,2,6,1}, 2, -1.79769e+308},
            {{-3,2,6,2}, 5, -1.79769e+308},
            {{-4,2,6,3}, 6, -1.79769e+308},
            {{-5,2,6,4}, 7, -1.79769e+308},
            {{-8,2,6,5}, 10, -1.79769e+308},
            {{-9,2,6,6}, 11, -1.79769e+308},
            {{-10,2,6,8}, 12, -1.79769e+308},
            {{-10,2,6,9}, 12, -1.79769e+308},
            {{-11,2,6,11}, 13, -1.79769e+308},
            {{-15,2,6,12}, 17, -1.79769e+308},
            {{-16,2,6,13}, 18, -1.79769e+308},
            {{-17,2,6,14}, 19, -1.79769e+308},
            {{-11,2,6,15}, 13, -1.79769e+308},
            {{-16,2,6,16}, 18, -1.79769e+308},
            {{-15,2,6,17}, 17, -1.79769e+308},
            {{3,2,7,0}, 1, -1.79769e+308},
            {{2,2,7,1}, 2, -1.79769e+308},
            {{-1,2,7,2}, 5, -1.79769e+308},
            {{-2,2,7,3}, 6, -1.79769e+308},
            {{-3,2,7,4}, 7, -1.79769e+308},
            {{-6,2,7,5}, 10, -1.79769e+308},
            {{-7,2,7,6}, 11, -1.79769e+308},
            {{-8,2,7,8}, 12, -1.79769e+308},
            {{-8,2,7,9}, 12, -1.79769e+308},
            {{-9,2,7,11}, 13, -1.79769e+308},
            {{-13,2,7,12}, 17, -1.79769e+308},
            {{-14,2,7,13}, 18, -1.79769e+308},
            {{-15,2,7,14}, 19, -1.79769e+308},
            {{-9,2,7,15}, 13, -1.79769e+308},
            {{-14,2,7,16}, 18, -1.79769e+308},
            {{-13,2,7,17}, 17, -1.79769e+308},
            {{2,3,0,1}, 3, -1.79769e+308},
            {{-3,3,0,2}, 8, -1.79769e+308},
            {{-15,3,0,6}, 20, -1.79769e+308},
            {{6,3,2,1}, 3, -1.79769e+308},
            {{1,3,2,2}, 8, -1.79769e+308},
            {{-11,3,2,6}, 20, -1.79769e+308},
            {{10,3,4,1}, 3, -1.79769e+308},
            {{5,3,4,2}, 8, -1.79769e+308},
            {{-7,3,4,6}, 20, -1.79769e+308},
            {{10,3,5,1}, 3, -1.79769e+308},
            {{5,3,5,2}, 8, -1.79769e+308},
            {{-7,3,5,6}, 20, -1.79769e+308},
            {{15,3,6,1}, 3, -1.79769e+308},
            {{10,3,6,2}, 8, -1.79769e+308},
            {{-2,3,6,6}, 20, -1.79769e+308},
            {{19,3,7,1}, 3, -1.79769e+308},
            {{14,3,7,2}, 8, -1.79769e+308},
            {{2,3,7,6}, 20, -1.79769e+308},
            {{19,3,8,1}, 3, -1.79769e+308},
            {{14,3,8,2}, 8, -1.79769e+308},
            {{2,3,8,6}, 20, -1.79769e+308},
            {{-2,4,1,0}, 4, -1.79769e+308},
            {{-7,4,1,1}, 9, -1.79769e+308},
            {{-14,4,1,4}, 16, -1.79769e+308},
            {{-19,4,1,5}, 21, -1.79769e+308},
            {{0,4,1,6}, 2, -1.79769e+308},
            {{-7,4,1,7}, 9, -1.79769e+308},
            {{3,4,2,0}, 4, -1.79769e+308},
            {{-2,4,2,1}, 9, -1.79769e+308},
            {{-9,4,2,4}, 16, -1.79769e+308},
            {{-14,4,2,5}, 21, -1.79769e+308},
            {{5,4,2,6}, 2, -1.79769e+308},
            {{-2,4,2,7}, 9, -1.79769e+308},
            {{12,4,3,0}, 4, -1.79769e+308},
            {{7,4,3,1}, 9, -1.79769e+308},
            {{0,4,3,4}, 16, -1.79769e+308},
            {{-5,4,3,5}, 21, -1.79769e+308},
            {{14,4,3,6}, 2, -1.79769e+308},
            {{7,4,3,7}, 9, -1.79769e+308},
            {{12,4,4,0}, 4, -1.79769e+308},
            {{7,4,4,1}, 9, -1.79769e+308},
            {{0,4,4,4}, 16, -1.79769e+308},
            {{-5,4,4,5}, 21, -1.79769e+308},
            {{14,4,4,6}, 2, -1.79769e+308},
            {{7,4,4,7}, 9, -1.79769e+308},
            {{17,4,5,0}, 4, -1.79769e+308},
            {{12,4,5,1}, 9, -1.79769e+308},
            {{5,4,5,4}, 16, -1.79769e+308},
            {{0,4,5,5}, 21, -1.79769e+308},
            {{19,4,5,6}, 2, -1.79769e+308},
            {{12,4,5,7}, 9, -1.79769e+308},
            {{7,4,6,0}, 4, -1.79769e+308},
            {{2,4,6,1}, 9, -1.79769e+308},
            {{-5,4,6,4}, 16, -1.79769e+308},
            {{-10,4,6,5}, 21, -1.79769e+308},
            {{9,4,6,6}, 2, -1.79769e+308},
            {{2,4,6,7}, 9, -1.79769e+308}
        };
        vector<tuple<key_t, size_t, double>> updates{
            {{1,2,6,0}, 1, 0.00694444},
            {{0,2,6,1}, 2, 0.00694444},
            {{-3,2,6,2}, 5, 0.00694444},
            {{-4,2,6,3}, 6, 0.00694444},
            {{-5,2,6,4}, 7, 0.00694444},
            {{-8,2,6,5}, 10, 0.00694444},
            {{-9,2,6,6}, 11, 0.00694444},
            {{-10,2,6,8}, 12, 0.00694444},
            {{-10,2,6,9}, 12, 0.00694444},
            {{-11,2,6,11}, 13, 0.00694444},
            {{-15,2,6,12}, 17, 0.00694444},
            {{-16,2,6,13}, 18, 0.00694444},
            {{-17,2,6,14}, 19, 0.00694444},
            {{-11,2,6,15}, 13, 0.00694444},
            {{-16,2,6,16}, 18, 0.00694444},
            {{-15,2,6,17}, 17, 0.00694444},
            {{-2,4,1,0}, 4, 0.0178571},
            {{-7,4,1,1}, 9, 0.0178571},
            {{-14,4,1,4}, 16, 0.0178571},
            {{-19,4,1,5}, 21, 0.0178571},
            {{0,4,1,6}, 2, 0.0178571},
            {{-7,4,1,7}, 9, 0.0178571},
            {{2,2,0,0}, 1, 0.00694444},
            {{1,2,0,1}, 2, 0.00694444},
            {{-2,2,0,2}, 5, 0.00694444},
            {{-3,2,0,3}, 6, 0.0248016},
            {{-4,2,0,4}, 7, 0.00694444},
            {{-7,2,0,5}, 10, 0.00694444},
            {{-8,2,0,6}, 11, 0.0248016},
            {{-9,2,0,8}, 12, 0.00694444},
            {{-9,2,0,9}, 12, 0.00694444},
            {{-10,2,0,11}, 13, 0.00694444},
            {{-14,2,0,12}, 17, 0.00694444},
            {{-15,2,0,13}, 18, 0.0248016},
            {{-16,2,0,14}, 19, 0.00694444},
            {{-10,2,0,15}, 13, 0.0248016},
            {{-15,2,0,16}, 18, 0.0248016},
            {{-14,2,0,17}, 17, 0.00694444},
            {{3,2,1,0}, 1, 0.00694444},
            {{2,2,1,1}, 2, 0.00694444},
            {{-1,2,1,2}, 5, 0.0248016},
            {{-2,2,1,3}, 6, 0.0248016},
            {{-3,2,1,4}, 7, 0.0248016},
            {{-6,2,1,5}, 10, 0.00694444},
            {{-7,2,1,6}, 11, 0.0248016},
            {{-8,2,1,8}, 12, 0.0248016},
            {{-8,2,1,9}, 12, 0.00694444},
            {{-9,2,1,11}, 13, 0.0138889},
            {{-13,2,1,12}, 17, 0.00694444},
            {{-14,2,1,13}, 18, 0.0248016},
            {{-15,2,1,14}, 19, 0.0248016},
            {{-9,2,1,15}, 13, 0.0248016},
            {{-14,2,1,16}, 18, 0.0248016},
            {{-13,2,1,17}, 17, 0.0138889},
            {{3,2,7,0}, 1, 0.00694444},
            {{2,2,7,1}, 2, 0.00694444},
            {{-1,2,7,2}, 5, 0.00694444},
            {{-2,2,7,3}, 6, 0.00694444},
            {{-3,2,7,4}, 7, 0.00694444},
            {{-6,2,7,5}, 10, 0.00694444},
            {{-7,2,7,6}, 11, 0.00694444},
            {{-8,2,7,8}, 12, 0.00694444},
            {{-8,2,7,9}, 12, 0.00694444},
            {{-9,2,7,11}, 13, 0.00694444},
            {{-13,2,7,12}, 17, 0.00694444},
            {{-14,2,7,13}, 18, 0.00694444},
            {{-15,2,7,14}, 19, 0.00694444},
            {{-9,2,7,15}, 13, 0.00694444},
            {{-14,2,7,16}, 18, 0.00694444},
            {{-13,2,7,17}, 17, 0.00694444},
            {{2,3,0,1}, 3, 0.0228175},
            {{-3,3,0,2}, 8, 0.0406746},
            {{-15,3,0,6}, 20, 0.0406746},
            {{3,4,2,0}, 4, 0.0248016},
            {{-2,4,2,1}, 9, 0.0426587},
            {{-9,4,2,4}, 16, 0.0426587},
            {{-14,4,2,5}, 21, 0.0426587},
            {{5,4,2,6}, 2, 0.0178571},
            {{-2,4,2,7}, 9, 0.0426587},
            {{7,2,2,0}, 1, 0.00694444},
            {{6,2,2,1}, 2, 0.00694444},
            {{3,2,2,2}, 5, 0.0138889},
            {{2,2,2,3}, 6, 0.0297619},
            {{1,2,2,4}, 7, 0.0138889},
            {{-2,2,2,5}, 10, 0.031746},
            {{-3,2,2,6}, 11, 0.047619},
            {{-4,2,2,8}, 12, 0.0138889},
            {{-4,2,2,9}, 12, 0.031746},
            {{-5,2,2,11}, 13, 0.031746},
            {{-9,2,2,12}, 17, 0.031746},
            {{-10,2,2,13}, 18, 0.0406746},
            {{-11,2,2,14}, 19, 0.031746},
            {{-5,2,2,15}, 13, 0.047619},
            {{-10,2,2,16}, 18, 0.047619},
            {{-9,2,2,17}, 17, 0.031746},
            {{6,3,2,1}, 3, 0.015873},
            {{1,3,2,2}, 8, 0.0228175},
            {{-11,3,2,6}, 20, 0.0406746},
            {{-4,1,0,0}, 14, 1.00694},
            {{7,4,6,0}, 4, 0.0337302},
            {{2,4,6,1}, 9, 0.047619},
            {{-5,4,6,4}, 16, 0.0654762},
            {{-10,4,6,5}, 21, 0.0654762},
            {{9,4,6,6}, 2, 0.0178571},
            {{2,4,6,7}, 9, 0.047619},
            {{11,2,3,0}, 1, 0.00694444},
            {{10,2,3,1}, 2, 0.00694444},
            {{7,2,3,2}, 5, 0.0228175},
            {{6,2,3,3}, 6, 0.0228175},
            {{5,2,3,4}, 7, 0.0248016},
            {{2,2,3,5}, 10, 0.0367063},
            {{1,2,3,6}, 11, 0.0367063},
            {{0,2,3,8}, 12, 0.0297619},
            {{0,2,3,9}, 12, 0.0545635},
            {{-1,2,3,11}, 13, 0.0545635},
            {{-5,2,3,12}, 17, 1.01389},
            {{-6,2,3,13}, 18, 0.0634921},
            {{-7,2,3,14}, 19, 0.0634921},
            {{-1,2,3,15}, 13, 0.0634921},
            {{-6,2,3,16}, 18, 0.0634921},
            {{-5,2,3,17}, 17, 0.0545635},
            {{10,3,4,1}, 3, 0.015873},
            {{5,3,4,2}, 8, 0.0228175},
            {{-7,3,4,6}, 20, 0.0406746},
            {{10,3,5,1}, 3, 0.015873},
            {{5,3,5,2}, 8, 0.0337302},
            {{-7,3,5,6}, 20, 0.0724206},
        };
        
        OrthogonalMaxSearchTree<key_t, size_t, double> tree(data);
        
        std::vector<int64_t> q_values;
        for (auto& v : data) {
            q_values.push_back(std::get<0>(std::get<0>(v)));
        }
        q_values.resize(unique(q_values.begin(), q_values.end()) - q_values.begin());
        
        auto do_queries = [&]() {
            size_t m = 0, M = 20;
            for (auto q : q_values) {
                auto it = tree.range_max(key_t(q, 0, 0, 0), key_t(q+1, 0, 0, 0), m, M);
                double rm = mininf;
                size_t max_idx = -1;
                for (size_t i = 0; i < data.size(); ++i) {
                    auto d = data[i];
                    if (std::get<0>(std::get<0>(d)) == q && std::get<1>(d) >= m && std::get<1>(d) < M) {
                        if (get<2>(d) > rm) {
                            rm = get<2>(d);
                            max_idx = i;
                        }
                        rm = max(rm, std::get<2>(d));
                    }
                }
                if (rm == mininf) {
                    assert(it == tree.end());
                }
                else {
                    if (it == tree.end() || rm != std::get<2>(*it)) {
                        std::cerr << "fail on query " << q << ", got ";
                        if (it == tree.end()) {
                            std::cerr << ".";
                        }
                        else {
                            cerr << "(" << get<0>(get<0>(*it)) << "," << get<1>(get<0>(*it)) << "," << get<2>(get<0>(*it)) << "," << get<3>(get<0>(*it)) << "), " << get<1>(*it) << ": " << get<2>(*it);
                        }
                        cerr << ", expected " << rm << '\n';
                        std::cerr << "contents:\n";
                        for (auto v : tree) {
                            cerr << "(" << get<0>(get<0>(v)) << "," << get<1>(get<0>(v)) << "," << get<2>(get<0>(v)) << "," << get<3>(get<0>(v)) << "), " << get<1>(v) << ": " << get<2>(v) << '\n';
                        }
                        cerr << "data:\n";
                        for (size_t i = 0; i < data.size(); ++i) {
                            if (i  == max_idx) {
                                std::cerr << "***";
                            }
                            auto v = data[i];
                            cerr << i << ": (" << get<0>(get<0>(v)) << "," << get<1>(get<0>(v)) << "," << get<2>(get<0>(v)) << "," << get<3>(get<0>(v)) << "), " << get<1>(v) << ": " << get<2>(v) << '\n';
                        }
                        exit(1);
                    }
                    
                }
            }
        };
        
        auto do_update = [&](key_t k, size_t i, double v) {
            for (auto& d : data) {
                if (std::get<0>(d) == k && std::get<1>(d) == i) {
                    std::get<2>(d) = v;
                }
            }
            tree.update(tree.find(k, i), v);
        };
        for (auto update : updates) {
            do_update(std::get<0>(update), std::get<1>(update), std::get<2>(update));
            do_queries();
        }
    }
    
    
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
