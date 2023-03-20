#include <cstdio>
#include <cstdlib>
#include <cstdlib>

#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <utility>

#include "centrolign/utility.hpp"
#include "centrolign/integer_sort.hpp"

using namespace std;
using namespace centrolign;

int main(int argc, char* argv[]) {
    
    random_device rd;
    default_random_engine gen(rd());
    
    // integer logarithm test
    {
        vector<pair<int, int>> hi_bit_problems{
            {1, 0},
            {2, 1},
            {3, 1},
            {4, 2},
            {5, 2},
            {6, 2},
            {7, 2},
            {8, 3},
            {9, 3},
            {10, 3},
            {100, 6},
            {1000, 9},
            {10000, 13},
            {100000, 16},
            {1000000, 19},
            {10000000, 23},
            {100000000, 26},
            {1000000000, 29}
        };
        for (auto& problem : hi_bit_problems) {
            if (hi_bit(problem.first) != problem.second) {
                cerr << "hi bit test failed on problem " << problem.first << " -> " << problem.second << " with result " << hi_bit(problem.first) << '\n';
                exit(1);
            }
        }
    }
    
    // integer sort tests
    {
        size_t num_reps = 5;
        for (int size : {0, 1, 2, 4, 8, 16}) {
            for (size_t i = 0; i < num_reps; ++i) {
                std::uniform_int_distribution<int> distr(0, size);
                std::vector<std::pair<int, int>> to_sort;
                for (int j = 0; j < size; ++j) {
                    to_sort.emplace_back(distr(gen), distr(gen));
                }
                auto indexes = range_vector(to_sort.size());
                indexes = integer_sort(indexes, [&](size_t idx) { return to_sort[idx].second; });
                indexes = integer_sort(indexes, [&](size_t idx) { return to_sort[idx].first; });
                for (size_t j = 1; j < indexes.size(); ++j) {
                    if (to_sort[indexes[j - 1]] > to_sort[indexes[j]]) {
                        cerr << "integer sort failed on input" << endl;
                        for (auto p : to_sort) {
                            cerr << p.first << ", " << p.second << '\n';
                        }
                        exit(1);
                    }
                }
            }
        }
    }
    
    // test reverse adapter
    {
        auto range = range_vector(100);
        auto copy = range;
        std::reverse(copy.begin(), copy.end());
        auto it = copy.begin();
        for (auto x : ReverseForEachAdapter<std::vector<size_t>>(range)) {
            if (x != *it) {
                cerr << "reverse for each failed at " << x << '\n';
                exit(1);
            }
            ++it;
        }
    }
    
    cerr << "passed all tests!" << endl;
    return 0;
}
