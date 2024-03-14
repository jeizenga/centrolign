#include <cstdio>
#include <cstdlib>

#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <unordered_set>

#include "centrolign/range_unique_query.hpp"

using namespace std;
using namespace centrolign;

vector<unsigned int> random_vector(size_t size, default_random_engine& gen) {
    
    uniform_int_distribution<unsigned int> distr(0, size-1);
    
    vector<unsigned int> vec(size);
    for (size_t i = 0; i < vec.size(); ++i) {
        vec[i] = distr(gen);
    }
    return vec;
}

template<size_t N>
void do_test_internal(vector<unsigned int>& vec) {
    
    RUQ<N> ruq(vec);
    
    for (size_t i = 0; i <= vec.size(); ++i) {
        for (size_t j = i; j <= vec.size(); ++j) {
            
            unordered_set<unsigned int> s(vec.begin() + i, vec.begin() + j);
            
            auto got = ruq.range_unique(i, j);
            if (got != s.size()) {
                
                cerr << "failed test for interval " << i << " " << j << " on data:\n";
                for (auto v : vec) {
                    cerr << '\t' << v << '\n';
                }
                cerr << "expected: " << s.size() << ", got: " << got << '\n';
                exit(1);
            }
        }
    }
}

void do_test(vector<unsigned int>& vec) {
    do_test_internal<2>(vec);
    do_test_internal<3>(vec);
    do_test_internal<4>(vec);
    do_test_internal<5>(vec);
}

int main(int argc, char* argv[]) {
    
    
    vector<vector<unsigned int>> vecs{
        {1, 1, 2, 2, 1, 2, 2, 1},
        {2, 2, 22, 10, 10, 11, 16, 10, 12, 18, 24, 7, 21, 21, 5, 5, 29, 14, 30, 4, 21, 1, 26, 27, 15, 12, 18, 11, 27, 11, 27, 15},
        {7, 27, 0, 16, 12, 28, 2, 8, 25, 7, 3, 0, 28, 9, 23, 29, 16, 15, 23, 19, 19, 24, 7, 13, 23, 21, 6, 15, 27, 9, 24, 11},
        {0, 1, 2, 3, 4, 6, 5, 7, 8, 9}
    };
    for (auto& vec : vecs) {
        do_test(vec);
    }
    
    random_device rd;
    default_random_engine gen(rd());
    
    uniform_int_distribution<size_t> size_distr(1, 60);
    size_t num_reps = 20;
    for (size_t rep = 0; rep < num_reps; ++rep) {
        size_t size = size_distr(gen);
        auto rvec = random_vector(size, gen);
        do_test(rvec);
    }
    
    cerr << "passed all tests!" << endl;
}
