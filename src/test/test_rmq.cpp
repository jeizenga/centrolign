#include <cstdio>
#include <cstdlib>

#include <iostream>
#include <vector>
#include <algorithm>
#include <random>

#include "centrolign/range_min_query.hpp"

using namespace std;
using namespace centrolign;

template<class T>
class TestRMQ : public RMQ<T> {
public:
    using typename RMQ<T>::ExhaustiveRMQ;
    using typename RMQ<T>::SparseTable;
};

vector<int> random_vector(size_t size, default_random_engine& gen) {
    
    uniform_int_distribution<int> distr(0, size-1);
    
    vector<int> vec(size);
    for (size_t i = 0; i < vec.size(); ++i) {
        vec[i] = distr(gen);
    }
    return vec;
}

template<class RangeMinQuery>
bool check_rmq(const RangeMinQuery& rmq, const vector<int>& vec,
               size_t& begin, size_t& end, size_t& result, size_t& arg_min) {
    
    for (size_t i = 0; i < vec.size(); ++i) {
        for (size_t j = i + 1; j <= vec.size(); ++j) {
            
            arg_min = i;
            for (size_t k = i + 1; k < j; ++k) {
                if (vec[k] < vec[arg_min]) {
                    arg_min = k;
                }
            }
            result = rmq.range_arg_min(i, j);
            if (result != arg_min) {
                begin = i;
                end = j;
                return false;
            }
        }
    }
    return true;
}

void do_test(vector<int>& vec) {
    TestRMQ<int>::ExhaustiveRMQ exhaustive_rmq;
    exhaustive_rmq.initialize(vec.begin(), vec.end());
    TestRMQ<int>::SparseTable sparse_table(vec);
    RMQ<int> rmq(vec);
    
    string failed_on = "";
    size_t begin = 0, end = 0, result = 0, expected = 0;
    if (!check_rmq(rmq, vec, begin, end, result, expected)) {
        failed_on = "RMQ";
    }
    else if (!check_rmq(sparse_table, vec, begin, end, result, expected)) {
        failed_on = "sparse table";
    }
    else  if (!check_rmq(exhaustive_rmq, vec, begin, end, result, expected)) {
        failed_on = "exhaustive table";
    }
    if (!failed_on.empty()) {
        cerr << "failed tests on " << failed_on << ", query " << begin << ":" << end << ", result " << result << ", expected " << expected << ", data" << endl;
        for (int i = 0; i < vec.size(); ++i) {
            if (i) {
                cerr << ", ";
            }
            cerr << vec[i];
        }
        cerr << endl;
        exit(1);
    }
}

int main(int argc, char* argv[]) {
    
    
    vector<vector<int>> vecs{
        {2, 2, 22, 10, 10, 11, 16, 10, 12, 18, 24, 7, 21, 21, 5, 5, 29, 14, 30, 4, 21, 1, 26, 27, 15, 12, 18, 11, 27, 11, 27, 15},
        {7, 27, 0, 16, 12, 28, 2, 8, 25, 7, 3, 0, 28, 9, 23, 29, 16, 15, 23, 19, 19, 24, 7, 13, 23, 21, 6, 15, 27, 9, 24, 11},
        {0, 1, 2, 3, 4, 6, 5, 7, 8, 9}
    };
    for (auto& vec : vecs) {
        do_test(vec);
    }
    
    random_device rd;
    default_random_engine gen(rd());
    
    size_t num_reps = 3;
    for (size_t size : {4, 16, 32, 64, 256, 512}) {
        for (size_t rep = 0; rep < num_reps; ++rep) {
            
            auto rvec = random_vector(size, gen);
            do_test(rvec);
        }
    }
    
    cerr << "passed all tests!" << endl;
}
