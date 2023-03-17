#include <cstdio>
#include <cstdlib>
#include <cstdlib>

#include <iostream>
#include <vector>
#include <algorithm>

#include "centrolign/range_min_query.hpp"

using namespace std;
using namespace centrolign;

template<class T>
class TestRMQ : public RMQ<T> {
public:
    using typename RMQ<T>::ExhaustiveRMQ;
    using typename RMQ<T>::SparseTable;
};

template<class RangeMinQuery>
bool check_rmq(const RangeMinQuery& rmq, const vector<int>& vec) {
    
    for (size_t i = 0; i < vec.size(); ++i) {
        for (size_t j = i + 1; j <= vec.size(); ++j) {
            
            size_t arg_min = i;
            for (size_t k = i + 1; k < j; ++k) {
                if (vec[k] < vec[arg_min]) {
                    arg_min = k;
                }
            }
            
            if (rmq.range_arg_min(i, j) != arg_min) {
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
    if (!check_rmq(sparse_table, vec)) {
        failed_on = "sparse table";
    }
    if (!check_rmq(exhaustive_rmq, vec)) {
        failed_on = "exhaustive table";
    }
    if (!check_rmq(rmq, vec)) {
        failed_on = "RMQ";
    }
    if (!failed_on.empty()) {
        cerr << "failed tests on " << failed_on << " with data:" << endl;
        for (int i = 0; i < vec.size(); ++i) {
            if (i) {
                cerr << ", ";
            }
            cerr << i;
        }
        cerr << endl;
        exit(1);
    }
}

int main(int argc, char* argv[]) {
    
    vector<int> vec{0, 1, 2, 3, 4, 6, 5, 7, 8, 9};
    do_test(vec);
    
    cerr << "passed all tests!" << endl;
}
