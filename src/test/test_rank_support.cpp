#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <vector>
#include <string>
#include <random>
#include <sstream>
#include <algorithm>
#include <unordered_set>

#include "centrolign/rank_support.hpp"
#include "centrolign/packed_vector.hpp"
#include "centrolign/test_util.hpp"

using namespace std;
using namespace centrolign;

void test_rank(const std::vector<int>& input) {
    
    PackedVector bit_vec(input.size());
    for (size_t i = 0; i < input.size(); ++i) {
        if (input[i]) {
            bit_vec[i] = 1;
        }
    }
    
    RankSupport rank_support(bit_vec);
    
    size_t rank = 0;
    for (size_t i = 0; i < input.size(); ++i) {
        size_t got = rank_support.rank(i);
        assert((int) rank_support.at(i) == bit_vec.at(i));
        if (got != rank) {
            std::cerr << "did not get expected result for rank of index " << i << ". expected " << rank << ", got " << got << '\n';
            std::cerr << "input:\n";
            for (size_t j = 0; j < input.size(); ++j) {
                std::cerr << j << '\t' << input[j] << '\n';
            }
            exit(1);
        }
        if (input[i]) {
            ++rank;
        }
    }
}


int main(int argc, char* argv[]) {
     
    random_device rd;
    default_random_engine gen(rd());
    
    {
        vector<int> bv{0};
        test_rank(bv);
    }
    
    {
        vector<int> bv{1};
        test_rank(bv);
    }
    
    {
        vector<int> bv;
        test_rank(bv);
    }

    {
        vector<int> bv{0, 0, 1, 0};
        test_rank(bv);
    }

    {
        vector<int> bv{0, 1, 1, 0, 1, 0, 0, 0, 1};
        test_rank(bv);
    }


    uniform_int_distribution<int> size_distr(0, 500);
    uniform_int_distribution<int> val_distr(0, 1);
    size_t num_reps = 5000;
    for (size_t rep = 0; rep < num_reps; ++rep) {
        std::vector<int> bv(size_distr(gen));
        for (size_t i = 0; i < bv.size(); ++i) {
            bv[i] = val_distr(gen);
        }
        test_rank(bv);
    }
    
    
    cerr << "passed all tests!" << endl;
}
