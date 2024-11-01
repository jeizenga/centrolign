#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <vector>
#include <string>
#include <random>
#include <sstream>
#include <algorithm>
#include <unordered_set>

#include "centrolign/bit_packed_vector.hpp"
#include "centrolign/test_util.hpp"

using namespace std;
using namespace centrolign;

void test_packed_vector(const std::vector<uint64_t>& input) {
    
    PackedVector vec(input.size());
    for (size_t i = 0; i < input.size(); ++i) {
        vec.set(i, input[i]);
    }
    for (size_t i = 0; i < input.size(); ++i) {
        uint64_t got = vec.at(i);
        uint64_t expected = input.at(i);
        if (got != expected) {
            std::cerr << "got " << got << " at position " << i << " with input vector:\n";
            std::cerr << '{';
            for (size_t i = 0; i < input.size(); ++i) {
                if (i) {
                    std::cerr << ", ";
                }
                std::cerr << input[i];
            }
            std::cerr << "}\n";
            exit(1);
        }
    }
}

int main(int argc, char* argv[]) {
     
    random_device rd;
    default_random_engine gen(rd());
    
    {
        PackedVector vec;
        assert(vec.size() == 0);
        assert(vec.empty());
    }

    {
        std::unordered_set<size_t> ones{2, 4, 78, 92, 56, 6, 14, 19, 55, 64, 63};

        PackedVector vec(256);
        assert(vec.size() == 256);
        assert(!vec.empty());

        for (auto i : ones) {
            vec.set(i, 1);
        }
        for (size_t i = 0; i < vec.size(); ++i) {
            if (ones.count(i)) {
                assert(vec.at(i) == 1);
            }
            else {
                assert(vec.at(i) == 0);
            }
        }
    }

    test_packed_vector({0, 1, 2, 3, 4, 5, 6, 7, 8});
    
    test_packed_vector({111, 250, 250, 98, 157, 203, 31, 110, 5, 123, 122, 28, 169, 85, 74, 59, 227, 226, 198, 236});

    uniform_int_distribution<int> size_distr(1, 50);
    uniform_int_distribution<uint64_t> val_distr(0, 512);
    size_t num_reps = 1000;
    for (size_t rep = 0; rep < num_reps; ++rep) {
        std::vector<uint64_t> input(size_distr(gen));
        for (size_t i = 0; i < input.size(); ++i) {
            input[i] = val_distr(gen);
        }
        test_packed_vector(input);
    }
    
    
    cerr << "passed all tests!" << endl;
}
