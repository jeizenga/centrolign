#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <unordered_set>
#include <vector>
#include <list>
#include <string>
#include <random>
#include <algorithm>
#include <iostream>

#include "centrolign/threading.hpp"

using namespace std;
using namespace centrolign;

int main(int argc, char* argv[]) {
     
    random_device rd;
    default_random_engine gen(rd());
    
    {
        Threading::set_num_threads(4);
        
        size_t size = 1000;
        std::vector<int> results(size, -1);
        
        function<void(size_t)> lambda = [&](size_t i) {
            results[i] = i * i;
        };
        
        Threading::parallel_for(size, lambda);
        
        for (size_t i = 0; i < size; ++i) {
            assert(results[i] == i * i);
        }
    }
    
    
    
    cerr << "passed all tests!" << endl;
}
