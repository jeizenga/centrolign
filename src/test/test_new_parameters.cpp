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

#include "centrolign/new_parameters.hpp"
#include "centrolign/test_util.hpp"

using namespace std;
using namespace centrolign;

int main(int argc, char* argv[]) {
     
    random_device rd;
    default_random_engine gen(rd());
    
    NewParameters params;
    cerr << params.generate_config() << '\n';
    
    cerr << "passed all tests!" << endl;
}
