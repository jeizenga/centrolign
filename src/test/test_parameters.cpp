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

#include "centrolign/parameters.hpp"
#include "centrolign/test_util.hpp"

using namespace std;
using namespace centrolign;

int main(int argc, char* argv[]) {
     
    random_device rd;
    default_random_engine gen(rd());
    
    int num_reps = 100;
    for (int rep = 0; rep < num_reps; ++rep) {
        
        Parameters params;
        
        params.simplify_window = uniform_int_distribution<int>(1, 10)(gen);
        params.max_walk_count = uniform_int_distribution<int>(1, 10)(gen);
        params.blocking_allele_size = uniform_int_distribution<int>(1, 10)(gen);
        params.max_count = uniform_int_distribution<int>(1, 10)(gen);
        params.max_num_match_pairs = uniform_int_distribution<int>(1, 10)(gen);
        params.pair_count_power = uniform_real_distribution<double>(0.001, 10.0)(gen);
        params.length_scale = uniform_int_distribution<int>(0, 1)(gen);
        params.count_penalty_threshold = uniform_real_distribution<double>(0.0, 10.0)(gen);
        params.logging_level = (logging::LoggingLevel) uniform_int_distribution<int>(0, 4)(gen);
        params.chaining_algorithm = (Anchorer::ChainAlgorithm) uniform_int_distribution<int>(0, 2)(gen);
        params.preserve_subproblems = uniform_int_distribution<int>(0, 1)(gen);
        params.subproblems_prefix = random_sequence(5, gen);
        if (uniform_int_distribution<int>(0, 1)(gen)) {
            params.subproblems_prefix = "";
        }
        params.tree_name = random_sequence(5, gen);
        if (uniform_int_distribution<int>(0, 1)(gen)) {
            params.tree_name = "";
        }
        params.all_pairs_prefix = random_sequence(5, gen);
        if (uniform_int_distribution<int>(0, 1)(gen)) {
            params.all_pairs_prefix = "";
        }
        params.fasta_name = random_sequence(5, gen);
        
        params.validate();
        
        auto config_text = params.generate_config();
        
        stringstream config_stream(config_text);
        
        Parameters reloaded(config_stream);
        
        if (params != reloaded) {
            cerr << "config did not replicate params:\n";
            cerr << config_text << '\n';
            exit(1);
        }
        
    }
    
    
    
    cerr << "passed all tests!" << endl;
}
