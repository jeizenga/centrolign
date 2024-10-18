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
    
    int num_reps = 20;
    for (int rep = 0; rep < num_reps; ++rep) {
        
        Parameters params;
        
        params.max_count = uniform_int_distribution<int>(1, 10)(gen);
        params.anchor_score_function = (ScoreFunction::AnchorScore) uniform_int_distribution<int>(0, 3)(gen);
        params.max_num_match_pairs = uniform_int_distribution<int>(1, 10)(gen);
        params.pair_count_power = uniform_real_distribution<double>(0.001, 10.0)(gen);
        params.length_intercept = uniform_real_distribution<double>(2.0, 10.0)(gen);
        params.length_decay_power = uniform_real_distribution<double>(1.0001, 10.0)(gen);
        for (size_t i = 0; i < 3; ++i) {
            params.anchor_gap_open[i] = uniform_real_distribution<double>(0.001, 10.0)(gen);
            params.anchor_gap_extend[i] = uniform_real_distribution<double>(0.001, 10.0)(gen);
        }
        sort(params.anchor_gap_open.begin(), params.anchor_gap_open.end());
        sort(params.anchor_gap_extend.begin(), params.anchor_gap_extend.end(), std::greater<double>());
        params.do_fill_in_anchoring = uniform_int_distribution<int>(0, 1)(gen);
        params.logging_level = (logging::LoggingLevel) uniform_int_distribution<int>(0, 4)(gen);
        params.chaining_algorithm = (Anchorer::ChainAlgorithm) uniform_int_distribution<int>(0, 2)(gen);
        params.constraint_method = (Partitioner::ConstraintMethod) uniform_int_distribution<int>(0, 3)(gen);
        params.minimum_segment_score = uniform_real_distribution<double>(0.001, 10.0)(gen);
        params.minimum_segment_average = uniform_real_distribution<double>(0.001, 10.0)(gen);
        params.window_length = uniform_real_distribution<double>(0.001, 10.0)(gen);
        params.generalized_length_mean = uniform_real_distribution<double>(0.001, 10.0)(gen);
        params.stitch_match = uniform_int_distribution<int>(1, 20)(gen);
        params.stitch_mismatch = uniform_int_distribution<int>(1, 20)(gen);
        for (size_t i = 0; i < 3; ++i) {
            params.stitch_gap_open[i] = uniform_int_distribution<int>(1, 20)(gen);
            params.stitch_gap_extend[i] = uniform_int_distribution<int>(1, 20)(gen);
        }
        sort(params.stitch_gap_open.begin(), params.stitch_gap_open.end());
        sort(params.stitch_gap_extend.begin(), params.stitch_gap_extend.end(), std::greater<int64_t>());
        params.max_trivial_size = uniform_int_distribution<int>(2, 20)(gen);
        params.min_wfa_size = uniform_int_distribution<int>(1, 20)(gen);
        params.max_wfa_size = params.min_wfa_size + uniform_int_distribution<int>(1, 20)(gen);
        params.max_wfa_ratio = uniform_real_distribution<double>(1.01, 20)(gen);
        params.wfa_pruning_dist = uniform_int_distribution<int>(1, 20)(gen);
        params.deletion_alignment_ratio = uniform_int_distribution<int>(1, 20)(gen);
        params.deletion_alignment_short_max_size = uniform_int_distribution<int>(1, 20)(gen);
        params.deletion_alignment_long_min_size = uniform_int_distribution<int>(1, 20)(gen);
        params.cyclize_tandem_duplications = uniform_int_distribution<int>(0, 1)(gen);
        params.max_tandem_duplication_search_rounds = uniform_int_distribution<int>(1, 20)(gen);
        params.min_cyclizing_length = uniform_int_distribution<int>(1000, 200000)(gen);
        params.preserve_subproblems = uniform_int_distribution<int>(0, 1)(gen);
        params.subproblems_prefix = random_sequence(5, gen);
        params.subalignments_filepath = random_sequence(5, gen);
        if (uniform_int_distribution<int>(0, 1)(gen)) {
            params.subproblems_prefix = "";
        }
        if (uniform_int_distribution<int>(0, 1)(gen)) {
            params.subalignments_filepath = "";
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
