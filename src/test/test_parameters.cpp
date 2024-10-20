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
        
        params.set("max_count", uniform_int_distribution<int>(1, 10)(gen));
        params.set("anchor_score_function", (ScoreFunction::AnchorScore) uniform_int_distribution<int>(0, 3)(gen));
        params.set("max_num_match_pairs", uniform_int_distribution<int>(1, 10)(gen));
        params.set("pair_count_power", uniform_real_distribution<double>(0.001, 10.0)(gen));
        params.set("length_intercept", uniform_real_distribution<double>(2.0, 10.0)(gen));
        params.set("length_decay_power", uniform_real_distribution<double>(1.0001, 10.0)(gen));
        std::array<double, 3> anchor_gap_open, anchor_gap_extend;
        for (size_t i = 0; i < 3; ++i) {
            anchor_gap_open[i] = uniform_real_distribution<double>(0.001, 10.0)(gen);
            anchor_gap_extend[i] = uniform_real_distribution<double>(0.001, 10.0)(gen);
        }
        sort(anchor_gap_open.begin(), anchor_gap_open.end());
        sort(anchor_gap_extend.begin(), anchor_gap_extend.end(), std::greater<double>());
        params.set("anchor_gap_open", anchor_gap_open);
        params.set("anchor_gap_extend", anchor_gap_extend);
        params.set("do_fill_in_anchoring", uniform_int_distribution<int>(0, 1)(gen));
        params.set("logging_level", (logging::LoggingLevel) uniform_int_distribution<int>(0, 4)(gen));
        params.set("chaining_algorithm", (Anchorer::ChainAlgorithm) uniform_int_distribution<int>(0, 2)(gen));
        params.set("constraint_method", (Partitioner::ConstraintMethod) uniform_int_distribution<int>(0, 3)(gen));
        params.set("minimum_segment_score", uniform_real_distribution<double>(0.001, 10.0)(gen));
        params.set("minimum_segment_average", uniform_real_distribution<double>(0.001, 10.0)(gen));
        params.set("window_length", uniform_real_distribution<double>(0.001, 10.0)(gen));
        params.set("generalized_length_mean", uniform_real_distribution<double>(0.001, 10.0)(gen));
        params.set("stitch_match", uniform_int_distribution<int>(1, 20)(gen));
        params.set("stitch_mismatch", uniform_int_distribution<int>(1, 20)(gen));
        std::array<int64_t, 3> stitch_gap_open, stitch_gap_extend;
        for (size_t i = 0; i < 3; ++i) {
            stitch_gap_open[i] = uniform_int_distribution<int>(1, 20)(gen);
            stitch_gap_extend[i] = uniform_int_distribution<int>(1, 20)(gen);
        }
        sort(stitch_gap_open.begin(), stitch_gap_open.end());
        sort(stitch_gap_extend.begin(), stitch_gap_extend.end(), std::greater<int64_t>());
        params.set("stitch_gap_open", stitch_gap_open);
        params.set("stitch_gap_extend", stitch_gap_extend);
        params.set("max_trivial_size", uniform_int_distribution<int>(2, 20)(gen));
        int min_wfa_size = uniform_int_distribution<int>(1, 20)(gen);
        params.set("min_wfa_size", min_wfa_size);
        params.set("max_wfa_size", min_wfa_size + uniform_int_distribution<int>(1, 20)(gen));
        params.set("max_wfa_ratio", uniform_real_distribution<double>(1.01, 20)(gen));
        params.set("wfa_pruning_dist", uniform_int_distribution<int>(1, 20)(gen));
        params.set("deletion_alignment_ratio", uniform_int_distribution<int>(1, 20)(gen));
        params.set("deletion_alignment_short_max_size", uniform_int_distribution<int>(1, 20)(gen));
        params.set("deletion_alignment_long_min_size", uniform_int_distribution<int>(1, 20)(gen));
        params.set("cyclize_tandem_duplications", uniform_int_distribution<int>(0, 1)(gen));
        params.set("max_tandem_duplication_search_rounds", uniform_int_distribution<int>(1, 20)(gen));
        params.set("min_cyclizing_length", uniform_int_distribution<int>(1000, 200000)(gen));
        params.set("preserve_subproblems", uniform_int_distribution<int>(0, 1)(gen));
        params.set("subproblems_prefix", random_sequence(5, gen));
        params.set("subalignments_filepath", random_sequence(5, gen));
        if (uniform_int_distribution<int>(0, 1)(gen)) {
            params.set("subproblems_prefix", "");
        }
        if (uniform_int_distribution<int>(0, 1)(gen)) {
            params.set("subalignments_filepath", "");
        }
        params.set("tree_name", random_sequence(5, gen));
        if (uniform_int_distribution<int>(0, 1)(gen)) {
            params.set("tree_name", "");
        }
        params.set("all_pairs_prefix", random_sequence(5, gen));
        if (uniform_int_distribution<int>(0, 1)(gen)) {
            params.set("all_pairs_prefix", "");
        }
        params.set("fasta_name", random_sequence(5, gen));
        
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
