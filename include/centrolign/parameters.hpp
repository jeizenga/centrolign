#ifndef centrolign_parameters_hpp
#define centrolign_parameters_hpp

#include <vector>
#include <string>

#include "centrolign/logging.hpp"
#include "centrolign/core.hpp"

namespace centrolign {

/*
 * The command line parameters and their defaults
 */
struct Parameters {
    Parameters() = default;
    Parameters(std::istream& in);
    ~Parameters() = default;
    
    // generate a config YAML that would recreate this set of parameters
    std::string generate_config() const;
    
    // check logical validity of parameters, throw an error if anything is off
    void validate() const;
    
    // pass all of the parameters to configurable modules
    void apply(Core& core) const;
        
    int64_t simplify_window = 10000;
    int64_t max_walk_count = 8;
    int64_t blocking_allele_size = 32;
    bool path_matches = true;
    int64_t max_count = 3000;
    int64_t max_num_match_pairs = 1250000;
    ScoreFunction::AnchorScore anchor_score_function = ScoreFunction::ConcaveLengthScaleInverseCount;
    double pair_count_power = 0.5;
    double length_intercept = 2250.0;
    double length_decay_power = 2.0;
    std::vector<double> anchor_gap_open = {1.25, 50.0, 5000.0};
    std::vector<double> anchor_gap_extend = {2.5, 0.1, 0.0015};
    bool do_fill_in_anchoring = true;
    logging::LoggingLevel logging_level = logging::Basic;
    Anchorer::ChainAlgorithm chaining_algorithm = Anchorer::SparseAffine;
    Partitioner::ConstraintMethod constraint_method = Partitioner::MinWindowAverage;
    double minimum_segment_score = 15000.0;
    double minimum_segment_average = 0.1;
    double window_length = 10000.0;
    double generalized_length_mean = -0.5;
    int64_t stitch_match = 20;
    int64_t stitch_mismatch = 40;
    std::vector<int64_t> stitch_gap_open = {30, 800, 2500};
    std::vector<int64_t> stitch_gap_extend = {20, 5, 1};
    int64_t max_trivial_size = 30000;
    int64_t min_wfa_size = 10000000;
    int64_t max_wfa_size = 50000000;
    double max_wfa_ratio = 1.05;
    int64_t wfa_pruning_dist = 25;
    int64_t deletion_alignment_ratio = 4;
    int64_t deletion_alignment_short_max_size = 4000;
    int64_t deletion_alignment_long_min_size = 2000;

    bool preserve_subproblems = false;
    bool skip_calibration = false;
    
    // input or output files, no default
    std::string subproblems_prefix;
    std::string subalignments_filepath;
    std::string tree_name;
    std::string all_pairs_prefix;
    std::string fasta_name;
    
    bool operator==(const Parameters& other) const;
    bool operator!=(const Parameters& other) const;
    
};

}

#endif /* centrolign_parameters_hpp */
