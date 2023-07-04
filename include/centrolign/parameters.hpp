#ifndef centrolign_parameters_hpp
#define centrolign_parameters_hpp

#include "centrolign/logging.hpp"
#include "centrolign/anchorer.hpp"
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
    int64_t max_count = 50;
    int64_t max_num_match_pairs = 1000000;
    double pair_count_power = 1.0;
    bool length_scale = true;
    double count_penalty_threshold = 16.0;
    logging::LoggingLevel logging_level = logging::Basic;
    Anchorer::ChainAlgorithm chaining_algorithm = Anchorer::SparseAffine;
    bool preserve_subproblems = false;
    
    // input or output files, no default
    std::string subproblems_prefix;
    std::string tree_name;
    std::string all_pairs_prefix;
    std::string fasta_name;
    
    bool operator==(const Parameters& other) const;
    bool operator!=(const Parameters& other) const;
};

}

#endif /* centrolign_parameters_hpp */
