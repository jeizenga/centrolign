#include "centrolign/parameters.hpp"

#include "centrolign/utility.hpp"

#include <stdexcept>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <limits>

namespace centrolign {

using namespace std;

inline std::string string_or_null(const std::string& str) {
    return str.empty() ? "\"\"" : str;
}

// TODO: repetitive
inline std::vector<double> parse_list_double(const std::string& value) {
    std::vector<double> values;
    for (auto& token : tokenize(value, ',')) {
        values.push_back(parse_double(token));
    }
    return values;
}
inline std::vector<int64_t> parse_list_int(const std::string& value) {
    std::vector<int64_t> values;
    for (auto& token : tokenize(value, ',')) {
        values.push_back(parse_int(token));
    }
    return values;
}

template<typename T>
inline std::string list_to_string(const std::vector<T>& values) {
    std::stringstream strm;
    strm << std::setprecision(17);
    for (size_t i = 0; i < values.size(); ++i) {
        if (i) {
            strm << ',';
        }
        strm << values[i];
    }
    return strm.str();
}

Parameters::Parameters(std::istream& in) {
    
    string config;
    {
        stringstream strm;
        strm << in.rdbuf();
        config = std::move(strm.str());
    }
    
    auto start = config.find("---");
    if (start >= config.size()) {
        throw runtime_error("Config file is missing '---' delimiter");
    }
    
    config = config.substr(start + 3, config.size());
    
    stringstream scan_strm(config);
    string line;
    while (scan_strm) {
        getline(scan_strm, line);
        
        if (all_of(line.begin(), line.end(), [](char c) { return isspace(c); })) {
            continue;
        }
        
        auto colon_pos = line.find(':');
        if (colon_pos >= line.size()) {
            throw runtime_error("Config has line missing ':' delimiter in line: " + line);
        }
        int name_start = 0, name_end = colon_pos, value_start = colon_pos + 1, value_end = line.size();
        
        while (name_start < name_end && isspace(line[name_start])) {
            ++name_start;
        }
        while (name_start < name_end && isspace(line[name_end - 1])) {
            --name_end;
        }
        while (value_start < value_end && isspace(line[value_start])) {
            ++value_start;
        }
        while (value_start < value_end && isspace(line[value_end - 1])) {
            --value_end;
        }
        if (name_start >= name_end) {
            throw runtime_error("Config has line with no variable name in line: " + line);
        }
        if (value_start >= value_end) {
            throw runtime_error("Config has line with no variable value in line: " + line);
        }
        if (line[value_start] == '"' &&
            (value_end == value_start + 1 || line[value_end - 1] != '"')) {
            throw runtime_error("Unmatched \" in config file line: " + line);
        }
        
        if (line[value_start] == '"') {
            ++value_start;
            --value_end;
        }
        
        string name = line.substr(name_start, name_end - name_start);
        string value = line.substr(value_start, value_end - value_start);
        
        if (name == "simplify_window") {
            simplify_window = parse_int(value);
        }
        else if (name == "max_walk_count") {
            max_walk_count = parse_int(value);
        }
        else if (name == "blocking_allele_size") {
            blocking_allele_size = parse_int(value);
        }
        else if (name == "path_matches") {
            path_matches = parse_bool(value);
        }
        else if (name == "use_color_set_size") {
            use_color_set_size = parse_bool(value);
        }
        else if (name == "max_count") {
            max_count = parse_int(value);
        }
        else if (name == "max_num_match_pairs") {
            max_num_match_pairs = parse_int(value);
        }
        else if (name == "anchor_score_function") {
            anchor_score_function = (ScoreFunction::AnchorScore) parse_int(value);
        }
        else if (name == "pair_count_power") {
            pair_count_power = parse_double(value);
        }
        else if (name == "length_intercept") {
            length_intercept = parse_double(value);
        }
        else if (name == "length_decay_power") {
            length_decay_power = parse_double(value);
        }
        else if (name == "anchor_gap_open") {
            anchor_gap_open = parse_list_double(value);
        }
        else if (name == "anchor_gap_extend") {
            anchor_gap_extend = parse_list_double(value);
        }
        else if (name == "do_fill_in_anchoring") {
            do_fill_in_anchoring = parse_bool(value);
        }
        else if (name == "logging_level") {
            logging_level = (logging::LoggingLevel) parse_int(value);
        }
        else if (name == "chaining_algorithm") {
            chaining_algorithm = (Anchorer::ChainAlgorithm) parse_int(value);
        }
        else if (name == "constraint_method") {
            constraint_method = (Partitioner::ConstraintMethod) parse_int(value);
        }
        else if (name == "minimum_segment_score") {
            minimum_segment_score = parse_double(value);
        }
        else if (name == "minimum_segment_average") {
            minimum_segment_average = parse_double(value);
        }
        else if (name == "window_length") {
            window_length = parse_double(value);
        }
        else if (name == "generalized_length_mean") {
            generalized_length_mean = parse_double(value);
        }
        else if (name == "stitch_match") {
            stitch_match = parse_int(value);
        }
        else if (name == "stitch_mismatch") {
            stitch_mismatch = parse_int(value);
        }
        else if (name == "stitch_gap_open") {
            stitch_gap_open = parse_list_int(value);
        }
        else if (name == "stitch_gap_extend") {
            stitch_gap_extend = parse_list_int(value);
        }
        else if (name == "max_trivial_size") {
            max_trivial_size = parse_int(value);
        }
        else if (name == "min_wfa_size") {
            min_wfa_size = parse_int(value);
        }
        else if (name == "max_wfa_size") {
            max_wfa_size = parse_int(value);
        }
        else if (name == "max_wfa_ratio") {
            max_wfa_ratio = parse_double(value);
        }
        else if (name == "wfa_pruning_dist") {
            wfa_pruning_dist = parse_int(value);
        }
        else if (name == "deletion_alignment_ratio") {
            deletion_alignment_ratio = parse_int(value);
        }
        else if (name == "deletion_alignment_short_max_size") {
            deletion_alignment_short_max_size = parse_int(value);
        }
        else if (name == "deletion_alignment_long_min_size") {
            deletion_alignment_long_min_size = parse_int(value);
        }
        else if (name == "cyclize_tandem_duplications") {
            cyclize_tandem_duplications = parse_int(value);
        }
        else if (name == "max_tandem_duplication_search_rounds") {
            max_tandem_duplication_search_rounds = parse_int(value);
        }
        else if (name == "min_cyclizing_length") {
            min_cyclizing_length = parse_int(value);
        }
        else if (name == "preserve_subproblems") {
            preserve_subproblems = parse_bool(value);
        }
        else if (name == "skip_calibration") {
            skip_calibration = parse_bool(value);
        }
        else if (name == "subproblems_prefix") {
            subproblems_prefix = value;
        }
        else if (name == "subalignments_filepath") {
            subalignments_filepath = value;
        }
        else if (name == "tree_name") {
            tree_name = value;
        }
        else if (name == "all_pairs_prefix") {
            all_pairs_prefix = value;
        }
        else if (name == "fasta_name") {
            fasta_name = value;
        }
        else {
            throw runtime_error("Config contains unrecognized variable: " + name);
        }
    }
}
 
string Parameters::generate_config() const {
    
    stringstream strm;
    strm << setprecision(17); // the number of significant digits in a double
    strm << "---\n";
    strm << ' ' << "simplify_window" << ": " << simplify_window << '\n';
    strm << ' ' << "max_walk_count" << ": " << max_walk_count << '\n';
    strm << ' ' << "blocking_allele_size" << ": " << blocking_allele_size << '\n';
    strm << ' ' << "path_matches" << ": " << path_matches << '\n';
    strm << ' ' << "use_color_set_size" << ": " << use_color_set_size << '\n';
    strm << ' ' << "max_count" << ": " << max_count << '\n';
    strm << ' ' << "max_num_match_pairs" << ": " << max_num_match_pairs << '\n';
    strm << ' ' << "anchor_score_function" << ": " << (int) anchor_score_function << '\n';
    strm << ' ' << "pair_count_power" << ": " << pair_count_power << '\n';
    strm << ' ' << "length_intercept" << ": " << length_intercept << '\n';
    strm << ' ' << "length_decay_power" << ": " << length_decay_power << '\n';
    strm << ' ' << "anchor_gap_open" << ": " << list_to_string(anchor_gap_open) << '\n';
    strm << ' ' << "anchor_gap_extend" << ": " << list_to_string(anchor_gap_extend) << '\n';
    strm << ' ' << "do_fill_in_anchoring" << ": " << do_fill_in_anchoring << '\n';
    strm << ' ' << "logging_level" << ": " << (int) logging_level << '\n';
    strm << ' ' << "chaining_algorithm" << ": " << (int) chaining_algorithm << '\n';
    strm << ' ' << "constraint_method" << ": " << (int) constraint_method << '\n';
    strm << ' ' << "minimum_segment_score" << ": " << minimum_segment_score << '\n';
    strm << ' ' << "minimum_segment_average" << ": " << minimum_segment_average << '\n';
    strm << ' ' << "window_length" << ": " << window_length << '\n';
    strm << ' ' << "generalized_length_mean" << ": " << generalized_length_mean << '\n';
    strm << ' ' << "stitch_match" << ": " << stitch_match << '\n';
    strm << ' ' << "stitch_mismatch" << ": " << stitch_mismatch << '\n';
    strm << ' ' << "stitch_gap_open" << ": " << list_to_string(stitch_gap_open) << '\n';
    strm << ' ' << "stitch_gap_extend" << ": " << list_to_string(stitch_gap_extend) << '\n';
    strm << ' ' << "max_trivial_size" << ": " << max_trivial_size << '\n';
    strm << ' ' << "min_wfa_size" << ": " << min_wfa_size << '\n';
    strm << ' ' << "max_wfa_size" << ": " << max_wfa_size << '\n';
    strm << ' ' << "max_wfa_ratio" << ": " << max_wfa_ratio << '\n';
    strm << ' ' << "wfa_pruning_dist" << ": " << wfa_pruning_dist << '\n';
    strm << ' ' << "deletion_alignment_ratio" << ": " << deletion_alignment_ratio << '\n';
    strm << ' ' << "deletion_alignment_short_max_size" << ": " << deletion_alignment_short_max_size << '\n';
    strm << ' ' << "deletion_alignment_long_min_size" << ": " << deletion_alignment_long_min_size << '\n';
    strm << ' ' << "cyclize_tandem_duplications" << ": " << cyclize_tandem_duplications << '\n';
    strm << ' ' << "max_tandem_duplication_search_rounds" << ": " << max_tandem_duplication_search_rounds << '\n';
    strm << ' ' << "min_cyclizing_length" << ": " << min_cyclizing_length << '\n';
    strm << ' ' << "preserve_subproblems" << ": " << preserve_subproblems << '\n';
    strm << ' ' << "skip_calibration" << ": " << skip_calibration << '\n';
    strm << ' ' << "subproblems_prefix" << ": " << string_or_null(subproblems_prefix) << '\n';
    strm << ' ' << "subalignments_filepath" << ": " << string_or_null(subalignments_filepath) << '\n';
    strm << ' ' << "tree_name" << ": " << string_or_null(tree_name) << '\n';
    strm << ' ' << "all_pairs_prefix" << ": " << string_or_null(all_pairs_prefix) << '\n';
    strm << ' ' << "fasta_name" << ": " << string_or_null(fasta_name) << '\n';
    
    return strm.str();
}


void Parameters::validate() const {
    if (simplify_window < 0) {
        throw runtime_error("Got negative value " + to_string(simplify_window) + " for simplification window");
    }
    if (max_walk_count < 0) {
        throw runtime_error("Got negative value " + to_string(max_walk_count) + " for simplification maximum walk count");
    }
    if (blocking_allele_size < 0) {
        throw runtime_error("Got negative value " + to_string(blocking_allele_size) + " for simplification blocking allele size");
    }
    if (max_count <= 0) {
        throw runtime_error("Got non-positive value " + to_string(max_count) + " for maximum anchor occurrence count");
    }
    if (max_num_match_pairs <= 0) {
        throw runtime_error("Got non-positive value " + to_string(max_num_match_pairs) + " for maximum number of anchor matches");
    }
    if (anchor_score_function < 0 || anchor_score_function > 3) {
        throw runtime_error("Got invalid value " + to_string(anchor_score_function) + " for anchor scoring function");
    }
    if (pair_count_power < 0.0) {
        throw runtime_error("Got negative value " + to_string(pair_count_power) + " for match count scoring power");
    }
    if (length_intercept < 1.0) {
        throw runtime_error("Got value " + to_string(length_intercept) + " below 1 for length intercept / maximum positive-scoring length");
    }
    if (length_decay_power < 0.0) {
        throw runtime_error("Got negative value " + to_string(length_decay_power) + " for length decay power");
    }
    if (anchor_gap_open.size() != 3 || anchor_gap_extend.size() != 3) {
        throw runtime_error("Anchor gap open/extend penalties must have length 3, got " + list_to_string(anchor_gap_open) + " and " + list_to_string(anchor_gap_extend));
    }
    for (size_t i = 1; i < anchor_gap_open.size(); ++i) {
        if (anchor_gap_open[i - 1] > anchor_gap_open[i] || anchor_gap_extend[i - 1] < anchor_gap_extend[i]) {
            throw runtime_error("Anchor gap open penalities must be provided in increasing order and gap extend penalties must be provided in decreasing order, got " + list_to_string(anchor_gap_open) + " and " + list_to_string(anchor_gap_extend));
        }
    }
    if (logging_level < 0 || logging_level > 4) {
        throw runtime_error("Got invalid value " + to_string(logging_level) + " for logging level");
    }
    if (chaining_algorithm < 0 || chaining_algorithm > 2) {
        throw runtime_error("Got invalid value " + to_string(chaining_algorithm) + " for chaining algorithm");
    }
    if (constraint_method < 0 || constraint_method > 3) {
        throw runtime_error("Got invalid value " + to_string(constraint_method) + " for constraint method");
    }
    if (minimum_segment_score < 0.0) {
        throw runtime_error("Got negative value " + to_string(minimum_segment_score) + " for minimum segment score");
    }
    if (minimum_segment_average < 0.0) {
        throw runtime_error("Got negative value " + to_string(minimum_segment_average) + " for minimum segment windowed-average score");
    }
    if (window_length <= 0.0) {
        throw runtime_error("Got non-positive value " + to_string(window_length) + " for windowed-average window length");
    }
    if (stitch_gap_open.size() != 3 || stitch_gap_extend.size() != 3) {
        throw runtime_error("Stitching gap open/extend penalties must have length 3, got " + list_to_string(stitch_gap_open) + " and " + list_to_string(stitch_gap_extend));
    }
    for (size_t i = 1; i < stitch_gap_open.size(); ++i) {
        if (stitch_gap_open[i - 1] > stitch_gap_open[i] || stitch_gap_extend[i - 1] < stitch_gap_extend[i]) {
            throw runtime_error("Stitching gap open penalities must be provided in increasing order and gap extend penalties must be provided in decreasing order, got " + list_to_string(stitch_gap_open) + " and " + list_to_string(stitch_gap_extend));
        }
    }
    int64_t stitch_max = numeric_limits<uint32_t>::max();
    if (stitch_match < 0) {
        throw runtime_error("Got negative stitching match bonus " + to_string(stitch_match));
    }
    if (stitch_match > stitch_max) {
        throw runtime_error("Stitching match bonus must be less than " + to_string(stitch_max) + ", got " + to_string(stitch_match));
    }
    if (stitch_mismatch < 0) {
        throw runtime_error("Got negative stitching mismatch penalty " + to_string(stitch_mismatch));
    }
    if (stitch_mismatch > stitch_max) {
        throw runtime_error("Stitching mismatch penalty must be less than " + to_string(stitch_max) + ", got " + to_string(stitch_mismatch));
    }
    for (size_t i = 0; i < stitch_gap_open.size(); ++i) {
        if (stitch_gap_open[i] < 0) {
            throw runtime_error("Got negative stitching gap open penalty " + to_string(stitch_gap_open[i]));
        }
        if (stitch_gap_open[i] > stitch_max) {
            throw runtime_error("Stitching gap open penalty must be less than " + to_string(stitch_max) + ", got " + to_string(stitch_gap_open[i]));
        }
        if (stitch_gap_extend[i] < 0) {
            throw runtime_error("Got negative stitching gap extend penalty " + to_string(stitch_gap_extend[i]));
        }
        if (stitch_gap_extend[i] > stitch_max) {
            throw runtime_error("Stitching gap extend penalty must be less than " + to_string(stitch_max) + ", got " + to_string(stitch_gap_extend[i]));
        }
    }
    if (max_trivial_size < 0) {
        throw runtime_error("Got negative maximum size for a trivial DP matrix " + to_string(max_trivial_size));
    }
    if (min_wfa_size < 0) {
        throw runtime_error("Got negative mininmum size of DP matrix for performing WFA " + to_string(min_wfa_size));
    }
    if (max_wfa_size < 0) {
        throw runtime_error("Got negative maximum size of DP matrix for performing WFA " + to_string(max_wfa_size));
    }
    if (max_wfa_size < min_wfa_size) {
        throw runtime_error("Got greater minimum than maximum for size of DP matrix for performing WFA: " + to_string(min_wfa_size) + ", " + to_string(max_wfa_size));
    }
    if (max_wfa_ratio < 1.0) {
        throw runtime_error("Got maximum ratio of longer:shorter graph size to perform WFA that is less than 1.0: " + to_string(max_wfa_ratio));
    }
    if (wfa_pruning_dist < 0) {
        throw runtime_error("Got negative pruning distance for WFA " + to_string(wfa_pruning_dist));
    }
    if (deletion_alignment_ratio < 1) {
        throw runtime_error("Got minimum ratio of graph sizes to be considered a deletion that is less than 1: " + to_string(deletion_alignment_ratio));
    }
    if (deletion_alignment_short_max_size < 0) {
        throw runtime_error("Got negative maximum size for shorter graph to be considered a deletion: " + to_string(deletion_alignment_short_max_size));
    }
    if (deletion_alignment_long_min_size < 0) {
        throw runtime_error("Got negative minimum size for longer graph to be considered a deletion: " + to_string(deletion_alignment_long_min_size));
    }
    if (cyclize_tandem_duplications && max_tandem_duplication_search_rounds == 0) {
        throw runtime_error("Cannot cyclize tandem duplications with 0 rounds of tandem duplication search");
    }
    if (max_tandem_duplication_search_rounds < 0) {
        throw runtime_error("Got number of tandem duplication search rounds that is negative: " + to_string(max_tandem_duplication_search_rounds));
    }
    if (min_cyclizing_length < 0) {
        throw runtime_error("Got minimum cyclizing tandem duplication size that is negative: " + to_string(min_cyclizing_length));
    }
    if (fasta_name.empty()) {
        throw runtime_error("FASTA input is missing");
    }
}


void Parameters::apply(Core& core) const {
    
    // pass through parameters
    core.simplifier.min_dist_window = simplify_window;
    core.simplifier.preserve_bubble_size = blocking_allele_size;
    core.simplifier.max_walks = max_walk_count;
    
    core.match_finder.path_matches = path_matches;
    core.match_finder.use_color_set_size = use_color_set_size;
    core.match_finder.max_count = max_count;
    
    core.score_function.anchor_score_function = anchor_score_function;
    core.score_function.pair_count_power = pair_count_power;
    core.score_function.length_intercept = length_intercept;
    core.score_function.length_decay_power = length_decay_power;
    
    core.anchorer.chaining_algorithm = chaining_algorithm;
    core.anchorer.max_num_match_pairs = max_num_match_pairs;
    for (size_t i = 0; i < core.anchorer.gap_open.size(); ++i) {
        core.anchorer.gap_open[i] = anchor_gap_open[i];
        core.anchorer.gap_extend[i] = anchor_gap_extend[i];
    }
    core.anchorer.do_fill_in_anchoring = do_fill_in_anchoring;
    
    core.partitioner.constraint_method = constraint_method;
    core.partitioner.minimum_segment_score = minimum_segment_score;
    core.partitioner.minimum_segment_average = minimum_segment_average;
    core.partitioner.window_length = window_length;
    core.partitioner.generalized_length_mean = generalized_length_mean;
    
    core.stitcher.alignment_params.match = stitch_match;
    core.stitcher.alignment_params.mismatch = stitch_mismatch;
    for (size_t i = 0; i < core.stitcher.alignment_params.gap_open.size(); ++i) {
        core.stitcher.alignment_params.gap_open[i] = stitch_gap_open[i];
        core.stitcher.alignment_params.gap_extend[i] = stitch_gap_extend[i];
    }
    core.stitcher.max_trivial_size = max_trivial_size;
    core.stitcher.min_wfa_size = min_wfa_size;
    core.stitcher.max_wfa_size = max_wfa_size;
    core.stitcher.max_wfa_ratio = max_wfa_ratio;
    core.stitcher.wfa_pruning_dist = wfa_pruning_dist;
    core.stitcher.deletion_alignment_ratio = deletion_alignment_ratio;
    core.stitcher.deletion_alignment_short_max_size = deletion_alignment_short_max_size;
    core.stitcher.deletion_alignment_long_min_size = deletion_alignment_long_min_size;
    
    core.bonder.min_length = min_cyclizing_length;
    
    core.cyclize_tandem_duplications = cyclize_tandem_duplications;
    core.max_tandem_duplication_search_rounds = max_tandem_duplication_search_rounds;
    core.preserve_subproblems = preserve_subproblems;
    core.skip_calibration = skip_calibration;
    
    core.subproblems_prefix = subproblems_prefix;
    core.subalignments_filepath = subalignments_filepath;
}

bool Parameters::operator==(const Parameters& other) const {
    
    return (simplify_window == other.simplify_window &&
            max_walk_count == other.max_walk_count &&
            blocking_allele_size == other.blocking_allele_size &&
            path_matches == other.path_matches &&
            use_color_set_size == other.use_color_set_size &&
            max_count == other.max_count &&
            max_num_match_pairs == other.max_num_match_pairs &&
            chaining_algorithm == other.chaining_algorithm &&
            constraint_method == other.constraint_method &&
            minimum_segment_score == other.minimum_segment_score &&
            minimum_segment_average == other.minimum_segment_average &&
            window_length == other.window_length &&
            generalized_length_mean == other.generalized_length_mean &&
            stitch_match == other.stitch_match &&
            stitch_mismatch == other.stitch_mismatch &&
            stitch_gap_open == other.stitch_gap_open &&
            stitch_gap_extend == other.stitch_gap_extend &&
            max_trivial_size == other.max_trivial_size &&
            min_wfa_size == other.min_wfa_size &&
            max_wfa_size == other.max_wfa_size &&
            max_wfa_ratio == other.max_wfa_ratio &&
            wfa_pruning_dist == other.wfa_pruning_dist &&
            deletion_alignment_ratio == other.deletion_alignment_ratio &&
            deletion_alignment_short_max_size == other.deletion_alignment_short_max_size &&
            deletion_alignment_long_min_size == other.deletion_alignment_long_min_size &&
            pair_count_power == other.pair_count_power &&
            length_intercept == other.length_intercept &&
            length_decay_power == other.length_decay_power &&
            anchor_gap_open == other.anchor_gap_open &&
            anchor_gap_extend == other.anchor_gap_extend &&
            do_fill_in_anchoring == other.do_fill_in_anchoring &&
            logging_level == other.logging_level &&
            preserve_subproblems == other.preserve_subproblems &&
            skip_calibration == other.skip_calibration &&
            subproblems_prefix == other.subproblems_prefix &&
            subalignments_filepath == other.subalignments_filepath &&
            tree_name == other.tree_name &&
            all_pairs_prefix == other.all_pairs_prefix &&
            fasta_name == other.fasta_name);
}

bool Parameters::operator!=(const Parameters& other) const {
    return !(*this == other);
}

}
