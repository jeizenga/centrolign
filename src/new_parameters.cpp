#include "centrolign/new_parameters.hpp"
#include "centrolign/utility.hpp"


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

NewParameters::NewParameters() {
    
    initialize_submodule(IO, "Parameters related to file I/O and logging");
    initialize_submodule(MatchFinding, "Parameters related to identifying matches between graphs");
    initialize_submodule(Anchoring, "Parameters related to identifying high-scoring chains of matches to anchor alignment");
    initialize_submodule(IdentifyingAlignability, "Parameters related to determining whether a graph region is alignable");
    initialize_submodule(Aligning, "Parameters related to constructing a base-level alignment");
    initialize_submodule(InducingCycles, "Parameters related to inducing cycles at tandem duplications");
    initialize_submodule(DeveloperTools, "Parameters that were designed only to facilitate developers");
    
    add_parameter(IO, "fasta_name", String, std::string(), "The path to a FASTA file containing all of input sequences");
    add_parameter(IO, "tree_name", String, std::string(), "The path to a guide tree for the alignment in Newick format (sample names must match sequence names from the input FASTA)");
    add_parameter(IO, "logging_level", Enum, logging::Basic, "The path to a guide tree for the alignment in Newick format (sample names must match sequence names from the input FASTA)");
    add_parameter(IO, "all_pairs_prefix", String, std::string(), "If provided, save the induced pairwise alignment for each pair of sequences in CIGAR format to files with this prefix");
    add_parameter(IO, "subproblems_prefix", String, std::string(), "If provided, save the results of the intermediate subproblems in GFA format to files with this prefix");
    add_parameter(IO, "subalignments_filepath", String, std::string(), "If provided, save the path-to-path alignment from each subproblem to files with this prefix");
    
    add_parameter(MatchFinding, "max_count", Integer, 3000, "Only query matches that occur at most this many times on either of the two graphs");
    add_parameter(MatchFinding, "use_color_set_size", Bool, true, "Use Hui's (1992) color set size index instead of a merge sort tree (CSS is generally faster)");
    
    add_parameter(Anchoring, "max_num_match_pairs", Integer, 1250000, "The maximum number of matches between two graphs that will be considered during chaining");
    add_parameter(Anchoring, "do_fill_in_anchoring", Bool, true, "Attempt to fill in the anchor chain using matches that were not considered due to the limit on the maximum number of matches");
    add_parameter(Anchoring, "chaining_algorithm", Enum, Anchorer::SparseAffine, "The chaining algorithm used:\n"
                  "- " + std::to_string((int) Anchorer::Exhaustive) + ": Simple exhaustive algorithm (slow)\n"
                  "- " + std::to_string((int) Anchorer::Sparse) + ": Sparse algorithm with free gaps\n"
                  "- " + std::to_string((int) Anchorer::SparseAffine) + ": Sparse algorithm with affine gap penalties");
    add_parameter(Anchoring, "anchor_gap_open", DoubleArray3, std::array<double, 3>{1.25, 50.0, 5000.0}, "The gap open penalties used for anchoring with affine gap penalties");
    add_parameter(Anchoring, "anchor_gap_extend", DoubleArray3, std::array<double, 3>{2.5, 0.1, 0.0015}, "The gap extend penalties used for anchoring with affine gap penalties");
    add_parameter(Anchoring, "anchor_score_function", Enum, ScoreFunction::ConcaveLengthScaleInverseCount, "The scoring function used to prioritize anchors during chaining:\n"
                  "- " + std::to_string((int) ScoreFunction::InverseCount) + ": Inverse of count\n"
                  "- " + std::to_string((int) ScoreFunction::LengthScaleInverseCount) + ": Length of match scaled by inverse of count\n"
                  "- " + std::to_string((int) ScoreFunction::ConcaveLengthScaleInverseCount) + ": Length scaled by inverse of count with a subtracted convex monomial term based on length\n"
                  "- " + std::to_string((int) ScoreFunction::ConcaveLengthScaleCountDifference) + ": Length with a subtracted convex monomial term based on length and count");
    add_parameter(Anchoring, "pair_count_power", Double, 0.5, "The power that the count is raised to when used as an inverse factor to the anchor scoring function");
    add_parameter(Anchoring, "length_intercept", Double, 2250.0, "When using an anchoring scoring function with a convex subtracted term, the longest possible postively-scoring match");
    add_parameter(Anchoring, "length_decay_power", Double, 2.0, "When using an anchoring scoring function with a convex subtracted term, the power of the subtracted monomial");
    
    add_parameter(IdentifyingAlignability, "constraint_method", Enum, Partitioner::MinWindowAverage, "The method used to partition the anchor chain into alignable and unalignable regions:\n"
                  "- " + std::to_string((int) Partitioner::Null) + ": Do not attempt to partition; consider all sequences alignable\n"
                  "- " + std::to_string((int) Partitioner::Unconstrained) + ": Choose the highest scoring set of anchors\n"
                  "- " + std::to_string((int) Partitioner::MinAverage) + ": Choose the highest scoring set of anchors, with each alignable segment having limit on its average value (score/length)\n"
                  "- " + std::to_string((int) Partitioner::MinWindowAverage) + ": Choose the highest scoring set of anchors, with each alignable segment having limit on a windowed average value (score/window size) across all windows inside the segment");
    add_parameter(IdentifyingAlignability, "minimum_segment_score", Double, 15000.0, "The minimum total score that an alignable segment must have");
    add_parameter(IdentifyingAlignability, "minimum_segment_average", Double, 0.1, "The minimum average score that an alignable segment must have");
    add_parameter(IdentifyingAlignability, "window_length", Double, 10000.0, "The length of the window used in the windowed average");
    add_parameter(IdentifyingAlignability, "generalized_length_mean", Double, -0.5, "Parameter of the Holder mean used to combine length on the two graphs into a single length measurement");
    
    add_parameter(Aligning, "stitch_match", Integer, 20, "Match value when stitching anchors into a base-level alignment");
    add_parameter(Aligning, "stitch_mismatch", Integer, 80, "Mismatch penalty when stitching anchors into a base-level alignment");
    add_parameter(Aligning, "stitch_gap_open", IntegerArray3, std::array<int64_t, 3>{60, 800, 2500}, "Piecewise affine gap open penalties when stitching anchors into a base-level alignment");
    add_parameter(Aligning, "stitch_gap_extend", IntegerArray3, std::array<int64_t, 3>{30, 5, 1}, "Piecewise affine gap extend penalties when stitching anchors into a base-level alignment");
    add_parameter(Aligning, "max_trivial_size", Integer, 30000, "Maximum size of a dynamic programming matrix that will be aligned even if was identified as unalignable");
    add_parameter(Aligning, "min_wfa_size", Integer, 40000000, "Minimum size of a dynamic programming matrix that will be aligned using graph-graph WFA");
    add_parameter(Aligning, "max_wfa_size", Integer, 75000000, "Maximum size of a dynamic programming matrix that will be aligned using graph-graph WFA");
    add_parameter(Aligning, "max_wfa_ratio", Double, 1.05, "Maximum ratio of long-to-short side of the dynamic programming matrix for graph-graph WFA to be used");
    add_parameter(Aligning, "wfa_pruning_dist", Integer, 25, "The lagging distance for a diagonal to be pruned in graph-graph WFA");
    add_parameter(Aligning, "deletion_alignment_ratio", Integer, 8, "The minimum ratio of long-to-short side of the dynamic programming matrix to use WFA-based implied deletion algorithm");
    add_parameter(Aligning, "deletion_alignment_short_max_size", Integer, 1500, "The maximum size of the short side of the dynamic programming matrix to use WFA-based implied deletion algorithm");
    add_parameter(Aligning, "deletion_alignment_long_min_size", Integer, 2000, "The minimum size of the long side of the dynamic programming matrix to use WFA-based implied deletion algorithm");
    
    add_parameter(InducingCycles, "cyclize_tandem_duplications", Bool, false, "Identify tandem duplications in the sequences and use them to induce cycles in the final graph");
    add_parameter(InducingCycles, "max_tandem_duplication_search_rounds", Integer, 3, "The maximum number of nested tandem duplications to attempt finding for any given subsequence");
    add_parameter(InducingCycles, "min_cyclizing_length", Integer, 100000, "The maximum number of nested tandem duplications to attempt finding for any given subsequence");
    
    add_parameter(DeveloperTools, "bonds_filepath", String, std::string(), "If provided, save the alignments of all tandem duplications identified in the cyclization process to files with this prefix");
    add_parameter(DeveloperTools, "preserve_subproblems", Bool, false, "Do not clear out data from completed subproblems as the algorithm goes");
    add_parameter(DeveloperTools, "skip_calibration", Bool, false, "Do not calibrate the scoring parameters to the sequence repetitiveness");
    
}


//NewParameters::NewParameters(std::istream& in) {
//
//    string config;
//    {
//        stringstream strm;
//        strm << in.rdbuf();
//        config = std::move(strm.str());
//    }
//
//    auto start = config.find("---");
//    if (start >= config.size()) {
//        throw runtime_error("Config file is missing '---' delimiter");
//    }
//
//    config = config.substr(start + 3, config.size());
//
//    stringstream scan_strm(config);
//    string line;
//    while (scan_strm) {
//        getline(scan_strm, line);
//
//        if (all_of(line.begin(), line.end(), [](char c) { return isspace(c); })) {
//            continue;
//        }
//
//        auto colon_pos = line.find(':');
//        if (colon_pos >= line.size()) {
//            throw runtime_error("Config has line missing ':' delimiter in line: " + line);
//        }
//        int name_start = 0, name_end = colon_pos, value_start = colon_pos + 1, value_end = line.size();
//
//        while (name_start < name_end && isspace(line[name_start])) {
//            ++name_start;
//        }
//        while (name_start < name_end && isspace(line[name_end - 1])) {
//            --name_end;
//        }
//        while (value_start < value_end && isspace(line[value_start])) {
//            ++value_start;
//        }
//        while (value_start < value_end && isspace(line[value_end - 1])) {
//            --value_end;
//        }
//        if (name_start >= name_end) {
//            throw runtime_error("Config has line with no variable name in line: " + line);
//        }
//        if (value_start >= value_end) {
//            throw runtime_error("Config has line with no variable value in line: " + line);
//        }
//        if (line[value_start] == '"' &&
//            (value_end == value_start + 1 || line[value_end - 1] != '"')) {
//            throw runtime_error("Unmatched \" in config file line: " + line);
//        }
//
//        if (line[value_start] == '"') {
//            ++value_start;
//            --value_end;
//        }
//
//        string name = line.substr(name_start, name_end - name_start);
//        string value = line.substr(value_start, value_end - value_start);
//
//        if (name == "use_color_set_size") {
//            use_color_set_size = parse_bool(value);
//        }
//        else if (name == "max_count") {
//            max_count = parse_int(value);
//        }
//        else if (name == "max_num_match_pairs") {
//            max_num_match_pairs = parse_int(value);
//        }
//        else if (name == "anchor_score_function") {
//            anchor_score_function = (ScoreFunction::AnchorScore) parse_int(value);
//        }
//        else if (name == "pair_count_power") {
//            pair_count_power = parse_double(value);
//        }
//        else if (name == "length_intercept") {
//            length_intercept = parse_double(value);
//        }
//        else if (name == "length_decay_power") {
//            length_decay_power = parse_double(value);
//        }
//        else if (name == "anchor_gap_open") {
//            anchor_gap_open = parse_list_double(value);
//        }
//        else if (name == "anchor_gap_extend") {
//            anchor_gap_extend = parse_list_double(value);
//        }
//        else if (name == "do_fill_in_anchoring") {
//            do_fill_in_anchoring = parse_bool(value);
//        }
//        else if (name == "logging_level") {
//            logging_level = (logging::LoggingLevel) parse_int(value);
//        }
//        else if (name == "chaining_algorithm") {
//            chaining_algorithm = (Anchorer::ChainAlgorithm) parse_int(value);
//        }
//        else if (name == "constraint_method") {
//            constraint_method = (Partitioner::ConstraintMethod) parse_int(value);
//        }
//        else if (name == "minimum_segment_score") {
//            minimum_segment_score = parse_double(value);
//        }
//        else if (name == "minimum_segment_average") {
//            minimum_segment_average = parse_double(value);
//        }
//        else if (name == "window_length") {
//            window_length = parse_double(value);
//        }
//        else if (name == "generalized_length_mean") {
//            generalized_length_mean = parse_double(value);
//        }
//        else if (name == "stitch_match") {
//            stitch_match = parse_int(value);
//        }
//        else if (name == "stitch_mismatch") {
//            stitch_mismatch = parse_int(value);
//        }
//        else if (name == "stitch_gap_open") {
//            stitch_gap_open = parse_list_int(value);
//        }
//        else if (name == "stitch_gap_extend") {
//            stitch_gap_extend = parse_list_int(value);
//        }
//        else if (name == "max_trivial_size") {
//            max_trivial_size = parse_int(value);
//        }
//        else if (name == "min_wfa_size") {
//            min_wfa_size = parse_int(value);
//        }
//        else if (name == "max_wfa_size") {
//            max_wfa_size = parse_int(value);
//        }
//        else if (name == "max_wfa_ratio") {
//            max_wfa_ratio = parse_double(value);
//        }
//        else if (name == "wfa_pruning_dist") {
//            wfa_pruning_dist = parse_int(value);
//        }
//        else if (name == "deletion_alignment_ratio") {
//            deletion_alignment_ratio = parse_int(value);
//        }
//        else if (name == "deletion_alignment_short_max_size") {
//            deletion_alignment_short_max_size = parse_int(value);
//        }
//        else if (name == "deletion_alignment_long_min_size") {
//            deletion_alignment_long_min_size = parse_int(value);
//        }
//        else if (name == "cyclize_tandem_duplications") {
//            cyclize_tandem_duplications = parse_int(value);
//        }
//        else if (name == "max_tandem_duplication_search_rounds") {
//            max_tandem_duplication_search_rounds = parse_int(value);
//        }
//        else if (name == "min_cyclizing_length") {
//            min_cyclizing_length = parse_int(value);
//        }
//        else if (name == "preserve_subproblems") {
//            preserve_subproblems = parse_bool(value);
//        }
//        else if (name == "skip_calibration") {
//            skip_calibration = parse_bool(value);
//        }
//        else if (name == "subproblems_prefix") {
//            subproblems_prefix = value;
//        }
//        else if (name == "subalignments_filepath") {
//            subalignments_filepath = value;
//        }
//        else if (name == "tree_name") {
//            tree_name = value;
//        }
//        else if (name == "all_pairs_prefix") {
//            all_pairs_prefix = value;
//        }
//        else if (name == "fasta_name") {
//            fasta_name = value;
//        }
//        else if (name == "bonds_prefix") {
//            bonds_prefix = value;
//        }
//        else {
//            throw runtime_error("Config contains unrecognized variable: " + name);
//        }
//    }
//}

void NewParameters::initialize_submodule(submodule_t module, const std::string& description) {
    params[module].first = description;
}

 
string NewParameters::generate_config() const {
    
    stringstream strm;
    strm << setprecision(17); // the number of significant digits in a double
    strm << "---\n";
    for (const auto& param_set : params) {
        // the header for this group of parameters
        strm << ' ' << '\n';
        strm << ' ' << "##########\n";
        strm << ' ' << "#" << param_set.second.first << "\n";
        strm << ' ' << "##########\n";
        strm << ' ' << '\n';
        for (const auto& param : param_set.second.second) {
            for (auto help_line : tokenize(param.get_help(), '\n')) {
                strm << ' ' << "# " << help_line << '\n';
            }
            strm << ' ' << param.get_name() << ": " << param.value_str() << '\n';
        }
    }
    
    return strm.str();
}


//void NewParameters::validate() const {
//    if (max_count <= 0) {
//        throw runtime_error("Got non-positive value " + to_string(max_count) + " for maximum anchor occurrence count");
//    }
//    if (max_num_match_pairs <= 0) {
//        throw runtime_error("Got non-positive value " + to_string(max_num_match_pairs) + " for maximum number of anchor matches");
//    }
//    if (anchor_score_function < 0 || anchor_score_function > 3) {
//        throw runtime_error("Got invalid value " + to_string(anchor_score_function) + " for anchor scoring function");
//    }
//    if (pair_count_power < 0.0) {
//        throw runtime_error("Got negative value " + to_string(pair_count_power) + " for match count scoring power");
//    }
//    if (length_intercept < 1.0) {
//        throw runtime_error("Got value " + to_string(length_intercept) + " below 1 for length intercept / maximum positive-scoring length");
//    }
//    if (length_decay_power < 0.0) {
//        throw runtime_error("Got negative value " + to_string(length_decay_power) + " for length decay power");
//    }
//    if (anchor_gap_open.size() != 3 || anchor_gap_extend.size() != 3) {
//        throw runtime_error("Anchor gap open/extend penalties must have length 3, got " + list_to_string(anchor_gap_open) + " and " + list_to_string(anchor_gap_extend));
//    }
//    for (size_t i = 1; i < anchor_gap_open.size(); ++i) {
//        if (anchor_gap_open[i - 1] > anchor_gap_open[i] || anchor_gap_extend[i - 1] < anchor_gap_extend[i]) {
//            throw runtime_error("Anchor gap open penalities must be provided in increasing order and gap extend penalties must be provided in decreasing order, got " + list_to_string(anchor_gap_open) + " and " + list_to_string(anchor_gap_extend));
//        }
//    }
//    if (logging_level < 0 || logging_level > 4) {
//        throw runtime_error("Got invalid value " + to_string(logging_level) + " for logging level");
//    }
//    if (chaining_algorithm < 0 || chaining_algorithm > 2) {
//        throw runtime_error("Got invalid value " + to_string(chaining_algorithm) + " for chaining algorithm");
//    }
//    if (constraint_method < 0 || constraint_method > 3) {
//        throw runtime_error("Got invalid value " + to_string(constraint_method) + " for constraint method");
//    }
//    if (minimum_segment_score < 0.0) {
//        throw runtime_error("Got negative value " + to_string(minimum_segment_score) + " for minimum segment score");
//    }
//    if (minimum_segment_average < 0.0) {
//        throw runtime_error("Got negative value " + to_string(minimum_segment_average) + " for minimum segment windowed-average score");
//    }
//    if (window_length <= 0.0) {
//        throw runtime_error("Got non-positive value " + to_string(window_length) + " for windowed-average window length");
//    }
//    if (stitch_gap_open.size() != 3 || stitch_gap_extend.size() != 3) {
//        throw runtime_error("Stitching gap open/extend penalties must have length 3, got " + list_to_string(stitch_gap_open) + " and " + list_to_string(stitch_gap_extend));
//    }
//    for (size_t i = 1; i < stitch_gap_open.size(); ++i) {
//        if (stitch_gap_open[i - 1] > stitch_gap_open[i] || stitch_gap_extend[i - 1] < stitch_gap_extend[i]) {
//            throw runtime_error("Stitching gap open penalities must be provided in increasing order and gap extend penalties must be provided in decreasing order, got " + list_to_string(stitch_gap_open) + " and " + list_to_string(stitch_gap_extend));
//        }
//    }
//    int64_t stitch_max = numeric_limits<uint32_t>::max();
//    if (stitch_match < 0) {
//        throw runtime_error("Got negative stitching match bonus " + to_string(stitch_match));
//    }
//    if (stitch_match > stitch_max) {
//        throw runtime_error("Stitching match bonus must be less than " + to_string(stitch_max) + ", got " + to_string(stitch_match));
//    }
//    if (stitch_mismatch < 0) {
//        throw runtime_error("Got negative stitching mismatch penalty " + to_string(stitch_mismatch));
//    }
//    if (stitch_mismatch > stitch_max) {
//        throw runtime_error("Stitching mismatch penalty must be less than " + to_string(stitch_max) + ", got " + to_string(stitch_mismatch));
//    }
//    for (size_t i = 0; i < stitch_gap_open.size(); ++i) {
//        if (stitch_gap_open[i] < 0) {
//            throw runtime_error("Got negative stitching gap open penalty " + to_string(stitch_gap_open[i]));
//        }
//        if (stitch_gap_open[i] > stitch_max) {
//            throw runtime_error("Stitching gap open penalty must be less than " + to_string(stitch_max) + ", got " + to_string(stitch_gap_open[i]));
//        }
//        if (stitch_gap_extend[i] < 0) {
//            throw runtime_error("Got negative stitching gap extend penalty " + to_string(stitch_gap_extend[i]));
//        }
//        if (stitch_gap_extend[i] > stitch_max) {
//            throw runtime_error("Stitching gap extend penalty must be less than " + to_string(stitch_max) + ", got " + to_string(stitch_gap_extend[i]));
//        }
//    }
//    if (max_trivial_size < 0) {
//        throw runtime_error("Got negative maximum size for a trivial DP matrix " + to_string(max_trivial_size));
//    }
//    if (min_wfa_size < 0) {
//        throw runtime_error("Got negative mininmum size of DP matrix for performing WFA " + to_string(min_wfa_size));
//    }
//    if (max_wfa_size < 0) {
//        throw runtime_error("Got negative maximum size of DP matrix for performing WFA " + to_string(max_wfa_size));
//    }
//    if (max_wfa_size < min_wfa_size) {
//        throw runtime_error("Got greater minimum than maximum for size of DP matrix for performing WFA: " + to_string(min_wfa_size) + ", " + to_string(max_wfa_size));
//    }
//    if (max_wfa_ratio < 1.0) {
//        throw runtime_error("Got maximum ratio of longer:shorter graph size to perform WFA that is less than 1.0: " + to_string(max_wfa_ratio));
//    }
//    if (wfa_pruning_dist < 0) {
//        throw runtime_error("Got negative pruning distance for WFA " + to_string(wfa_pruning_dist));
//    }
//    if (deletion_alignment_ratio < 1) {
//        throw runtime_error("Got minimum ratio of graph sizes to be considered a deletion that is less than 1: " + to_string(deletion_alignment_ratio));
//    }
//    if (deletion_alignment_short_max_size < 0) {
//        throw runtime_error("Got negative maximum size for shorter graph to be considered a deletion: " + to_string(deletion_alignment_short_max_size));
//    }
//    if (deletion_alignment_long_min_size < 0) {
//        throw runtime_error("Got negative minimum size for longer graph to be considered a deletion: " + to_string(deletion_alignment_long_min_size));
//    }
//    if (cyclize_tandem_duplications && max_tandem_duplication_search_rounds == 0) {
//        throw runtime_error("Cannot cyclize tandem duplications with 0 rounds of tandem duplication search");
//    }
//    if (max_tandem_duplication_search_rounds < 0) {
//        throw runtime_error("Got number of tandem duplication search rounds that is negative: " + to_string(max_tandem_duplication_search_rounds));
//    }
//    if (min_cyclizing_length < 0) {
//        throw runtime_error("Got minimum cyclizing tandem duplication size that is negative: " + to_string(min_cyclizing_length));
//    }
//    if (fasta_name.empty()) {
//        throw runtime_error("FASTA input is missing");
//    }
//}


//void NewParameters::apply(Core& core) const {
//
//    // pass through parameters
//    core.path_match_finder.use_color_set_size = use_color_set_size;
//    core.path_match_finder.max_count = max_count;
//
//    core.score_function.anchor_score_function = anchor_score_function;
//    core.score_function.pair_count_power = pair_count_power;
//    core.score_function.length_intercept = length_intercept;
//    core.score_function.length_decay_power = length_decay_power;
//
//    core.anchorer.chaining_algorithm = chaining_algorithm;
//    core.anchorer.max_num_match_pairs = max_num_match_pairs;
//    for (size_t i = 0; i < core.anchorer.gap_open.size(); ++i) {
//        core.anchorer.gap_open[i] = anchor_gap_open[i];
//        core.anchorer.gap_extend[i] = anchor_gap_extend[i];
//    }
//    core.anchorer.do_fill_in_anchoring = do_fill_in_anchoring;
//
//    core.partitioner.constraint_method = constraint_method;
//    core.partitioner.minimum_segment_score = minimum_segment_score;
//    core.partitioner.minimum_segment_average = minimum_segment_average;
//    core.partitioner.window_length = window_length;
//    core.partitioner.generalized_length_mean = generalized_length_mean;
//
//    core.stitcher.alignment_params.match = stitch_match;
//    core.stitcher.alignment_params.mismatch = stitch_mismatch;
//    for (size_t i = 0; i < core.stitcher.alignment_params.gap_open.size(); ++i) {
//        core.stitcher.alignment_params.gap_open[i] = stitch_gap_open[i];
//        core.stitcher.alignment_params.gap_extend[i] = stitch_gap_extend[i];
//    }
//    core.stitcher.max_trivial_size = max_trivial_size;
//    core.stitcher.min_wfa_size = min_wfa_size;
//    core.stitcher.max_wfa_size = max_wfa_size;
//    core.stitcher.max_wfa_ratio = max_wfa_ratio;
//    core.stitcher.wfa_pruning_dist = wfa_pruning_dist;
//    core.stitcher.deletion_alignment_ratio = deletion_alignment_ratio;
//    core.stitcher.deletion_alignment_short_max_size = deletion_alignment_short_max_size;
//    core.stitcher.deletion_alignment_long_min_size = deletion_alignment_long_min_size;
//
//    core.bonder.min_length = min_cyclizing_length;
//
//    core.cyclize_tandem_duplications = cyclize_tandem_duplications;
//    core.max_tandem_duplication_search_rounds = max_tandem_duplication_search_rounds;
//    core.preserve_subproblems = preserve_subproblems;
//    core.skip_calibration = skip_calibration;
//
//    core.subproblems_prefix = subproblems_prefix;
//    core.subalignments_filepath = subalignments_filepath;
//    core.induced_pairwise_prefix = all_pairs_prefix;
//    core.bonds_prefix = bonds_prefix;
//}

//bool NewParameters::operator==(const NewParameters& other) const {
//
//    return (use_color_set_size == other.use_color_set_size &&
//            max_count == other.max_count &&
//            max_num_match_pairs == other.max_num_match_pairs &&
//            chaining_algorithm == other.chaining_algorithm &&
//            constraint_method == other.constraint_method &&
//            minimum_segment_score == other.minimum_segment_score &&
//            minimum_segment_average == other.minimum_segment_average &&
//            window_length == other.window_length &&
//            generalized_length_mean == other.generalized_length_mean &&
//            stitch_match == other.stitch_match &&
//            stitch_mismatch == other.stitch_mismatch &&
//            stitch_gap_open == other.stitch_gap_open &&
//            stitch_gap_extend == other.stitch_gap_extend &&
//            max_trivial_size == other.max_trivial_size &&
//            min_wfa_size == other.min_wfa_size &&
//            max_wfa_size == other.max_wfa_size &&
//            max_wfa_ratio == other.max_wfa_ratio &&
//            wfa_pruning_dist == other.wfa_pruning_dist &&
//            deletion_alignment_ratio == other.deletion_alignment_ratio &&
//            deletion_alignment_short_max_size == other.deletion_alignment_short_max_size &&
//            deletion_alignment_long_min_size == other.deletion_alignment_long_min_size &&
//            pair_count_power == other.pair_count_power &&
//            length_intercept == other.length_intercept &&
//            length_decay_power == other.length_decay_power &&
//            anchor_gap_open == other.anchor_gap_open &&
//            anchor_gap_extend == other.anchor_gap_extend &&
//            do_fill_in_anchoring == other.do_fill_in_anchoring &&
//            logging_level == other.logging_level &&
//            preserve_subproblems == other.preserve_subproblems &&
//            skip_calibration == other.skip_calibration &&
//            subproblems_prefix == other.subproblems_prefix &&
//            subalignments_filepath == other.subalignments_filepath &&
//            tree_name == other.tree_name &&
//            all_pairs_prefix == other.all_pairs_prefix &&
//            fasta_name == other.fasta_name &&
//            bonds_prefix == other.bonds_prefix);
//}
//
//bool NewParameters::operator!=(const NewParameters& other) const {
//    return !(*this == other);
//}

NewParameters::Parameter::Parameter(const Parameter& other) : type(other.type), name(other.name), help(other.help) {
    switch (type) {
        case Integer:
        case Enum:
            value.i = other.value.i;
            break;
        case Double:
            value.d = other.value.d;
            break;
        case Bool:
            value.b = other.value.b;
            break;
        case String:
            value.s = new std::string(*other.value.s);
            break;
        case DoubleArray3:
            value.da = other.value.da;
            break;
        case IntegerArray3:
            value.ia = other.value.ia;
            break;
        default:
            throw std::runtime_error("Unrecognized parameter data type " + std::to_string((int) type));
            break;
    }
}

NewParameters::Parameter::Parameter(Parameter&& other) : type(other.type), name(std::move(other.name)), help(std::move(other.help)) {
    switch (type) {
        case Integer:
        case Enum:
            value.i = other.value.i;
            break;
        case Double:
            value.d = other.value.d;
            break;
        case Bool:
            value.b = other.value.b;
            break;
        case String:
            value.s = other.value.s;
            other.value.s = nullptr;
            break;
        case DoubleArray3:
            value.da = other.value.da;
            break;
        case IntegerArray3:
            value.ia = other.value.ia;
            break;
        default:
            throw std::runtime_error("Unrecognized parameter data type " + std::to_string((int) type));
            break;
    }
}

NewParameters::Parameter::~Parameter() {
    
    static const std::map<type_t, std::string> type_names = {
        {Integer, "integer"},
        {Double, "double"},
        {Bool, "bool"},
        {Enum, "enum"},
        {String, "string"},
        {DoubleArray3, "double[3]"},
        {IntegerArray3, "integer[3]"}
    };
    if (type == String) {
        delete value.s;
    }
}

NewParameters::type_t NewParameters::Parameter::get_type() const {
    return type;
}

const std::string& NewParameters::Parameter::get_name() const {
    return name;
}

const std::string& NewParameters::Parameter::get_help() const {
    return help;
}

std::string NewParameters::Parameter::value_str() const {
    std::stringstream strm;
    switch (type) {
        case Integer:
        case Enum:
            strm << value.i;
            break;
        case Bool:
            strm << value.b;
            break;
        case Double:
            strm << value.d;
            break;
        case String:
            strm << *value.s;
            break;
        case DoubleArray3:
            for (size_t i = 0; i < value.da.size(); ++i) {
                if (i) {
                    strm << ',';
                }
                strm << value.da[i];
            }
            break;
        case IntegerArray3:
            for (size_t i = 0; i < value.ia.size(); ++i) {
                if (i) {
                    strm << ',';
                }
                strm << value.ia[i];
            }
            break;
            
        default:
            throw std::runtime_error("Unrecognized parameter value type " + std::to_string(type));
            break;
    }
    return strm.str();
}

}
