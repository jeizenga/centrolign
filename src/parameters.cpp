#include "centrolign/parameters.hpp"

#include <cassert>
#include <functional>
#include <iomanip>

#include "centrolign/utility.hpp"


namespace centrolign {

Parameters::Parameters() {
    
    initialize_submodule(IO, "Parameters related to file I/O and logging");
    initialize_submodule(MatchFinding, "Parameters related to identifying matches between graphs");
    initialize_submodule(Anchoring, "Parameters related to identifying high-scoring chains of matches to anchor alignments");
    initialize_submodule(IdentifyingAlignability, "Parameters related to determining whether a graph region is alignable");
    initialize_submodule(Aligning, "Parameters related to constructing a base-level alignment");
    initialize_submodule(InducingCycles, "Parameters related to inducing cycles at tandem duplications");
    initialize_submodule(DeveloperTools, "Parameters that were designed only to facilitate software development");
    
    add_parameter(IO, "fasta_name", String, std::string(), "The path to a FASTA file containing all of input sequences");
    add_parameter(IO, "tree_name", String, std::string(), "The path to a guide tree for the alignment in Newick format (sample names must match sequence names from the input FASTA)");
    add_parameter(IO, "logging_level", Enum, logging::Basic, "The level of verbosity of logging to stderr during execution:\n"
                  "- " + std::to_string((int) logging::Silent) + ": Silent\n"
                  "- " + std::to_string((int) logging::Minimal) + ": Minimal\n"
                  "- " + std::to_string((int) logging::Basic) + ": Basic\n"
                  "- " + std::to_string((int) logging::Verbose) + ": Verbose\n"
                  "- " + std::to_string((int) logging::Debug) + ": Debug");
    add_parameter(IO, "subproblems_prefix", String, std::string(), "If provided, save the results of the intermediate subproblems in GFA format to files with this prefix");
    add_parameter(IO, "restart", Bool, false, "Attempt to restart mid-execution using the saved partial results from 'subproblems_prefix'");
    add_parameter(IO, "all_pairs_prefix", String, std::string(), "If provided, save the induced pairwise alignment for each pair of sequences in CIGAR format to files with this prefix");
    add_parameter(IO, "subalignments_filepath", String, std::string(), "If provided, save the path-to-path alignment from each subproblem to files with this prefix");
    
    add_parameter(MatchFinding, "max_count", Integer, 3000, "Only query matches that occur at most this many times on either of the two graphs");
    add_parameter(MatchFinding, "use_color_set_size", Bool, true, "Use Hui's (1992) color set size index instead of a merge sort tree (CSS is generally faster and uses less memory)");
    
    add_parameter(Anchoring, "max_num_match_pairs", Integer, 1250000, "The maximum number of matches between two graphs that will be considered during chaining");
    add_parameter(Anchoring, "do_fill_in_anchoring", Bool, true, "Attempt to fill in gaps in the anchor chain using matches that were not considered due to the limit on the maximum number of matches");
    add_parameter(Anchoring, "global_anchoring", Bool, true, "Identify chains that cover the whole sequence, as opposed to local regions");
    add_parameter(Anchoring, "split_matches_at_branchpoints", Bool, true, "Allow the chaining algorithm to split anchors at forking paths in the graph to avoid reachability artifacts");
    add_parameter(Anchoring, "anchor_split_limit", Integer, 5, "If splitting at branch points, how close to the end of the anchor must the split be");
    add_parameter(Anchoring, "min_split_length", Integer, 128, "If splitting at branch points, only split anchors that are at least this long");
    add_parameter(Anchoring, "min_path_length_spread", Integer, 50, "If splitting at branch points, only split anchors at forks whose paths differ by at least this much in length");
    add_parameter(Anchoring, "max_split_match_set_size", Integer, 16, "If splitting at branch points, only split anchors with at most this many matching sequences");
    add_parameter(Anchoring, "chaining_algorithm", Enum, Anchorer::SparseAffine, "The chaining algorithm used:\n"
                  "- " + std::to_string((int) Anchorer::Exhaustive) + ": Simple exhaustive algorithm (slow)\n"
                  "- " + std::to_string((int) Anchorer::Sparse) + ": Sparse algorithm with no gap penalties\n"
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
    add_parameter(IdentifyingAlignability, "generalized_length_mean", Double, -0.5, "Parameter of the Holder mean used to combine the lengths on the two graphs into a single length measurement");
    add_parameter(IdentifyingAlignability, "boundary_score_factor", Double, 0.95, "When realigning regions after inducing cycles, treat the boundaries of the realignment as having score equal to this proportion times the minimum segment score");
    
    add_parameter(Aligning, "stitch_match", Integer, 20, "Match value when stitching anchors into a base-level alignment");
    add_parameter(Aligning, "stitch_mismatch", Integer, 80, "Mismatch penalty when stitching anchors into a base-level alignment");
    add_parameter(Aligning, "stitch_gap_open", IntegerArray3, std::array<int64_t, 3>{60, 800, 2500}, "Piecewise affine gap open penalties when stitching anchors into a base-level alignment");
    add_parameter(Aligning, "stitch_gap_extend", IntegerArray3, std::array<int64_t, 3>{30, 5, 1}, "Piecewise affine gap extend penalties when stitching anchors into a base-level alignment");
    add_parameter(Aligning, "max_trivial_size", Integer, 30000, "Maximum size of a dynamic programming matrix that will be aligned even if it was identified as unalignable");
    add_parameter(Aligning, "min_wfa_size", Integer, 40000000, "Minimum size of a dynamic programming matrix that will be aligned using graph-graph WFA");
    add_parameter(Aligning, "max_wfa_size", Integer, 75000000, "Maximum size of a dynamic programming matrix that will be aligned using graph-graph WFA");
    add_parameter(Aligning, "max_wfa_ratio", Double, 1.05, "Maximum ratio of long-to-short side of the dynamic programming matrix for graph-graph WFA to be used");
    add_parameter(Aligning, "wfa_pruning_dist", Integer, 25, "The lagging distance for a diagonal to be pruned in graph-graph WFA");
    add_parameter(Aligning, "deletion_alignment_ratio", Integer, 8, "The minimum ratio of long-to-short side of the dynamic programming matrix to use WFA-based implied deletion algorithm");
    add_parameter(Aligning, "deletion_alignment_short_max_size", Integer, 1500, "The maximum size of the short side of the dynamic programming matrix to use WFA-based implied deletion algorithm");
    add_parameter(Aligning, "deletion_alignment_long_min_size", Integer, 2000, "The minimum size of the long side of the dynamic programming matrix to use WFA-based implied deletion algorithm");
    add_parameter(Aligning, "indel_fuzz_score_proportion", Double, 0.001, "Remove low-scoring anchors that are restricting the location of large indels when their score is worth at most this proportion of their neighboring anchors");
    add_parameter(Aligning, "min_indel_fuzz_length", Integer, 50, "When removing low-scoring anchors to de-specify the location of a indel, require the indel to be at least this long");
    
    add_parameter(InducingCycles, "cyclize_tandem_duplications", Bool, false, "Identify tandem duplications in the sequences and use them to induce cycles in the final graph");
    add_parameter(InducingCycles, "max_tandem_duplication_search_rounds", Integer, 3, "The maximum number of nested tandem duplications to attempt finding for any given subsequence");
    add_parameter(InducingCycles, "min_cyclizing_length", Integer, 100000, "The minimum size of a tandem duplication to look for");
    add_parameter(InducingCycles, "tandem_dup_score_proportion", Double, 0.2, "Require tandem duplication anchor chains to have at least this proportion of the score of the corresponding section of a self-to-self anchor chain");
    add_parameter(InducingCycles, "include_tandem_dup_gap_scores", Bool, true, "When computing the score of tandem duplication chains, include the gap scores");
    add_parameter(InducingCycles, "deviation_drift_factor", Double, 150.0, "When identifying tandem duplications, allow the chain to have indel deviations of this much times sqrt(length)");
    add_parameter(InducingCycles, "separation_drift_factor", Double, 50.0, "When identifying tandem duplications, require the chain to be separated from the main diagonal by the length minus this much times sqrt(length)");
    add_parameter(InducingCycles, "trim_window_proportion", Double, 0.1, "Trim off the ends of tandem duplications until they meet the minimum score requirement using only a window on each end of length equal to this proportion times 'min_cyclizing_length'");
    add_parameter(InducingCycles, "deduplication_slosh_proportion", Double, 0.1, "Consider two tandem duplications to be the same if they differ by at most this much times 'min_cyclizing_length'");
    add_parameter(InducingCycles, "max_realignment_cycle_size", Integer, 10000, "After cyclizing, realign cycles shorter than this length");
    add_parameter(InducingCycles, "inconsistent_indel_window", Integer, 100, "After cyclizing, look for inconsistently-placed indels to realign that are separated by at most this length");
    add_parameter(InducingCycles, "min_inconsistency_disjoint_length", Integer, 8, "Require inconsistently-placed indels to have disjoint un-merged sequences of at least this length from two segments of the same input sequence");
    add_parameter(InducingCycles, "min_inconsistency_total_length", Integer, 50, "Require inconsistently-placed indels to have total un-merged sequences of at least this length from two segments of the same input sequence");
    add_parameter(InducingCycles, "realignment_min_padding", Integer, 1000, "When realigning after cyclizing, try to pad alignment problems with this much padding sequence from every path");
    add_parameter(InducingCycles, "realignment_max_padding", Integer, 10000, "When realigning after cyclizing, stop adding padding if it would require any path to add this much sequence");
    
    
    add_parameter(DeveloperTools, "bonds_prefix", String, std::string(), "If provided, save the alignments of all tandem duplications identified in the cyclization process to files with this prefix");
    add_parameter(DeveloperTools, "preserve_subproblems", Bool, false, "Do not clear out data from completed subproblems as the algorithm progresses");
    add_parameter(DeveloperTools, "skip_calibration", Bool, false, "Do not calibrate the scoring parameters to the input sequences' repetitiveness");
    
}

void Parameters::apply(Core& core) const {
    
    // pass through parameters
    
    // note: fasta and tree handled in main function
    core.subproblems_prefix = parameter("subproblems_prefix").get<std::string>();
    core.subalignments_filepath = parameter("subalignments_filepath").get<std::string>();
    core.induced_pairwise_prefix = parameter("all_pairs_prefix").get<std::string>();
    core.bonds_prefix = parameter("bonds_prefix").get<std::string>();
    
    core.path_match_finder.use_color_set_size = parameter("use_color_set_size").get<bool>();
    core.path_match_finder.max_count = parameter("max_count").get<int64_t>();
    
    core.score_function.anchor_score_function = parameter("anchor_score_function").get<ScoreFunction::AnchorScore>();
    core.score_function.pair_count_power = parameter("pair_count_power").get<double>();
    core.score_function.length_intercept = parameter("length_intercept").get<double>();
    core.score_function.length_decay_power = parameter("length_decay_power").get<double>();
    
    core.anchorer.chaining_algorithm = parameter("chaining_algorithm").get<Anchorer::ChainAlgorithm>();
    core.anchorer.do_fill_in_anchoring = parameter("do_fill_in_anchoring").get<bool>();
    core.anchorer.max_num_match_pairs = parameter("max_num_match_pairs").get<int64_t>();
    core.anchorer.global_anchoring = parameter("global_anchoring").get<bool>();
    core.anchorer.split_matches_at_branchpoints = parameter("split_matches_at_branchpoints").get<bool>();
    core.anchorer.anchor_split_limit = parameter("anchor_split_limit").get<int64_t>();
    core.anchorer.min_split_length = parameter("min_split_length").get<int64_t>();
    core.anchorer.min_path_length_spread = parameter("min_path_length_spread").get<int64_t>();
    core.anchorer.max_split_match_set_size = parameter("max_split_match_set_size").get<int64_t>();
    auto anchor_gap_open = parameter("anchor_gap_open").get<std::array<double, 3>>();
    auto anchor_gap_extend = parameter("anchor_gap_extend").get<std::array<double, 3>>();
    for (size_t i = 0; i < core.anchorer.gap_open.size(); ++i) {
        core.anchorer.gap_open[i] = anchor_gap_open[i];
        core.anchorer.gap_extend[i] = anchor_gap_extend[i];
    }
    
    core.partitioner.constraint_method = parameter("constraint_method").get<Partitioner::ConstraintMethod>();
    core.partitioner.minimum_segment_score = parameter("minimum_segment_score").get<double>();
    core.partitioner.minimum_segment_average = parameter("minimum_segment_average").get<double>();
    core.partitioner.window_length = parameter("window_length").get<double>();
    core.partitioner.generalized_length_mean = parameter("generalized_length_mean").get<double>();
    core.partitioner.boundary_score_factor = parameter("boundary_score_factor").get<double>();
    
    core.stitcher.alignment_params.match = parameter("stitch_match").get<int64_t>();
    core.stitcher.alignment_params.mismatch = parameter("stitch_mismatch").get<int64_t>();
    auto stitch_gap_open = parameter("stitch_gap_open").get<std::array<int64_t, 3>>();
    auto stitch_gap_extend = parameter("stitch_gap_extend").get<std::array<int64_t, 3>>();
    for (size_t i = 0; i < core.stitcher.alignment_params.gap_open.size(); ++i) {
        core.stitcher.alignment_params.gap_open[i] = stitch_gap_open[i];
        core.stitcher.alignment_params.gap_extend[i] = stitch_gap_extend[i];
    }
    core.stitcher.max_trivial_size = parameter("max_trivial_size").get<int64_t>();
    core.stitcher.min_wfa_size = parameter("min_wfa_size").get<int64_t>();
    core.stitcher.max_wfa_size = parameter("max_wfa_size").get<int64_t>();
    core.stitcher.max_wfa_ratio = parameter("max_wfa_ratio").get<double>();
    core.stitcher.wfa_pruning_dist = parameter("wfa_pruning_dist").get<int64_t>();
    core.stitcher.deletion_alignment_ratio = parameter("deletion_alignment_ratio").get<int64_t>();
    core.stitcher.deletion_alignment_short_max_size = parameter("deletion_alignment_short_max_size").get<int64_t>();
    core.stitcher.deletion_alignment_long_min_size = parameter("deletion_alignment_long_min_size").get<int64_t>();
    core.stitcher.indel_fuzz_score_proportion = parameter("indel_fuzz_score_proportion").get<double>();
    core.stitcher.min_indel_fuzz_length = parameter("min_indel_fuzz_length").get<int64_t>();
    
    core.bonder.min_length = parameter("min_cyclizing_length").get<int64_t>();
    core.bonder.min_opt_proportion = parameter("tandem_dup_score_proportion").get<double>();
    core.bonder.include_gap_scores = parameter("include_tandem_dup_gap_scores").get<bool>();
    core.bonder.deviation_drift_factor = parameter("deviation_drift_factor").get<double>();
    core.bonder.separation_drift_factor = parameter("separation_drift_factor").get<double>();
    core.bonder.deduplication_slosh_proportion = parameter("deduplication_slosh_proportion").get<double>();
    core.bonder.trim_window_proportion = parameter("trim_window_proportion").get<double>();
    
    core.inconsistency_identifier.max_tight_cycle_size = parameter("max_realignment_cycle_size").get<int64_t>();
    core.inconsistency_identifier.max_bond_inconsistency_window = parameter("inconsistent_indel_window").get<int64_t>();
    core.inconsistency_identifier.min_inconsistency_disjoint_length = parameter("min_inconsistency_disjoint_length").get<int64_t>();
    core.inconsistency_identifier.min_inconsistency_total_length = parameter("min_inconsistency_total_length").get<int64_t>();
    core.inconsistency_identifier.padding_target_min_length = parameter("realignment_min_padding").get<int64_t>();
    core.inconsistency_identifier.padding_max_length_limit = parameter("realignment_max_padding").get<int64_t>();
    
    core.cyclize_tandem_duplications = parameter("cyclize_tandem_duplications").get<bool>();
    core.max_tandem_duplication_search_rounds = parameter("max_tandem_duplication_search_rounds").get<int64_t>();
    
    core.preserve_subproblems = parameter("preserve_subproblems").get<bool>();
    core.skip_calibration = parameter("skip_calibration").get<bool>();
}

// TODO: repetitive
template<typename T, size_t N>
std::array<T, N> parse_array(const std::string& value, const std::function<T(const std::string&)>& parser) {
    std::array<T, N> values;
    auto tokens = tokenize(value, ',');
    if (tokens.size() != values.size()) {
        throw std::runtime_error("Expected " + std::to_string(N) + " comma-separated values in input string '" + value + "'");
    }
    for (size_t i = 0; i < values.size(); ++i) {
        values[i] = parser(tokens[i]);
    }
    return values;
}

template<size_t N>
inline std::array<double, N> parse_array_double(const std::string& value) {
    return parse_array<double, N>(value, [](const std::string& str) { return parse_double(str); });
}
template<size_t N>
inline std::array<int64_t, N> parse_array_int(const std::string& value) {
    return parse_array<int64_t, N>(value, [](const std::string& str) { return parse_int(str); });
}

Parameters::Parameters(std::istream& in) : Parameters() {

    std::string config;
    {
        std::stringstream strm;
        strm << in.rdbuf();
        config = std::move(strm.str());
    }

    auto yaml_delim = config.find("---");
    if (yaml_delim < config.size()) {
        // skip ahead of the new file delimiter
        config = config.substr(yaml_delim + 3, config.size());
    }


    std::stringstream scan_strm(config);
    std::string line;
    while (scan_strm) {
        getline(scan_strm, line);
        
        size_t comment_start = line.find('#');
        auto end = comment_start < line.size() ? line.begin() + comment_start : line.end();
        if (all_of(line.begin(), end, [](char c) { return isspace(c); })) {
            // all whitespace or comment in this line
            continue;
        }

        auto colon_pos = find(line.begin(), end, ':');
        if (colon_pos >= end) {
            throw std::runtime_error("Config has line missing ':' delimiter in line '" + line + "'");
        }
        int name_start = 0;
        int name_end = (colon_pos - line.begin());
        int value_start = (colon_pos - line.begin()) + 1;
        int value_end = end - line.begin();

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
            throw std::runtime_error("Config has line with no variable name in line '" + line + "'");
        }
        if (value_start >= value_end) {
            throw std::runtime_error("Config has line with no variable value in line '" + line + "'");
        }
        if (line[value_start] == '"' &&
            (value_end == value_start + 1 || line[value_end - 1] != '"')) {
            throw std::runtime_error("Unmatched \" in config file line: " + line);
        }

        if (line[value_start] == '"') {
            ++value_start;
            --value_end;
        }

        std::string name = line.substr(name_start, name_end - name_start);
        std::string value = line.substr(value_start, value_end - value_start);
        
        auto& param = parameter(name);
        switch (param.get_type()) {
            case Integer:
            case Enum:
                param.set<int64_t>(parse_int(value));
                break;
            case Bool:
                param.set<bool>(parse_bool(value));
                break;
            case Double:
                param.set<double>(parse_double(value));
                break;
            case String:
                param.set<std::string>(value);
                break;
            case DoubleArray3:
                param.set<std::array<double, 3>>(parse_array_double<3>(value));
                break;
            case IntegerArray3:
                param.set<std::array<int64_t, 3>>(parse_array_int<3>(value));
                break;
            default:
                throw std::runtime_error("Unrecognized parameter type");
                break;
        }
    }
}

void Parameters::initialize_submodule(submodule_t module, const std::string& description) {
    params[module].first = description;
}

const Parameters::Parameter& Parameters::parameter(const std::string& name) const {
    auto it = param_position.find(name);
    if (it == param_position.end()) {
        throw std::runtime_error("No parameter with name " + name);
    }
    return params.at(it->second.first).second.at(it->second.second);
}

Parameters::Parameter& Parameters::parameter(const std::string& name) {
    auto it = param_position.find(name);
    if (it == param_position.end()) {
        throw std::runtime_error("No parameter with name " + name);
    }
    return params.at(it->second.first).second.at(it->second.second);
}
 
std::string Parameters::generate_config() const {
    
    std::stringstream strm;
    strm << std::setprecision(17); // the number of significant digits in a double
    strm << "---\n";
    for (const auto& param_set : params) {
        // the header for this group of parameters
        strm << ' ' << '\n';
        strm << ' ' << "##########\n";
        strm << ' ' << "# " << param_set.second.first << "\n";
        strm << ' ' << "##########\n";
        strm << ' ' << '\n';
        for (const auto& param : param_set.second.second) {
            for (auto help_line : tokenize(param.get_help(), '\n')) {
                strm << ' ' << "# " << help_line << '\n';
            }
            auto value = param.value_str();
            bool has_whitespace = std::find_if(value.begin(), value.end(), [](char c) { return isspace(c); }) != value.end();
            bool needs_quotes = has_whitespace || value.empty();
            strm << ' ' << param.get_name() << ": " << (needs_quotes ? "\"" : "") << value << (needs_quotes ? "\"" : "") << '\n';
        }
    }
    
    return strm.str();
}


void Parameters::validate() const {
    
    enforce_gt<int64_t>("max_count", 0);
    enforce_gt<int64_t>("max_num_match_pairs", 0);
    enforce_range<int64_t>("anchor_score_function", 0,
                           container_max(std::vector<int>{ScoreFunction::InverseCount, ScoreFunction::LengthScaleInverseCount, ScoreFunction::ConcaveLengthScaleInverseCount, ScoreFunction::ConcaveLengthScaleCountDifference}));
    enforce_gt<double>("pair_count_power", 0.0);
    enforce_geq<double>("length_intercept", 1.0);
    enforce_geq<double>("length_decay_power", 0.0);
    enforce_geq<int64_t>("anchor_split_limit", 0);
    enforce_geq<int64_t>("min_path_length_spread", 0);
    enforce_geq<int64_t>("min_split_length", 0);
    enforce_geq<int64_t>("max_split_match_set_size", 0);
    enforce_range<int64_t>("logging_level", 0,
                           container_max(std::vector<int>{logging::Silent, logging::Minimal, logging::Basic, logging::Verbose, logging::Debug}));
    enforce_range<int64_t>("chaining_algorithm", 0,
                           container_max(std::vector<int>{Anchorer::Exhaustive, Anchorer::Sparse, Anchorer::SparseAffine}));
    enforce_range<int64_t>("constraint_method", 0,
                           container_max(std::vector<int>{Partitioner::Null, Partitioner::Unconstrained, Partitioner::MinAverage, Partitioner::MinWindowAverage}));
    enforce_geq<double>("minimum_segment_score", 0.0);
    enforce_geq<double>("minimum_segment_average", 0.0);
    enforce_gt<double>("window_length", 0.0);
    enforce_geq<double>("boundary_score_factor", 0.0);
    enforce_geq<int64_t>("stitch_match", 0);
    enforce_leq<int64_t>("stitch_match", std::numeric_limits<int32_t>::max());
    enforce_geq<int64_t>("stitch_mismatch", 0);
    enforce_leq<int64_t>("stitch_mismatch", std::numeric_limits<int32_t>::max());
    enforce_geq<int64_t>("max_trivial_size", 0);
    enforce_geq<int64_t>("min_wfa_size", 0);
    enforce_geq<int64_t>("max_wfa_size", 0);
    enforce_geq<double>("max_wfa_ratio", 1.0);
    enforce_geq<int64_t>("wfa_pruning_dist", 0);
    enforce_geq<int64_t>("deletion_alignment_ratio", 1);
    enforce_geq<int64_t>("deletion_alignment_short_max_size", 0);
    enforce_geq<int64_t>("deletion_alignment_long_min_size", 0);
    enforce_range<double>("indel_fuzz_score_proportion", 0.0, 1.0);
    enforce_geq<int64_t>("min_indel_fuzz_length", 0);
    if (parameter("cyclize_tandem_duplications").get<bool>()) {
        enforce_gt<int64_t>("max_tandem_duplication_search_rounds", 0);
    }
    enforce_geq<int64_t>("min_cyclizing_length", 0);
    enforce_range<double>("tandem_dup_score_proportion", 0.0, 1.0);
    enforce_geq<double>("deviation_drift_factor", 0.0);
    enforce_geq<double>("separation_drift_factor", 0.0);
    enforce_geq<double>("deduplication_slosh_proportion", 0.0);
    enforce_geq<double>("trim_window_proportion", 0.0);
    enforce_geq<int64_t>("max_realignment_cycle_size", 0);
    enforce_geq<int64_t>("inconsistent_indel_window", 0);
    enforce_geq<int64_t>("min_inconsistency_disjoint_length", 0);
    enforce_geq<int64_t>("min_inconsistency_total_length", 0);
    enforce_geq<int64_t>("realignment_min_padding", 0);
    enforce_geq<int64_t>("realignment_max_padding", 0);
    
    if (parameter("restart").get<bool>() && parameter("subproblems_prefix").get<std::string>().empty()) {
        throw std::runtime_error("Cannot restart mid-execution without setting 'subproblems_prefix'");
    }
    if (parameter("fasta_name").get<std::string>().empty()) {
        throw std::runtime_error("FASTA input is missing");
    }
    
    if (parameter("max_wfa_size").get<int64_t>() < parameter("min_wfa_size").get<int64_t>()) {
        throw std::runtime_error("Got greater minimum than maximum for size of DP matrix for performing WFA: " + std::to_string(parameter("max_wfa_size").get<int64_t>()) + ", " + std::to_string(parameter("min_wfa_size").get<int64_t>()));
    }
    auto anchor_gap_open = parameter("anchor_gap_open").get<std::array<double, 3>>();
    auto anchor_gap_extend = parameter("anchor_gap_extend").get<std::array<double, 3>>();
    for (size_t i = 1; i < anchor_gap_open.size(); ++i) {
        if (anchor_gap_open[i - 1] > anchor_gap_open[i] || anchor_gap_extend[i - 1] < anchor_gap_extend[i]) {
            throw std::runtime_error("Anchor gap open penalities must be provided in increasing order and gap extend penalties must be provided in decreasing order, got " + parameter("anchor_gap_open").value_str()  + " and " + parameter("anchor_gap_extend").value_str());
        }
    }
    for (size_t i = 0; i < anchor_gap_open.size(); ++i) {
        if (anchor_gap_open[i] < 0) {
            throw std::runtime_error("Got negative anchoring gap open penalty " + std::to_string(anchor_gap_open[i]));
        }
        if (anchor_gap_extend[i] < 0) {
            throw std::runtime_error("Got negative anchoring gap extend penalty " + std::to_string(anchor_gap_extend[i]));
        }
    }
    auto stitch_gap_open = parameter("stitch_gap_open").get<std::array<int64_t, 3>>();
    auto stitch_gap_extend = parameter("stitch_gap_extend").get<std::array<int64_t, 3>>();
    for (size_t i = 1; i < stitch_gap_open.size(); ++i) {
        if (stitch_gap_open[i - 1] > stitch_gap_open[i] || stitch_gap_extend[i - 1] < stitch_gap_extend[i]) {
            throw std::runtime_error("Stitching gap open penalities must be provided in increasing order and gap extend penalties must be provided in decreasing order, got " + parameter("stitch_gap_open").value_str() + " and " + parameter("stitch_gap_extend").value_str());
        }
    }
    for (size_t i = 0; i < stitch_gap_open.size(); ++i) {
        if (stitch_gap_open[i] < 0) {
            throw std::runtime_error("Got negative stitching gap open penalty " + std::to_string(stitch_gap_open[i]));
        }
        if (stitch_gap_open[i] > std::numeric_limits<int32_t>::max()) {
            throw std::runtime_error("Stitching gap open penalty must be less than " + std::to_string(std::numeric_limits<int32_t>::max()) + ", got " + std::to_string(stitch_gap_open[i]));
        }
        if (stitch_gap_extend[i] < 0) {
            throw std::runtime_error("Got negative stitching gap extend penalty " + std::to_string(stitch_gap_extend[i]));
        }
        if (stitch_gap_extend[i] > std::numeric_limits<int32_t>::max()) {
            throw std::runtime_error("Stitching gap extend penalty must be less than " + std::to_string(std::numeric_limits<int32_t>::max()) + ", got " + std::to_string(stitch_gap_extend[i]));
        }
    }
}


bool Parameters::operator==(const Parameters& other) const {
    // note: we only check for equivalent parameter sets, types, and values, not help documentation or ordering
    if (params.size() != other.params.size()) {
        return false;
    }
    
    for (const auto& param_set : params) {
        if (!other.params.count(param_set.first)) {
            return false;
        }
        
        for (const auto& param : param_set.second.second) {
            auto it = other.param_position.find(param.get_name());
            if (it == other.param_position.end()) {
                return false;
            }
            if (it->second.first != param_set.first) {
                return false;
            }
            const auto& other_param = other.params.at(it->second.first).second.at(it->second.second);
            assert(param.get_name() == other_param.get_name());
            if (param.get_type() != other_param.get_type()) {
                return false;
            }
            if (param.value_str() != other_param.value_str()) {
                return false;
            }
        }
    }
    return true;
}

bool Parameters::operator!=(const Parameters& other) const {
    return !(*this == other);
}

Parameters::Parameter::Parameter(const Parameter& other) : type(other.type), name(other.name), help(other.help) {
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
            if (other.value.s == nullptr) {
                value.s = nullptr;
            }
            else {
                value.s = new std::string(*other.value.s);
            }
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

Parameters::Parameter::Parameter(Parameter&& other) : type(other.type), name(std::move(other.name)), help(std::move(other.help)) {
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


Parameters::Parameter& Parameters::Parameter::operator=(Parameter&& other) {
    if (other.name != name || other.type != type) {
        throw std::runtime_error("Cannot assign Parameter of mismatching type or name");
    }
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
    return *this;
}

Parameters::Parameter& Parameters::Parameter::operator=(const Parameter& other) {
    if (other.name != name || other.type != type) {
        throw std::runtime_error("Cannot assign Parameter of mismatching type or name");
    }
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
            if (other.value.s == nullptr) {
                value.s = nullptr;
            }
            else {
                value.s = new std::string(*other.value.s);
            }
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
    return *this;
}

Parameters::Parameter::~Parameter() {
    if (type == String) {
        delete value.s;
    }
}

Parameters::type_t Parameters::Parameter::get_type() const {
    return type;
}

const std::string& Parameters::Parameter::get_name() const {
    return name;
}

const std::string& Parameters::Parameter::get_help() const {
    return help;
}

std::string Parameters::Parameter::value_str() const {
    std::stringstream strm;
    switch (type) {
        case Integer:
        case Enum:
            strm << value.i;
            break;
        case Bool:
            strm << (value.b ? "true" : "false");
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
