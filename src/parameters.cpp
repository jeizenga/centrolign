#include "centrolign/parameters.hpp"

#include "centrolign/utility.hpp"

#include <stdexcept>
#include <iostream>
#include <sstream>
#include <iomanip>

namespace centrolign {

using namespace std;

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
        else if (name == "max_count") {
            max_count = parse_int(value);
        }
        else if (name == "max_num_match_pairs") {
            max_num_match_pairs = parse_int(value);
        }
        else if (name == "pair_count_power") {
            pair_count_power = parse_double(value);
        }
        else if (name == "length_scale") {
            length_scale = parse_bool(value);
        }
        else if (name == "count_penalty_threshold") {
            count_penalty_threshold = parse_double(value);
        }
        else if (name == "logging_level") {
            logging_level = (logging::LoggingLevel) parse_int(value);
        }
        else if (name == "chaining_algorithm") {
            chaining_algorithm = (Anchorer::ChainAlgorithm) parse_int(value);
        }
        else if (name == "preserve_subproblems") {
            preserve_subproblems = parse_bool(value);
        }
        else if (name == "subproblems_prefix") {
            subproblems_prefix = value;
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
    strm << setprecision(18); // the number of significant digits in a double
    strm << "---\n";
    strm << ' ' << "simplify_window" << ": " << simplify_window << '\n';
    strm << ' ' << "max_walk_count" << ": " << max_walk_count << '\n';
    strm << ' ' << "blocking_allele_size" << ": " << blocking_allele_size << '\n';
    strm << ' ' << "path_matches" << ": " << path_matches << '\n';
    strm << ' ' << "max_count" << ": " << max_count << '\n';
    strm << ' ' << "max_num_match_pairs" << ": " << max_num_match_pairs << '\n';
    strm << ' ' << "pair_count_power" << ": " << pair_count_power << '\n';
    strm << ' ' << "length_scale" << ": " << length_scale << '\n';
    strm << ' ' << "count_penalty_threshold" << ": " << count_penalty_threshold << '\n';
    strm << ' ' << "logging_level" << ": " << (int) logging_level << '\n';
    strm << ' ' << "chaining_algorithm" << ": " << (int) chaining_algorithm << '\n';
    strm << ' ' << "preserve_subproblems" << ": " << preserve_subproblems << '\n';
    strm << ' ' << "subproblems_prefix" << ": " << string_or_null(subproblems_prefix) << '\n';
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
    if (pair_count_power < 0.0) {
        throw runtime_error("Got negative value " + to_string(pair_count_power) + " for match count scoring power");
    }
    if (count_penalty_threshold < 0.0) {
        throw runtime_error("Got negative value " + to_string(count_penalty_threshold) + " for match count scoring penalty threshold");
    }
    if (logging_level < 0 || logging_level > 4) {
        throw runtime_error("Got invalid value " + to_string(logging_level) + " for logging level");
    }
    if (chaining_algorithm < 0 || chaining_algorithm > 2) {
        throw runtime_error("Got invalid value " + to_string(chaining_algorithm) + " for chaining algorithm");
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
    core.match_finder.max_count = max_count;
    
    core.anchorer.max_num_match_pairs = max_num_match_pairs;
    core.anchorer.pair_count_power = pair_count_power;
    core.anchorer.length_scale = length_scale;
    core.anchorer.chaining_algorithm = chaining_algorithm;
    core.anchorer.count_penalty_threshold = count_penalty_threshold;
    
    core.preserve_subproblems = preserve_subproblems;
    
    core.subproblems_prefix = subproblems_prefix;
}

bool Parameters::operator==(const Parameters& other) const {
    
    return (simplify_window == other.simplify_window &&
            max_walk_count == other.max_walk_count &&
            blocking_allele_size == other.blocking_allele_size &&
            path_matches == other.path_matches &&
            max_count == other.max_count &&
            max_num_match_pairs == other.max_num_match_pairs &&
            pair_count_power == other.pair_count_power &&
            length_scale == other.length_scale &&
            count_penalty_threshold == other.count_penalty_threshold &&
            logging_level == other.logging_level &&
            chaining_algorithm == other.chaining_algorithm &&
            preserve_subproblems == other.preserve_subproblems &&
            subproblems_prefix == other.subproblems_prefix &&
            tree_name == other.tree_name &&
            all_pairs_prefix == other.all_pairs_prefix &&
            fasta_name == other.fasta_name);
}

bool Parameters::operator!=(const Parameters& other) const {
    return !(*this == other);
}

}
