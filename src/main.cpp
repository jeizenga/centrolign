#include <getopt.h>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <cstdint>
#include <limits>
#include <algorithm>

#include "centrolign/utility.hpp"
#include "centrolign/core.hpp"
#include "centrolign/logging.hpp"
#include "centrolign/gfa.hpp"
#include "centrolign/parameters.hpp"
#include "centrolign/version.hpp"

using namespace std;
using namespace centrolign;

void print_help() {
    
    Parameters defaults;
    
    cerr << "\n";
    cerr << "Usage: centrolign [options] sequences.fasta\n";
    cerr << "\n";
    cerr << "Options:\n";
    cerr << " --tree / -T FILE            Newick format guide tree for alignment [in FASTA order]\n";
    cerr << " --all-pairs / -A PREFIX     Output all induced pairwise alignments as files starting with PREFIX\n";
    cerr << " --all-subprobs / -S PREFIX  Output all subtree multiple sequence alignments as files starting with PREFIX\n";
    cerr << " --subalignments / -s FILE   Output a file containing the aligned path for each subproblem to FILE\n";
    cerr << " --cyclize / -c              Merge long tandem duplications into cycles in the final graph\n";
    cerr << " --cyclizing-size / -y INT   If cyclizing, merge tandem duplications larger than INT bases [" << defaults.min_cyclizing_length << "]\n";
    //cerr << " --simplify-window / -w      Size window used for graph simplification [" << defaults.simplify_window << "]\n";
    //cerr << " --simplify-count / -c       Number of walks through window that trigger simplification [" << defaults.max_walk_count << "]\n";
    //cerr << " --blocking-size / -b        Minimum size allele to block simplification [" << defaults.blocking_allele_size << "]\n";
    //cerr << " --non-path-matches / -P     Query matches on all walks through graph instead of only input sequences\n";
    cerr << " --max-count / -m INT        The maximum number of times an anchor can occur [" << defaults.max_count << "]\n";
    cerr << " --max-anchors / -a INT      The maximum number of anchors [" << defaults.max_num_match_pairs << "]\n";
    cerr << " --count-power / -p FLOAT    Scale anchor weights by the count raised to this power [" << defaults.pair_count_power << "]\n";
    //cerr << " --chain-alg / -g INT        Select from: 0 (exhaustive), 1 (sparse), 2 (sparse affine) [" << (int) defaults.chaining_algorithm << "]\n";
    cerr << " --no-unaln / -u             Do not attempt to identify unalignable regions\n";
    cerr << " --verbosity / -v INT        Select from: 0 (silent), 1 (minimal), 2 (basic), 3 (verbose), 4 (debug) [" << (int) defaults.logging_level << "]\n";
    cerr << " --config / -C FILE          Config file of parameters (overrides all other command line input)\n";
    cerr << " --restart / -R              Restart from a previous incomplete run (requires -S in first run)\n";
    cerr << " --help / -h                 Print this message and exit\n";
}

// make a dummy Newick string for in-order alignment
string in_order_newick_string(const vector<pair<string, string>>& sequences) {
    
    for (const auto& seq : sequences) {
        if (find(seq.first.begin(), seq.first.end(), '"') != seq.first.end()) {
            throw runtime_error("Sequence names cannot have internal quotation marks: " + seq.first);
        }
    }
    
    stringstream strm;
    for (size_t i = 1; i < sequences.size(); ++i) {
        strm << '(';
    }
    strm << '"' << sequences.front().first << '"';
    for (size_t i = 1; i < sequences.size(); ++i) {
        strm << ',' << '"' << sequences[i].first << '"' << ')';
    }
    strm << ';';
    return strm.str();
}

int main(int argc, char** argv) {
        
    // init the local params with the defaults
    Parameters params;
    
    // params that live outside the parameter object
    std::string config_name;
    bool restart = false;
    bool force_gfa_output = false;
    
    // opts without a short option
    static const int opt_skip_calibration = 1000;
    static const int opt_force_gfa_output = 1001;
    static const int opt_bond_prefix = 1003;
    while (true)
    {
        static struct option options[] = {
            {"tree", required_argument, NULL, 'T'},
            {"all-pairs", required_argument, NULL, 'A'},
            {"all-subprobs", required_argument, NULL, 'S'},
            {"subalignments", required_argument, NULL, 's'},
            {"cyclize", no_argument, NULL, 'c'},
            {"cyclizing-size", required_argument, NULL, 'y'},
            {"simplify-window", required_argument, NULL, 'w'},
            {"blocking-size", required_argument, NULL, 'b'},
            {"non-path-matches", no_argument, NULL, 'P'},
            {"max-count", required_argument, NULL, 'm'},
            {"max-anchors", required_argument, NULL, 'a'},
            {"count-power", required_argument, NULL, 'p'},
            {"chain-alg", required_argument, NULL, 'g'},
            {"no-unaln", no_argument, NULL, 'u'},
            {"verbosity", required_argument, NULL, 'v'},
            {"config", required_argument, NULL, 'C'},
            {"restart", no_argument, NULL, 'R'},
            {"help", no_argument, NULL, 'h'},
            {"skip-calibration", no_argument, NULL, opt_skip_calibration},
            {"force-gfa-output", no_argument, NULL, opt_force_gfa_output},
            {"bond-prefix", required_argument, NULL, opt_bond_prefix},
            {NULL, 0, NULL, 0}
        };
        int o = getopt_long(argc, argv, "T:A:S:s:w:cy:b:Pm:a:p:g:uv:C:Rh", options, NULL);
        
        if (o == -1) {
            // end of options
            break;
        }
        switch (o)
        {
            case 'T':
                params.tree_name = optarg;
                break;
            case 'A':
                params.all_pairs_prefix = optarg;
                break;
            case 'S':
                params.subproblems_prefix = optarg;
                break;
            case 's':
                params.subalignments_filepath = optarg;
                break;
            case 'c':
                params.cyclize_tandem_duplications = true;
                break;
            case 'y':
                params.min_cyclizing_length = parse_int(optarg);
                break;
            case 'w':
                params.simplify_window = parse_int(optarg);
                break;
            case 'b':
                params.blocking_allele_size = parse_int(optarg);
                break;
            case 'P':
                params.path_matches = false;
                break;
            case 'm':
                params.max_count = parse_int(optarg);
                break;
            case 'a':
                params.max_num_match_pairs = parse_int(optarg);
                break;
            case 'p':
                params.pair_count_power = parse_double(optarg);
                break;
            case 'g':
                params.chaining_algorithm = (Anchorer::ChainAlgorithm) parse_int(optarg);
                break;
            case 'u':
                params.constraint_method = Partitioner::Null;
                break;
            case 'v':
                params.logging_level = (logging::LoggingLevel) parse_int(optarg);
                break;
            case 'C':
                config_name = optarg;
                break;
            case 'R':
                restart = true;
                break;
            case 'h':
                print_help();
                return 0;
            case opt_skip_calibration:
                params.skip_calibration = true;
                break;
            case opt_force_gfa_output:
                force_gfa_output = true;
                break;
            case opt_bond_prefix:
                params.bonds_prefix = optarg;
                break;
            default:
                print_help();
                return 1;
        }
    }
    
    const int min_num_positional = 0;
    const int max_num_positional = 1;
    if (argc - optind < min_num_positional || argc - optind > max_num_positional) {
        cerr << "error: expected between " << min_num_positional << " and " << max_num_positional << " positional arguments but got " << (argc - optind) << "\n";
        print_help();
        return 1;
    }
    
    if (argc - optind == 1) {
        params.fasta_name = argv[optind++];
    }
    
    if (!config_name.empty()) {
        
        Parameters defaults;
        if (params != defaults) {
            cerr << "warning: All other command-line arguments are being overridden by config file parameters.\n\n";
        }
        
        ifstream config_file;
        istream* config_stream = get_input(config_name, config_file);
        params = Parameters(*config_stream);
    }
    
    
    // make sure the parameters are all good
    try {
        params.validate();
    }
    catch (exception& ex) {
        cerr << "error: " << ex.what() << '\n';
        print_help();
        return 1;
    }
      
    logging::level = params.logging_level;
    
    if (restart && params.subproblems_prefix.empty()) {
        cerr << "error: cannot restart without prefix for saved subproblems\n";
        return 1;
    }
    
    {
        stringstream strm;
        strm << "Executing command:";
        for (int i = 0; i < argc; ++i) {
            strm << ' ' << argv[i];
        }
        logging::log(logging::Minimal, strm.str());
    }
    
    logging::log(logging::Minimal, "Centrolign version: " + to_string(Version::MAJOR) + "." + to_string(Version::MINOR) + "." + to_string(Version::PATCH));
    logging::log(logging::Minimal, "Centrolign commit: " + Version::GIT_HASH);
    
    {
        stringstream strm;
        strm << "config file:\n\n";
        strm << params.generate_config() << '\n';
        logging::log(logging::Debug, strm.str());
    }
    
    logging::log(logging::Debug, "Baseline memory usage is " + format_memory_usage(current_memory_usage()) + ".");
    
    ifstream fasta_file, tree_file;
    istream* fasta_stream = nullptr;
    istream* tree_stream = nullptr;
    fasta_stream = get_input(params.fasta_name, fasta_file);
    if (!params.tree_name.empty()) {
        tree_stream = get_input(params.tree_name, tree_file);
    }
    
    logging::log(logging::Basic, "Reading input.");
    
    vector<pair<string, string>> parsed = parse_fasta(*fasta_stream);
    
    if (parsed.size() < 2) {
        cerr << "error: FASTA input from " << (params.fasta_name == "-" ? string("STDIN") : params.fasta_name) << " contains " << parsed.size() << " sequence(s), cannot form an alignment\n";
        return 1;
    }
    
    // remember the sequence names, even though we'll give them away to the Core soon
    // TODO: ugly
    vector<string> seq_names;
    for (const auto& name_and_seq : parsed) {
        seq_names.push_back(name_and_seq.first);
    }
    
    string newick_string;
    if (tree_stream == nullptr) {
        // make a dummy Newick string
        if (seq_names.size() > 2) {
            cerr << "warning: it is highly recommended to provide a guide tree (-T) when aligning > 2 sequences\n";
        }
        newick_string = move(in_order_newick_string(parsed));
    }
    else {
        // read it from the file
        stringstream sstrm;
        sstrm << tree_stream->rdbuf();
        newick_string = move(sstrm.str());
    }
    Tree tree(newick_string);
    
    Core core(move(parsed), move(tree));
    
    // pass through parameters
    if (seq_names.size() == 2) {
        // we need this for the CIGAR output TODO: ugly
        params.preserve_subproblems = true;
    }
    params.apply(core);
    
    if (!core.subalignments_filepath.empty()) {
        if (ifstream(core.subalignments_filepath).good()) {
            throw runtime_error("Subalignment file already exists: " + core.subalignments_filepath);
        }
    }
    if (restart) {
        core.restart();
    }
            
    // do the alignment
    core.execute();
    
    if (seq_names.size() == 2 && !force_gfa_output) {
        // output a CIGAR
        const auto& root = core.root_subproblem();
        const auto& leaf1 = core.leaf_subproblem(seq_names.front());
        const auto& leaf2 = core.leaf_subproblem(seq_names.back());
        cout << explicit_cigar(root.alignment, leaf1.graph, leaf2.graph) << '\n';
    }
    else {
        // output a GFA
        const auto& root = core.root_subproblem();
        write_gfa(root.graph, root.tableau, cout);
    }
    
    int64_t max_mem = max_memory_usage();
    
    if (max_mem < 0) {
        logging::log(logging::Basic, "Failed to measure memory usage.");
    }
    else {
        logging::log(logging::Basic, "Maximum memory usage: " + format_memory_usage(max_mem) + ".");
    }
    
    logging::log(logging::Minimal, "Run completed successfully, exiting.");
    
    return 0;
}
