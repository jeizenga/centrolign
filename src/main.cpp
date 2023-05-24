#include <getopt.h>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>

#include "centrolign/utility.hpp"
#include "centrolign/core.hpp"
#include "centrolign/logging.hpp"
#include "centrolign/gfa.hpp"

using namespace std;
using namespace centrolign;

/*
 * the default parameters
 */
struct Defaults {
    static const int64_t simplify_window = 10000;
    static const int64_t max_walk_count = 8;
    static const int64_t blocking_allele_size = 32;
    static const int64_t max_count = 100;
    static const int64_t max_num_match_pairs = 5000000;
    static constexpr double pair_count_power = 1.0;
    static const bool length_scale = true;
    static const bool sparse_chaining = true;
    static const logging::LoggingLevel logging_level = logging::Basic;
};

void print_help() {
    cerr << "Usage: centrolign [options] sequences.fasta\n";
    cerr << "Options:\n";
    cerr << " --tree / -T FILE          Newick format guide tree for alignment [in FASTA order]\n";
    cerr << " --simplify-window / -w    Size window used for graph simplification [" << Defaults::simplify_window << "]\n";
    cerr << " --simplify-count / -c     Number of walks through window that trigger simplification [" << Defaults::max_walk_count << "]\n";
    cerr << " --blocking-size / -b      Minimum size allele to block simplification [" << Defaults::blocking_allele_size << "]\n";
    cerr << " --max-count / -m INT      The maximum number of times an anchor can occur [" << Defaults::max_count << "]\n";
    cerr << " --max-anchors / -a INT    The maximum number of anchors [" << Defaults::max_num_match_pairs << "]\n";
    cerr << " --count-power / -p FLOAT  Scale anchor weights by the count raised to this power [" << Defaults::pair_count_power << "]\n";
    cerr << " --no-length-scale / -l    Do not scale anchor weights by length\n";
    cerr << " --no-sparse-chain / -s    Do not use sparse chaining algorithm\n";
    cerr << " --verbosity / -v INT      Select from: 0 (silent), 1 (minimal), 2 (basic), 3 (verbose), 4 (debug) [" << (int) Defaults::logging_level << "]\n";
    cerr << " --help / -h               Print this message and exit\n";
}

// make a dummy Newick string for in-order alignment
string in_order_newick_string(const vector<pair<string, string>>& sequences) {
    stringstream strm;
    for (size_t i = 1; i < sequences.size(); ++i) {
        strm << '(';
    }
    strm << sequences.front().first;
    for (size_t i = 1; i < sequences.size(); ++i) {
        strm << ',' << sequences[i].first << ')';
    }
    strm << ';';
    return strm.str();
}

int main(int argc, char** argv) {
        
    // init the local params with the defaults
    int64_t max_count = Defaults::max_count;
    int64_t max_num_match_pairs = Defaults::max_num_match_pairs;
    double pair_count_power = Defaults::pair_count_power;
    bool length_scale = Defaults::length_scale;
    bool sparse_chaining = Defaults::sparse_chaining;
    logging::level = Defaults::logging_level;
    int64_t simplify_window = Defaults::simplify_window;
    int64_t max_walk_count = Defaults::max_walk_count;
    int64_t blocking_allele_size = Defaults::blocking_allele_size;
    
    // files provided by argument
    string tree_name;
    
    while (true)
    {
        static struct option options[] = {
            {"tree", required_argument, NULL, 'T'},
            {"simplify-window", required_argument, NULL, 'w'},
            {"simplify-count", required_argument, NULL, 'c'},
            {"blocking-size", required_argument, NULL, 'b'},
            {"max-count", required_argument, NULL, 'm'},
            {"max-anchors", required_argument, NULL, 'a'},
            {"count-power", required_argument, NULL, 'p'},
            {"no-length-scale", no_argument, NULL, 'l'},
            {"no-sparse-chain", no_argument, NULL, 's'},
            {"verbosity", required_argument, NULL, 'v'},
            {"help", no_argument, NULL, 'h'},
            {NULL, 0, NULL, 0}
        };
        int o = getopt_long(argc, argv, "T:w:c:b:m:a:p:lsv:h", options, NULL);
        
        if (o == -1) {
            // end of uptions
            break;
        }
        switch (o)
        {
            case 'T':
                tree_name = optarg;
                break;
            case 'w':
                simplify_window = parse_int(optarg);
                break;
            case 'c':
                max_walk_count = parse_int(optarg);
                break;
            case 'b':
                blocking_allele_size = parse_int(optarg);
                break;
            case 'm':
                max_count = parse_int(optarg);
                break;
            case 'a':
                max_num_match_pairs = parse_int(optarg);
                break;
            case 'p':
                pair_count_power = parse_double(optarg);
                break;
            case 'l':
                length_scale = false;
                break;
            case 's':
                sparse_chaining = false;
                break;
            case 'v':
                logging::level = (logging::LoggingLevel) parse_int(optarg);
                break;
            case 'h':
                print_help();
                return 0;
            default:
                print_help();
                return 1;
        }
    }
    
    const int num_positional = 1;
    if (optind + num_positional != argc) {
        cerr << "error: expected " << num_positional << " positional argument but got " << (argc - optind) << "\n\n";
        print_help();
        return 1;
    }
        
    {
        stringstream strm;
        strm << "executing command:";
        for (int i = 0; i < argc; ++i) {
            strm << ' ' << argv[i];
        }
        logging::log(logging::Minimal, strm.str());
    }
    
    string fasta_name = argv[optind++];
    
    ifstream fasta_file, tree_file;
    istream* fasta_stream = nullptr;
    istream* tree_stream = nullptr;
    fasta_stream = get_input(fasta_name, fasta_file);
    if (!tree_name.empty()) {
        tree_stream = get_input(tree_name, tree_file);
    }
    
    logging::log(logging::Basic, "Reading input");
    vector<pair<string, string>> parsed = parse_fasta(*fasta_stream);
    
    if (parsed.size() < 2) {
        cerr << "error: FASTA input from " << (fasta_name == "-" ? string("STDIN") : fasta_name) << " contains " << parsed.size() << " sequences, cannot form an alignment\n";
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
        if (seq_names.size() != 2) {
            cerr << "warning: it is *highly* recommended to provide a guide tree (-T) when aligning >= 2 sequences\n";
        }
        newick_string = in_order_newick_string(parsed);
    }
    else {
        // read it from the file
        stringstream sstrm;
        sstrm << tree_stream->rdbuf();
        newick_string = sstrm.str();
    }
    Tree tree(newick_string);
    
    Core core(move(parsed), move(tree));
    
    // pass through parameters
    core.simplifier.min_dist_window = simplify_window;
    core.simplifier.preserve_bubble_size = blocking_allele_size;
    core.simplifier.max_walks = max_walk_count;
    
    core.match_finder.max_count = max_count;
    core.match_finder.max_num_match_pairs = max_num_match_pairs;
    
    core.anchorer.pair_count_power = pair_count_power;
    core.anchorer.length_scale = length_scale;
    core.anchorer.sparse_chaining = sparse_chaining;
    
    core.preserve_subproblems = true;
    
    // do the alignment
    core.execute();
    
    if (seq_names.size() == 2) {
        // output a CIGAR
        const auto& root = core.root_subproblem();
        const auto& leaf1 = core.leaf_subproblem(seq_names.front());
        const auto& leaf2 = core.leaf_subproblem(seq_names.back());
        cout << explicit_cigar(root.alignment, leaf1.graph, leaf2.graph) << '\n';
    }
    else {
        // output a GFA
        const auto& root = core.root_subproblem();
        write_gfa(root.graph, cout);
    }
    
    logging::log(logging::Minimal, "Run completed successfully, exiting.");
    
    return 0;
}
