#include <getopt.h>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>

#include "centrolign/utility.hpp"
#include "centrolign/anchorer.hpp"
#include "centrolign/modify_graph.hpp"
#include "centrolign/chain_merge.hpp"
#include "centrolign/stitcher.hpp"
#include "centrolign/logging.hpp"

using namespace std;
using namespace centrolign;

/*
 * the default parameters
 */
struct Defaults {
    static const int64_t max_count = 100;
    static const int64_t max_num_match_pairs = 5000000;
    static constexpr double pair_count_power = 1.0;
    static const bool length_scale = true;
    static const bool sparse_chaining = true;
    static const bool output_anchors = false;
};

void print_help() {
    cerr << "Usage: centrolign [options] sequences.fasta\n";
    cerr << "Options:\n";
    cerr << " --max-count / -m INT      The maximum number of times an anchor can occur [" << Defaults::max_count << "]\n";
    cerr << " --max-anchors / -a INT    The maximum number of anchors [" << Defaults::max_num_match_pairs << "]\n";
    cerr << " --count-power / -p FLOAT  Scale anchor weights by the count raised to this power [" << Defaults::pair_count_power << "]\n";
    cerr << " --no-length-scale / -l    Do not scale anchor weights by length\n";
    cerr << " --no-sparse-chain / -s    Do not use sparse chaining algorithm\n";
    cerr << " --output-anchors / -A     Output the anchoring results as a table\n";
    cerr << " --help / -h               Print this message and exit\n";
}

int main(int argc, char** argv) {
    
    // init the local params with the defaults
    int64_t max_count = Defaults::max_count;
    int64_t max_num_match_pairs = Defaults::max_num_match_pairs;
    double pair_count_power = Defaults::pair_count_power;
    bool length_scale = Defaults::length_scale;
    bool sparse_chaining = Defaults::sparse_chaining;
    bool output_anchors = Defaults::output_anchors;
    
    while (true)
    {
        static struct option options[] = {
            {"max-count", required_argument, NULL, 'm'},
            {"max-anchors", required_argument, NULL, 'a'},
            {"count-power", required_argument, NULL, 'p'},
            {"no-length-scale", no_argument, NULL, 'l'},
            {"no-sparse-chain", no_argument, NULL, 's'},
            {"output-anchors", no_argument, NULL, 'A'},
            {"help", no_argument, NULL, 'h'},
            {NULL, 0, NULL, 0}
        };
        int o = getopt_long(argc, argv, "m:a:p:lsAh", options, NULL);
        
        if (o == -1) {
            // end of uptions
            break;
        }
        
        switch (o)
        {
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
            case 'A':
                output_anchors = true;
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
    
    logging::log(logging::Basic, "reading input...");
    vector<pair<string, string>> parsed;
    if (fasta_name == "-") {
        // read piped input
        parsed = parse_fasta(cin);
    }
    else {
        // open a file
        ifstream fasta_in(fasta_name);
        if (!fasta_in) {
            cerr << "error: could not open FASTA file " << fasta_name << '\n';
            return 1;
        }
        parsed = parse_fasta(fasta_in);
    }
    
    if (parsed.size() != 2) {
        cerr << "error: this prototype only supports pairwise sequence alignment\n";
        return 1;
    }
    
    logging::log(logging::Basic, "building graphs...");
    BaseGraph graph1 = make_base_graph(parsed.front().first, parsed.front().second);
    BaseGraph graph2 = make_base_graph(parsed.back().first, parsed.back().second);
    
    SentinelTableau sentinels1 = add_sentinels(graph1, 5, 6);
    SentinelTableau sentinels2 = add_sentinels(graph2, 7, 8);
    
    logging::log(logging::Basic, "computing reachability...");
    ChainMerge chain_merge1(graph1);
    ChainMerge chain_merge2(graph2);
    
    logging::log(logging::Basic, "anchoring...");
    Anchorer anchorer;
    anchorer.max_count = max_count;
    anchorer.pair_count_power = pair_count_power;
    anchorer.length_scale = length_scale;
    anchorer.max_num_match_pairs = max_num_match_pairs;
    anchorer.sparse_chaining = sparse_chaining;
    
    auto anchors = anchorer.anchor_chain(graph1, graph2, chain_merge1, chain_merge2);
    
    if (output_anchors) {
        cout << "idx1" << '\t' << "idx2" << '\t' << "len" << '\t' << "cnt1" << '\t' << "cnt2" << '\n';
        for (const auto& anchor : anchors) {
            cout << anchor.walk1.front() << '\t' << anchor.walk2.front() << '\t' << anchor.walk1.size() << '\t' << anchor.count1 << '\t' << anchor.count2 << '\n';
        }
    }
    else {
        
        logging::log(logging::Basic, "stitching anchors into an alignment...");
        
        Stitcher stitcher;
        
        Alignment alignment = stitcher.stitch(anchors, graph1, graph2,
                                              sentinels1, sentinels2,
                                              chain_merge1, chain_merge2);
        
        cout << explicit_cigar(alignment, graph1, graph2) << '\n';
        
    }
    
    logging::log(logging::Minimal, "run completed successfully, exiting.");
    
    return 0;
}
