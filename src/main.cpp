#include <getopt.h>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <fstream>

#include "centrolign/utility.hpp"
#include "centrolign/anchorer.hpp"
#include "centrolign/modify_graph.hpp"
#include "centrolign/chain_merge.hpp"

using namespace std;
using namespace centrolign;

/*
 * the default parameters
 */
struct Defaults {
    static const int64_t max_count = 50;
    static const int64_t max_num_match_pairs = 10000;
};

void print_help() {
    cerr << "Usage: centrolign [options] sequences.fasta\n";
    cerr << "Options:\n";
    cerr << " --max-count / -m INT     The maximum number of times an anchor can occur [" << Defaults::max_count << "]\n";
    cerr << " --max-anchors / -a INT   The maximum number of anchors [" << Defaults::max_num_match_pairs << "]\n";
    cerr << " --help / -h              Print this message and exit\n";
}

int main(int argc, char** argv) {
    
    int64_t max_count = Defaults::max_count;
    int64_t max_num_match_pairs = Defaults::max_num_match_pairs;
    
    while (true)
    {
        static struct option options[] = {
            {"max-count", required_argument, NULL, 'm'},
            {"max-anchors", required_argument, NULL, 'a'},
            {"help", no_argument, NULL, 'h'},
            {NULL, 0, NULL, 0}
        };
        int o = getopt_long(argc, argv, "hm:a:", options, NULL);
        
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
    
    string fasta_name = argv[optind++];
    
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
    
    BaseGraph graph1 = make_base_graph(parsed.front().first, parsed.front().second);
    BaseGraph graph2 = make_base_graph(parsed.back().first, parsed.back().second);
    
    SentinelTableau sentinels1 = add_sentinels(graph1, 5, 6);
    SentinelTableau sentinels2 = add_sentinels(graph2, 7, 8);
    
    ChainMerge chain_merge1(graph1);
    ChainMerge chain_merge2(graph2);
    
    Anchorer anchorer;
    anchorer.max_count = max_count;
    anchorer.max_num_match_pairs = max_num_match_pairs;
    
    auto anchors = anchorer.anchor_chain(graph1, graph2, chain_merge1, chain_merge2);
    
    cout << "idx1" << '\t' << "idx2" << '\t' << "len" << '\t' << "cnt1" << '\t' << "cnt2" << '\n';
    for (const auto& anchor : anchors) {
        cout << anchor.walk1.front() << '\t' << anchor.walk2.front() << '\t' << anchor.walk1.size() << '\t' << anchor.count1 << '\t' << anchor.count2 << '\n';
    }
    
    return 0;
}
