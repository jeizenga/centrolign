#include <getopt.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>

#include "centrolign/utility.hpp"
#include "centrolign/anchorer.hpp"
#include "centrolign/modify_graph.hpp"
#include "centrolign/chain_merge.hpp"

using namespace std;
using namespace centrolign;

struct Defaults {
    static const int max_count = 50;
};

void print_help() {
    cerr << "Usage: centrolign [options] sequences.fasta\n";
    cerr << "Options:\n";
    cerr << " --max-count / -m INT   The maximum number of times an anchor can occur [" << Defaults::max_count << "]\n";
    cerr << " --help / -h            Print this message and exit\n";
}

int main(int argc, char** argv) {
    
    int max_count = Defaults::max_count;
    
    int c;
    while (true)
    {
        static struct option options[] =
        {
            {"max-count", required_argument, 0, 'm'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };
        c = getopt_long(argc, argv, "hm:", options, NULL);
        
        if (c == -1) {
            // end of uptions
            break;
        }
        
        switch (c)
        {
            case 'm':
                max_count = parse_int(optarg);
                break;
            case 'h':
                print_help();
                return 0;
            default:
                print_help();
                return 1;
        }
    }
    
    int num_positional = 1;
    if (optind + num_positional != argc) {
        cerr << "error: expected " << num_positional << " positional argument but got " << (argc - optind) << '\n';
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
    
    return 0;
}
