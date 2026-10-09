#include <unordered_map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <getopt.h>

#include "centrolign/utility.hpp"
#include "centrolign/modify_graph.hpp"
#include "centrolign/match_finder.hpp"


using namespace std;
using namespace centrolign;

void print_help() {
    cerr << "usage:\n";
    cerr << "dump_mrms <fasta> > mrms.tsv\n";
}


int main(int argc, char* argv[]) {
    
    // Choose a score function that is always positive so nothing gets thrown out
    ScoreFunction score_function;
    score_function.anchor_score_function = ScoreFunction::InverseCount;
    
    PathMatchFinder match_finder(score_function);
    match_finder.max_count = 40;
    
    while (true)
    {
        static struct option options[] = {
            {"max-count", required_argument, NULL, 'm'},
            {"help", no_argument, NULL, 'h'},
            {NULL, 0, NULL, 0}
        };
        int o = getopt_long(argc, argv, "m:h", options, NULL);
        
        if (o == -1) {
            // end of uptions
            break;
        }
        switch (o)
        {
            case 'm':
                match_finder.max_count = atoi(optarg);
                break;
            case 'h':
                print_help();
                return 0;
            default:
                print_help();
                return 1;
        }
    }

    if (argc - optind != 1) {
        cerr << "error: expected 1 positional argument\n";
        print_help();
        return 1;
    }
    
    ifstream fasta_in(argv[optind]);
    if (!fasta_in) {
        cerr << "error: could not open FASTA file " << argv[optind] << '\n';
        return 1;
    }
    
    auto seqs = parse_fasta(fasta_in);
    
    if (seqs.size() != 2) {
        cerr << "FASTA must contain 2 sequences\n";
        return 1;
    }
    
    auto graph1 = make_base_graph(seqs.front().first, seqs.front().second);
    auto graph2 = make_base_graph(seqs.back().first, seqs.back().second);
    auto tableau1 = add_sentinels(graph1, 5, 6);
    auto tableau2 = add_sentinels(graph1, 7, 8);
    
    auto matches = match_finder.find_matches(graph1, graph2, tableau1, tableau2);
    
    std::stable_sort(matches.begin(), matches.end(), [](const match_set_t& a, const match_set_t& b) {
        auto c1 = a.walks1.size() * a.walks2.size();
        auto c2 = b.walks1.size() * b.walks2.size();
        return c1 < c2 || (c1 == c2 && a.walks1.front().size() > b.walks1.front().size());
    });
    
    for (const auto& match_set : matches) {
        for (const auto& walk1 : match_set.walks1) {
            for (const auto& walk2 : match_set.walks2) {
                std::cout << walk1.front() << '\t' << walk2.front() << '\t' << walk1.size() << '\n';
            }
        }
    }
}
