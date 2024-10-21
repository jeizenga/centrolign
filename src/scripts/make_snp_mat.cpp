#include <unordered_map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <tuple>
#include <utility>
#include <getopt.h>

#include "centrolign/gfa.hpp"
#include "centrolign/graph.hpp"
#include "centrolign/modify_graph.hpp"
#include "centrolign/superbubbles.hpp"


using namespace std;
using namespace centrolign;

void print_help() {
    cerr << "usage:\n";
    cerr << "make_snp_mat [options] graph.gfa > snp_mat.tsv\n\n";
    cerr << "options:\n";
    cerr << " --base / -b   Use bases in the output encoding\n";
    cerr << " --help / -h   Print this message and exit\n";
}

int main(int argc, char* argv[]) {

    bool output_base = false;
    
    while (true)
    {
        static struct option options[] = {
            {"base", no_argument, NULL, 'b'},
            {"help", no_argument, NULL, 'h'},
            {NULL, 0, NULL, 0}
        };
        int o = getopt_long(argc, argv, "hb", options, NULL);
        
        if (o == -1) {
            // end of uptions
            break;
        }
        switch (o)
        {
            case 'b':
                output_base = true;
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
        cerr << "error: expected 1 positional arguments but got " << (argc - optind) << "\n";
        print_help();
        return 1;
    }
    
    
    ifstream gfa_in(argv[optind]);
    if (!gfa_in) {
        cerr << "error: could not open GFA file " << argv[optind] << '\n';
        return 1;
    }
    
    BaseGraph graph = read_gfa(gfa_in, false);
    auto tableau = add_sentinels(graph, '^', '$');
    
    SuperbubbleTree bub_tree(graph, tableau);
    
    unordered_map<uint64_t, size_t> snp_starts;
    
    for (uint64_t bub_id = 0; bub_id < bub_tree.structure_size(); ++bub_id) {
        
        if (!bub_tree.chains_inside(bub_id).empty()) {
            // not a leaf bubble
            continue;
        }
        
        uint64_t src, snk;
        tie(src, snk) = bub_tree.structure_boundaries(bub_id);
        
        const auto& next = graph.next(src);
        if (graph.next(src).size() != 2) {
            // not biallelic
            continue;
        }
        
        bool skip = false;
        for (auto allele : next) {
            const auto& next_next = graph.next(allele);
            if (next_next.size() != 1 || next_next.front() != snk) {
                skip = true;
            }
        }
        
        if (skip) {
            // not a snp
            continue;
        }
        
        // record a snp
        snp_starts[src] = snp_starts.size();
    }
    
    // header row expected by phylip
    cout << graph.path_size() << '\t' << snp_starts.size() << '\n';
    for (uint64_t path_id = 0; path_id < graph.path_size(); ++path_id) {
        
        // init with missing values
        string row(snp_starts.size(), '?');
        
        const auto& path = graph.path(path_id);
        for (size_t i = 0; i < path.size(); ++i) {
            auto it = snp_starts.find(path[i]);
            if (it != snp_starts.end()) {
                // enter the allele of this snp
                if (output_base) {
                    row[it->second] = graph.label(path[i + 1]);
                }
                else {
                    row[it->second] = (path[i + 1] == graph.next(path[i]).front()) ? '0' : '1';
                }
            }
        }
        
        // output the row
        cout << graph.path_name(path_id);
        for (auto c : row) {
            cout << '\t' << c;
        }
        cout << '\n';
    }
}
