#include <unordered_map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <tuple>
#include <utility>

#include "centrolign/gfa.hpp"
#include "centrolign/graph.hpp"
#include "centrolign/modify_graph.hpp"
#include "centrolign/superbubbles.hpp"


using namespace std;
using namespace centrolign;

unordered_map<char, char> translator{
    {'0', 'z'},
    {'1', 'o'},
    {'2', 't'},
    {'3', 'h'},
    {'4', 'f'},
    {'5', 'i'},
    {'6', 's'},
    {'7', 'v'},
    {'8', 'e'},
    {'9', 'n'},
    {'.', 'p'}
};

string translate_sample_name(const string& name) {
    string translated = name;
    for (size_t i = 0; i < translated.size(); ++i) {
        auto it = translator.find(translated[i]);
        if (it != translator.end()) {
            translated[i] = it->second;
        }
    }
    return translated;
}

int main(int argc, char* argv[]) {
    
    if (argc != 2) {
        cerr << "usage:\n";
        cerr << "make_snp_mat graph.gfa > snp_mat.tsv\n";
        return 1;
    }
    
    ifstream gfa_in(argv[1]);
    if (!gfa_in) {
        cerr << "error: could not open GFA file " << argv[1] << '\n';
        return 1;
    }
    
    BaseGraph graph = read_gfa(gfa_in, false);
    auto tableau = add_sentinels(graph, '^', '$');
    
    SuperbubbleTree bub_tree(graph, tableau);
    
    unordered_map<uint64_t, size_t> snp_starts;
    
    for (uint64_t bub_id = 0; bub_id < bub_tree.superbubble_size(); ++bub_id) {
        
        if (!bub_tree.chains_inside(bub_id).empty()) {
            // not a leaf bubble
            continue;
        }
        
        uint64_t src, snk;
        tie(src, snk) = bub_tree.superbubble_boundaries(bub_id);
        
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
                row[it->second] = (path[i + 1] == graph.next(path[i]).front()) ? '0' : '1';
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
