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
#include "centrolign/structure_distances.hpp"
#include "centrolign/snarls.hpp"
#include "centrolign/utility.hpp"


using namespace std;
using namespace centrolign;

void print_help() {
    cerr << "usage:\n";
    cerr << "make_snp_mat [options] graph.gfa > snp_mat.tsv\n\n";
    cerr << "options:\n";
    cerr << " --base / -b        Use bases in the output encoding\n";
    cerr << " --indels / -i      Include point indels (< --sv-lim)\n";
    cerr << " --svs / -s         Include biallelic structural variants (>= --sv-lim)\n";
    cerr << " --sv-lim / -l INT  The size at which a variant is considered a structural variant [50]\n";
    cerr << " --no-header / -n   Do not include the Phyllip header line\n";
    cerr << " --help / -h        Print this message and exit\n";
}

int main(int argc, char* argv[]) {

    bool output_base = false;
    bool include_indels = false;
    bool include_svs = false;
    bool include_header = true;
    int64_t sv_lim = 50;
    
    while (true)
    {
        static struct option options[] = {
            {"base", no_argument, NULL, 'b'},
            {"indels", no_argument, NULL, 'i'}
            {"svs", no_argument, NULL, 's'},
            {"sv-lim", required_argument, NULL, 'l'},
            {"no-header", no_argument, NULL, 'n'},
            {"help", no_argument, NULL, 'h'},
            {NULL, 0, NULL, 0}
        };
        int o = getopt_long(argc, argv, "hbinsl:", options, NULL);
        
        if (o == -1) {
            // end of uptions
            break;
        }
        switch (o)
        {
            case 'b':
                output_base = true;
                break;
            case 'i':
                include_indels = true;
                break;
            case 's':
                include_svs = true;
                break;
            case 'l':
                sv_lim = parse_int(optarg);
                break;
            case 'n':
                include_header = false;
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
    
    SnarlTree snarl_tree(graph, tableau);
    SnarlDistances snarl_distances(snarl_tree, graph);
    
    typedef enum {
        SNP,
        PointIndel,
        SV
    } var_type_t;
    
    std::vector<std::pair<uint64_t, var_type_t>> variants;
    
    // trivial == contains only trivial snarls
    // simple == contains only point variants
    std::vector<bool> chain_is_trivial(snarl_tree.chain_size());
    std::vector<bool> chain_is_simple(snarl_tree.chain_size());
    std::vector<bool> snarl_is_trivial(snarl_tree.chain_size());
    std::vector<bool> snarl_is_simple(snarl_tree.structure_size());
    std::vector<bool> snarl_is_biallelic(snarl_tree.structure_size());
    
    for (auto feature : snarl_tree.postorder()) {
        if (feature.second) {
            // chain
            bool trivial = true;
            bool simple = true;
            for (auto snarl_id : snarl_tree.structures_inside(feature.first)) {
                trivial = trivial && snarl_is_trivial[snarl_id];
                simple = simple && snarl_is_simple[snarl_id];
            }
            chain_is_trivial[feature.first] = trivial;
            chain_is_simple[feature.first] = simple;
        }
        else {
            // snarl
            
            if (!snarl_tree.snarl_is_acyclic(feature.first)) {
                snarl_is_trivial[feature.first] = false;
                snarl_is_simple[feature.first] = false;
                snarl_is_biallelic[feature.first] = false;
                continue;
            }
            
            uint64_t src, snk;
            std::tie(src, snk) = snarl_tree.structure_boundaries(feature.first);
            if (graph.next_size(src) == 1 && graph.next(src).front() == snk) {
                snarl_is_trivial[feature.first] = true;
            }
            else {
                snarl_is_trivial[feature.first] = false;
            }
            
            size_t min_dist, max_dist;
            std::tie(min_dist, max_dist) = snarl_distances.structure_min_max_dist(feature.first);
            snarl_is_simple[feature.first] = (max_dist != -1 && max_dist < sv_lim);
            
            NetGraph net_graph(graph, snarl_tree, feature.first);
            // find the source and sink node
            uint64_t net_src_id = -1, net_snk_id = -1;
            for (uint64_t net_id = 0; net_id < net_graph.node_size(); ++net_id) {
                if (net_graph.label(net_id) == std::make_pair(src, false)) {
                    net_src_id = net_id;
                }
                if (net_graph.label(net_id) == std::make_pair(snk, false)) {
                    net_snk_id = net_id;
                }
            }
            assert(net_src_id != -1 && net_snk_id != -1);
            
            bool biallelic = true;
            if (net_graph.next_size(net_src_id) != 2) {
                biallelic = false;
            }
            else {
                for (auto next_id : net_graph.next(net_src_id)) {
                    if (next_id != net_snk_id) {
                        // this isn't a deletion allele
                        uint64_t next_feature_id;
                        bool is_chain;
                        std::tie(next_feature_id, is_chain) = net_graph.label(next_id);
                        if (is_chain) {
                            if (!chain_is_trivial[next_feature_id]) {
                                // there are variants on this chain, which create more than 2 alleles for this site
                                biallelic = false;
                            }
                        }
                        if (net_graph.next_size(next_id) != 1 || net_graph.next(next_id).front() != net_snk_id) {
                            // there is further branching beyond this feature
                            biallelic = false;
                        }
                    }
                }
            }
            snarl_is_biallelic[feature.first] = biallelic;
            
        }
    }
    
    for (uint64_t bub_id = 0; bub_id < snarl_tree.structure_size(); ++bub_id) {
        
        if (!snarl_tree.chains_inside(bub_id).empty()) {
            // not a leaf bubble
            continue;
        }
        
        uint64_t src, snk;
        tie(src, snk) = snarl_tree.structure_boundaries(bub_id);
        
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
    
    if (include_header) {
        // header row expected by phylip
        cout << graph.path_size() << '\t' << snp_starts.size() << '\n';
    }
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
