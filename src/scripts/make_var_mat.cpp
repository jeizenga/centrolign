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
    cerr << "make_var_mat [options] graph.gfa > var_mat.tsv\n\n";
    cerr << "options:\n";
    cerr << " --base / -b        Use bases in the output encoding\n";
    cerr << " --indels / -i      Include point indels (< --sv-lim)\n";
    cerr << " --mnvs / -m        Include multi-nucleotide variants (< --sv-lim)\n";
    cerr << " --svs / -s         Include structural variants (>= --sv-lim)\n";
    cerr << " --sv-lim / -l INT  The size at which a variant is considered a structural variant [50]\n";
    cerr << " --no-header / -n   Do not include the Phyllip header line\n";
    cerr << " --help / -h        Print this message and exit\n";
}

int main(int argc, char* argv[]) {

    bool output_base = false;
    bool include_indels = false;
    bool include_svs = false;
    bool include_mnvs = false;
    bool include_header = true;
    int64_t sv_lim = 50;
    
    while (true)
    {
        static struct option options[] = {
            {"base", no_argument, NULL, 'b'},
            {"indels", no_argument, NULL, 'i'},
            {"mnvs", no_argument, NULL, 'm'},
            {"svs", no_argument, NULL, 's'},
            {"sv-lim", required_argument, NULL, 'l'},
            {"no-header", no_argument, NULL, 'n'},
            {"help", no_argument, NULL, 'h'},
            {NULL, 0, NULL, 0}
        };
        int o = getopt_long(argc, argv, "hbinsml:", options, NULL);
        
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
            case 'm':
                include_mnvs = true;
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
    
    std::cerr << "Loading graph...\n";
    BaseGraph graph = read_gfa(gfa_in, false);
    auto tableau = add_sentinels(graph, '^', '$');
    
    std::cerr << "Finding snarls...\n";
    SnarlTree snarl_tree(graph, tableau);
    
    std::cerr << "Computing snarl sizes...\n";
    SnarlDistances snarl_distances(snarl_tree, graph);
    
    // trivial == contains only trivial snarls
    // simple == contains only point variants
    // biallelic == only two paths through it
    
    std::cerr << "Identifying snarl features...\n";
    std::vector<bool> chain_is_trivial(snarl_tree.chain_size());
    std::vector<bool> chain_is_simple(snarl_tree.chain_size());
    
    std::vector<bool> snarl_is_trivial(snarl_tree.structure_size());
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
    std::cerr << "Selecting variants...\n";
    
    typedef enum {
        SNP,
        PointIndel,
        MNV,
        SV
    } var_type_t;
    
    // snarl id and type
    std::vector<std::pair<uint64_t, var_type_t>> variants;
    for (uint64_t snarl_id = 0; snarl_id < snarl_tree.structure_size(); ++snarl_id) {
        
        if (snarl_is_biallelic[snarl_id] && !snarl_is_trivial[snarl_id]) {
                        
            size_t min_dist, max_dist;
            std::tie(min_dist, max_dist) = snarl_distances.structure_min_max_dist(snarl_id);
            
            if (min_dist == max_dist && min_dist == 3) {
                variants.emplace_back(snarl_id, SNP);
            }
            else if (min_dist == 2 && max_dist < sv_lim) {
                variants.emplace_back(snarl_id, PointIndel);
            }
            else if (max_dist < sv_lim) {
                variants.emplace_back(snarl_id, MNV);
            }
            else {
                variants.emplace_back(snarl_id, SV);
            }
        }
    }
    
    // choose the variants we want and assign them to a column
    std::unordered_map<uint64_t, std::pair<uint64_t, size_t>> source_to_column;
    for (const auto& var : variants) {
        if (var.second == SNP ||
            (var.second == PointIndel && include_indels) ||
            (var.second == MNV && include_mnvs) ||
            (var.second == SV && include_svs)) {
            
            uint64_t src, snk;
            std::tie(src, snk) = snarl_tree.structure_boundaries(var.first);
            source_to_column[src] = std::make_pair(snk, source_to_column.size());
        }
    }
    std::cerr << "Outputting table...\n";
    
    if (include_header) {
        // header row expected by phylip
        cout << graph.path_size() << '\t' << source_to_column.size() << '\n';
    }
//    for (size_t i = 0; i < source_to_column.size(); ++i) {
//        std::cout << '\t' << i;
//    }
//    std::cout << '\n';
    for (uint64_t path_id = 0; path_id < graph.path_size(); ++path_id) {
        
        // init with missing values
        std::vector<std::vector<std::string>> row(source_to_column.size());
        
        const auto& path = graph.path(path_id);
        std::vector<std::pair<uint64_t, size_t>> curr_vars;
        for (size_t i = 0; i < path.size(); ++i) {
            if (!curr_vars.empty() && curr_vars.back().first == path[i]) {
                // we're exiting this snarl
                curr_vars.pop_back();
            }
            
            for (const auto& var : curr_vars) {
                row[var.second].back().push_back(graph.label(path[i]));
            }
            
            auto it = source_to_column.find(path[i]);
            if (it != source_to_column.end()) {
                if (!output_base) {
                    // just record and ID for the allele
                    // note: this will break if we go to non-biallelic variants
                    row[it->second.second].push_back((path[i + 1] == graph.next(path[i]).front()) ? "0" : "1");
                }
                else {
                    // set up to start including the base sequence
                    row[it->second.second].emplace_back();
                    curr_vars.emplace_back(it->second);
                }
            }
        }
        
        // output the row
        cout << graph.path_name(path_id);
        for (const auto& alleles : row) {
            cout << '\t';
            if (alleles.empty()) {
                cout << '?';
            }
            else {
                for (size_t i = 0; i < alleles.size(); ++i) {
                    if (i) {
                        cout << ',';
                    }
                    if (alleles[i].empty()) {
                        cout << '-';
                    }
                    else {
                        cout << alleles[i];
                    }
                }
            }
        }
        cout << '\n';
    }
}
