#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <tuple>
#include <utility>
#include <sstream>
#include <limits>
#include <cassert>
#include <algorithm>
#include <getopt.h>

#include "centrolign/tree.hpp"
#include "centrolign/utility.hpp"
#include "centrolign/core.hpp"
#include "centrolign/modify_graph.hpp"
#include "centrolign/gfa.hpp"

using namespace std;
using namespace centrolign;

// being naughty
class ExposedCore : public Core {
public:
    ExposedCore(const vector<string>& samples,
                const string& subprob_prefix) : Core(std::move(prepare_dummy_seqs(samples)), std::move(prepare_dummy_tree(samples))) {
        
        this->subproblems_prefix = subprob_prefix;
    }
    
    string root_subproblem_name() const {
        return subproblem_file_name(root_subproblem());
    }
    
private:
    static vector<pair<string, string>> prepare_dummy_seqs(const vector<string>& samples) {
        vector<pair<string, string>> dummy_seqs(samples.size());
        for (size_t i = 0; i < samples.size(); ++i) {
            dummy_seqs[i].first = samples[i];
        }
        return dummy_seqs;
    }
    
    static Tree prepare_dummy_tree(const vector<string>& samples) {
        return Tree(in_order_newick_string(samples));
    }
};


void print_help() {
    cerr << "usage:\n";
    cerr << "remove_sample [options] -s sample [-s sample2 ...] -p output_prefix graph.gfa\n\n";
    cerr << "options:\n";
    cerr << " --prefix / -p PREF      Prefix for graph output (required)\n";
    cerr << " --sample / -s SAMP      Sample to remove from the graph (may repeat, >= 1 required)\n";
    cerr << " --tree-in / -t FILE     Guide tree for the graph in Newick format\n";
    cerr << " --tree-out / -T FILE    Output file for tree with sample(s) rotated to outer branches (requires --tree-in)\n";
    cerr << " --fasta-pref / -f PREF  Prefix for FASTAs containing the removed samples\n";
    cerr << " --help / -h             Print this message and exit\n";
}
int main(int argc, char* argv[]) {
    
    string graph_prefix;
    unordered_set<string> removed_samples;
    string tree_in_file;
    string tree_out_file;
    string fasta_prefix;
    
    while (true)
    {
        static struct option options[] = {
            {"prefix", required_argument, NULL, 'p'},
            {"sample", required_argument, NULL, 's'},
            {"tree-in", required_argument, NULL, 't'},
            {"tree-out", required_argument, NULL, 'T'},
            {"fasta-pref", required_argument, NULL, 'f'},
            {"help", no_argument, NULL, 'h'},
            {NULL, 0, NULL, 0}
        };
        int o = getopt_long(argc, argv, "p:s:t:T:f:h", options, NULL);
        
        if (o == -1) {
            // end of uptions
            break;
        }
        switch (o)
        {
            case 'p':
                graph_prefix = optarg;
                break;
            case 's':
                removed_samples.insert(optarg);
                break;
            case 't':
                tree_in_file = optarg;
                break;
            case 'T':
                tree_out_file = optarg;
                break;
            case 'f':
                fasta_prefix = optarg;
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
        cerr << "error: expected 1 positional argument but got " << (argc - optind) << "\n";
        print_help();
        return 1;
    }
    
    if (graph_prefix.empty()) {
        cerr << "error: --prefix is required\n";
        print_help();
        return 1;
    }
    
    if (!tree_out_file.empty() && tree_in_file.empty()) {
        cerr << "error: --tree-out requires --tree-in to be provided\n";
        print_help();
        return 1;
    }
    
    if (tree_out_file.empty() && !tree_in_file.empty()) {
        cerr << "warning: --tree-in is unused without --tree-out\n";
    }
    
    string graph_file = optarg;
    
    ifstream graph_in(graph_file);
    if (!graph_in) {
        cerr << "error: failed to open " << graph_file << '\n';
        return 1;
    }
    
    BaseGraph graph = read_gfa(graph_in, false);
    
    vector<string> retained_samples;
    BaseGraph pruned;
    for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
        pruned.add_node(graph.label(node_id));
    }
    for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
        for (uint64_t next_id : graph.next(node_id)) {
            pruned.add_edge(node_id, next_id);
        }
    }
    for (uint64_t path_id = 0; path_id < graph.path_size(); ++path_id) {
        
        string path_name = graph.path_name(path_id);
        if (removed_samples.count(path_name)) {
            if (!fasta_prefix.empty()) {
                
                string fasta_file = fasta_prefix + "_" + path_name + ".fasta";
                ofstream fasta_out(fasta_file);
                if (!fasta_out) {
                    cerr << "error: failed to write to " << fasta_file << '\n';
                    return 1;
                }
                fasta_out << '>' << path_name << '\n';
                const auto& path = graph.path(path_id);
                for (size_t i = 0; i < path.size(); ++i) {
                    fasta_out << graph.label(path[i]);
                    if (i % 80 == (80 - 1)) {
                        fasta_out << '\n';
                    }
                }
            }
        }
        else {
            uint64_t new_path_id = pruned.add_path(path_name);
            for (auto node_id : graph.path(path_id)) {
                pruned.extend_path(new_path_id, node_id);
            }
            retained_samples.emplace_back(std::move(path_name));
        }
    }
    
    auto tableau = add_sentinels(pruned, '^', '$');
    
    purge_uncovered_nodes(pruned, tableau);
    
    ExposedCore graph_namer(retained_samples, graph_prefix);
    
    string graph_out_file = graph_namer.root_subproblem_name();
    
    ofstream graph_out(graph_out_file);
    if (!graph_out) {
        cerr << "error: failed to write to " << graph_out_file << '\n';
        return 1;
    }
    
    write_gfa(pruned, tableau, graph_out, false);
    
    if (!tree_out_file.empty()) {
        
        ifstream tree_in(tree_in_file);
        if (!tree_in) {
            cerr << "error: failed to read from " << tree_in_file << '\n';
            return 1;
        }
        
        stringstream sstrm;
        sstrm << tree_in.rdbuf();
        string newick_string = std::move(sstrm.str());
        
        Tree tree(newick_string);
        
        vector<uint64_t> retained_leaf_ids;
        for (uint64_t node_id = 0; node_id < tree.node_size(); ++node_id) {
            if (tree.is_leaf(node_id) && !removed_samples.count(tree.label(node_id))) {
                retained_leaf_ids.push_back(node_id);
            }
        }
        
        tree.prune(retained_leaf_ids);
        tree.compact();
        
        string retained_newick_string = tree.to_newick();
        
        size_t num_parens = removed_samples.size();
        if (tree.node_size() == 0 && retained_samples.size() == 1) {
            num_parens = 0;
        }
        
        vector<string> ordered_removed_samples(removed_samples.begin(), removed_samples.end());
        sort(ordered_removed_samples.begin(), ordered_removed_samples.end());
        
        stringstream regrafted_newick_stream;
        for (size_t i = 0; i < num_parens; ++i) {
            regrafted_newick_stream << '(';
        }
        for (size_t i = 0; retained_newick_string[i] != ';'; ++i) {
            regrafted_newick_stream << retained_newick_string[i];
        }
        for (size_t i = 0; i < removed_samples.size(); ++i) {
            if (i < num_parens) {
                regrafted_newick_stream << ',';
            }
            regrafted_newick_stream << ordered_removed_samples[i] << ":0";
            if (i < num_parens) {
                regrafted_newick_stream << ')';
            }
        }
        regrafted_newick_stream << ";\n";
        
        ofstream tree_out(tree_out_file);
        if (!tree_out) {
            cerr << "error: failed to write to " << tree_out_file << '\n';
        }
        tree_out << regrafted_newick_stream.str();
    }
    
    return 0;
}
