#include <iostream>
#include <fstream>
#include <cassert>
#include <getopt.h>

#include "centrolign/gfa.hpp"
#include "centrolign/bridges.hpp"
#include "centrolign/modify_graph.hpp"
#include "centrolign/compacted_graph.hpp"
#include "centrolign/adjacency_graph.hpp"
#include "centrolign/step_index.hpp"
#include "centrolign/utility.hpp"

using namespace std;
using namespace centrolign;

void print_help() {
    cerr << "usage:\n";
    cerr << "find_universal_bridges [options] graph.gfa > bridge.txt\n\n";
    cerr << "options:\n";
    cerr << " --rightmost / -r  Return the rightmost rather than leftmost universal bridge\n";
    cerr << " --help / -h       Print this message and exit\n";
}

int main(int argc, char* argv[]) {
    
    bool leftmost = true;
    
    while (true)
    {
        static struct option options[] = {
            {"rightmost", no_argument, NULL, 'r'},
            {"help", no_argument, NULL, 'h'},
            {NULL, 0, NULL, 0}
        };
        int o = getopt_long(argc, argv, "rh", options, NULL);
        
        if (o == -1) {
            // end of uptions
            break;
        }
        switch (o)
        {
            case 'r':
                leftmost = false;
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
        std::cerr << "error: expected 1 positional argument but got " << (argc - optind) << "\n";
        print_help();
        return 1;
    }
    
    ifstream gfa_in(argv[optind]);
    if (!gfa_in) {
        std::cerr << "error: could not open GFA file " << argv[optind] << '\n';
        return 1;
    }
    
    gfa_in.seekg(0);
    if (!gfa_in) {
        std::cerr << "error: GFA file " << argv[optind] << " is not seekable\n";
        return 1;
    }
    
    std::cerr << "Loading graph...\n";
    BaseGraph graph = read_gfa(gfa_in, false);
    
    if (graph.node_size() == 0) {
        std::cerr << "error: Graph is empty.\n";
        return 1;
    }
    if (graph.path_size() == 0) {
        std::cerr << "error: Graph has no paths.\n";
        return 1;
    }

    std::cerr << "Finding bridges...\n";
    auto tableau = add_sentinels(graph, '^', '$');
    
    CompactedGraph compacted_graph(graph);
    AdjacencyGraph adj_graph(compacted_graph);
    
    auto adj_bridges = bridges(adj_graph);
    
    
    std::cerr << "Choosing between bridges...\n";
    
    StepIndex step_index(graph);
    
    uint64_t bridge_node_id = -1;
    size_t bridge_position = -1;
    for (const auto& bridge : adj_bridges) {
        uint64_t compacted_id = -1;
        for (const auto& edge : adj_graph.next_edges(bridge.first)) {
            if (edge.target == bridge.second) {
                compacted_id = edge.label;
                break;
            }
        }
        assert(compacted_id != -1);
        uint64_t node_id = leftmost ? compacted_graph.front(compacted_id) : compacted_graph.back(compacted_id);
        if (node_id == tableau.src_id && leftmost) {
            if (compacted_graph.back(compacted_id) != tableau.src_id) {
                node_id = graph.next(node_id).front();
            }
            else {
                continue;
            }
        }
        if (node_id == tableau.snk_id && !leftmost) {
            if (compacted_graph.front(compacted_id) != tableau.snk_id) {
                node_id = graph.previous(node_id).front();
            }
            else {
                continue;
            }
        }
        
        for (const auto& path_step : step_index.path_steps(node_id)) {
            if (path_step.first == 0) {
                if (bridge_position == -1 ||
                    (leftmost && bridge_position > path_step.second) ||
                    (!leftmost && bridge_position < path_step.second)) {
                    bridge_node_id = node_id;
                    bridge_position = path_step.second;
                }
                break;
            }
        }
    }
    
    if (bridge_node_id == -1) {
        std::cout << "Graph does not contain any bridges.\n";
    }
    else {
        gfa_in.clear();
        gfa_in.seekg(0);
        if (!gfa_in) {
            std::cerr << "error: Seek to beginning of GFA failed.\n";
            return 1;
        }

        std::vector<size_t> node_len;
        std::string line;
        size_t path_stream_pos = -1;
        while (gfa_in) {
            getline(gfa_in, line);

            if (line.empty()) {
                continue;
            }
            if (line.front() == 'S') {
                auto tokens = tokenize(line);
                assert(tokens.size() == 3);
                int64_t node_id = parse_int(tokens[1]);
                while (node_len.size() <= node_id) {
                    node_len.push_back(-1);
                }
                node_len[node_id] = tokens[2].size();
            }
            else if (line.front() == 'P') {
                std::string path_name = line.substr(2, line.find('\t', 2) - 2);
                if (path_name == graph.path_name(0)) {
                    size_t line_end = gfa_in.tellg();
                    path_stream_pos = line_end  - line.size() - 1;
                }
            }
        }

        assert(path_stream_pos != -1);
        
        gfa_in.clear();
        gfa_in.seekg(path_stream_pos);
        getline(gfa_in, line);
        auto steps = tokenize(line.substr(line.find('\t', 2) + 1), ',');
        size_t walked_pos = 0;
        size_t i = 0;
        while (walked_pos < bridge_position) {
            const auto& step = steps[i++];
            walked_pos += node_len[parse_int(step.substr(0, step.size() - 1))];
        }

        std::cout << "node:\t" << steps[i - 1].substr(0, steps[i - 1].size() - 1) << '\n';
        for (const auto& path_step : step_index.path_steps(bridge_node_id)) {
            std::cout << "position:\t" << graph.path_name(path_step.first) << '\t' << path_step.second << '\n';
        }
    }
    
    return 0;
}
