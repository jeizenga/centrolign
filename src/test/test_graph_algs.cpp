#include <cstdio>
#include <cstdlib>
#include <cstdlib>
#include <unordered_set>

#include "centrolign/utility.hpp"
#include "centrolign/random_graph.hpp"
#include "centrolign/determinize.hpp"
#include "centrolign/topological_order.hpp"

using namespace std;
using namespace centrolign;

bool is_reverse_deterministic(const BaseGraph& graph) {
    for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
        unordered_set<char> prev_bases;
        for (uint64_t prev_id : graph.previous(node_id)) {
            char base = graph.base(prev_id);
            if (prev_bases.count(base)) {
                return false;
            }
            prev_bases.insert(base);
        }
    }
    return true;
}

bool is_probably_equivalent(const BaseGraph& graph,
                            const BaseGraph& determinized,
                            default_random_engine& gen) {
    
    auto select_random_path = [&](const BaseGraph& graph) {
        
        std::vector<uint64_t> sources;
        for (uint64_t nid = 0; nid < graph.node_size(); ++nid) {
            if (graph.previous_size(nid) == 0) {
                sources.push_back(nid);
            }
        }
        
        uint64_t here = sources[uniform_int_distribution<size_t>(0, sources.size() - 1)(gen)];
        std::string sequence(1, graph.base(here));
        while (graph.next_size(here) != 0) {
            vector<uint64_t> nexts = graph.next(here);
            size_t i = uniform_int_distribution<size_t>(0, nexts.size() - 1)(gen);
            here = nexts[i];
            sequence.push_back(graph.base(here));
        }
        return sequence;
    };
    
    
    
    auto find_path = [&](const BaseGraph& graph, string& seq) {
        std::vector<uint64_t> sources;
        for (uint64_t nid = 0; nid < graph.node_size(); ++nid) {
            if (graph.previous_size(nid) == 0 && graph.base(nid) == seq.front()) {
                sources.push_back(nid);
            }
        }
        
        bool found_path = false;
        for (auto start_id : sources) {
            
            vector<pair<uint64_t, size_t>> stack(1, make_pair(start_id, (size_t) 0));
            while (!stack.empty()) {
                size_t i;
                uint64_t nid;
                tie(nid, i) = stack.back();
                stack.pop_back();
                if (i + 1 == seq.size()) {
                    if (graph.next_size(nid) == 0) {
                        found_path = true;
                        break;
                    }
                }
                else {
                    for (uint64_t next_id : graph.next(nid)) {
                        if (graph.base(next_id) == seq[i + 1]) {
                            stack.emplace_back(next_id, i + 1);
                        }
                    }
                }
            }
            if (found_path) {
                break;
            }
        }
        return found_path;
    };
    
    size_t num_tests = 5;
    for (size_t i = 0; i < num_tests; ++i) {
        string seq = select_random_path(graph);
        if (!find_path(determinized, seq)) {
            cerr << "failed to find sequence " << seq << " from original graph:\n";
            print_graph(graph, cerr);
            cerr << "in determinized graph:\n";
            print_graph(determinized, cerr);
            exit(1);
        }
    }
    
    return true;
}


void test_topological_order(const BaseGraph& graph) {
    auto order = topological_order(graph);
    vector<size_t> index(order.size());
    for (size_t i = 0; i < order.size(); ++i) {
        index[order[i]] = i;
    }
    for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
        for (auto next_id : graph.next(node_id)) {
            if (index[node_id] >= index[next_id]) {
                std::cerr << "topological ordering failed on edge " << node_id << " -> " << next_id << " on graph:\n";
                print_graph(graph, cerr);
                exit(1);
            }
        }
        for (auto prev_id : graph.previous(node_id)) {
            if (index[prev_id] >= index[node_id]) {
                std::cerr << "topological ordering failed on edge " << prev_id << " -> " << node_id << " on graph:\n";
                print_graph(graph, cerr);
                exit(1);
            }
        }
    }
}

void do_tests(const BaseGraph& graph, default_random_engine& gen) {
    BaseGraph determinized = determinize(graph);
    if (!is_reverse_deterministic(determinized)) {
        cerr << "determinization failed on graph:\n";
        print_graph(graph, cerr);
        cerr << "resulted in determinized graph:\n";
        print_graph(determinized, cerr);
        exit(1);
    }
    is_probably_equivalent(graph, determinized, gen);
    test_topological_order(graph);
    test_topological_order(determinized);
}

void add_sentinels(BaseGraph& graph) {
    assert(graph.node_size() != 0);
    uint64_t src_id = graph.add_node('#');
    uint64_t snk_id = graph.add_node('$');
    for (uint64_t node_id = 0; node_id + 2 < graph.node_size(); ++node_id) {
        if (graph.previous_size(node_id) == 0) {
            graph.add_edge(src_id, node_id);
        }
        if (graph.next_size(node_id) == 0) {
            graph.add_edge(node_id, snk_id);
        }
    }
}

int main(int argc, char* argv[]) {
     
    random_device rd;
    default_random_engine gen(rd());
    
    BaseGraph graph1;
    string graph1_labels = "GGGAT";
    for (size_t i = 0; i < 5; ++i) {
        graph1.add_node(graph1_labels[i]);
        for (size_t j = 0; j < i; ++j) {
            graph1.add_edge(j, i);
        }
    }
    // we have to have sentinels to make the determinize algorithm identify
    // the initial position of a sequenc
    add_sentinels(graph1);
    do_tests(graph1, gen);
    
    size_t num_reps = 10;
    vector<pair<size_t, size_t>> graph_sizes;
    graph_sizes.emplace_back(8, 15);
    graph_sizes.emplace_back(10, 12);
    graph_sizes.emplace_back(20, 35);
    for (auto& sizes : graph_sizes) {
        size_t num_nodes = sizes.first;
        size_t num_edges = sizes.second;
        for (size_t i = 0; i < num_reps; ++i) {
            BaseGraph graph = random_graph(num_nodes, num_edges, gen);
            add_sentinels(graph);
            do_tests(graph, gen);
        }
    }
    
    cerr << "passed all tests!" << endl;
}
