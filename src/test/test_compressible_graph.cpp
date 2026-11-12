#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <algorithm>
#include <stdexcept>

#include "centrolign/compressible_graph.hpp"
#include "centrolign/modify_graph.hpp"
#include "centrolign/test_util.hpp"

using namespace std;
using namespace centrolign;

void fail(const string& msg, const BaseGraph& graph) {
    cerr << "failed: " << msg << '\n';
    cerr << cpp_representation(graph, "graph") << '\n';
    exit(1);
}

// check that the decompressed graph is the same graph under the translation of node IDs
void check_translated(const BaseGraph& graph, const SentinelTableau& tableau,
                      const BaseGraph& decompressed, const SentinelTableau& decompressed_tableau,
                      const vector<uint64_t>& trans) {

    if (decompressed.node_size() != graph.node_size() || trans.size() != graph.node_size()) {
        fail("node sizes do not match", graph);
    }
    vector<bool> hit(graph.node_size(), false);
    for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
        if (trans[node_id] >= graph.node_size() || hit[trans[node_id]]) {
            fail("translation is not a permutation", graph);
        }
        hit[trans[node_id]] = true;
    }
    for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
        uint64_t new_id = trans[node_id];
        if (graph.label(node_id) != decompressed.label(new_id)) {
            fail("labels do not match on node " + to_string(node_id), graph);
        }
        // next edges must match in order
        vector<uint64_t> next_ids;
        for (auto next_id : graph.next(node_id)) {
            next_ids.push_back(trans[next_id]);
        }
        if (next_ids != decompressed.next(new_id)) {
            fail("next edges do not match on node " + to_string(node_id), graph);
        }
        // previous edges must match as multisets
        vector<uint64_t> prev_ids;
        for (auto prev_id : graph.previous(node_id)) {
            prev_ids.push_back(trans[prev_id]);
        }
        auto decompressed_prev_ids = decompressed.previous(new_id);
        sort(prev_ids.begin(), prev_ids.end());
        sort(decompressed_prev_ids.begin(), decompressed_prev_ids.end());
        if (prev_ids != decompressed_prev_ids) {
            fail("previous edges do not match on node " + to_string(node_id), graph);
        }
    }
    if (graph.path_size() != decompressed.path_size()) {
        fail("path sizes do not match", graph);
    }
    for (uint64_t path_id = 0; path_id < graph.path_size(); ++path_id) {
        if (graph.path_name(path_id) != decompressed.path_name(path_id) ||
            decompressed.path_id(graph.path_name(path_id)) != path_id) {
            fail("path names do not match", graph);
        }
        vector<uint64_t> path;
        for (auto node_id : graph.path(path_id)) {
            path.push_back(trans[node_id]);
        }
        if (path != decompressed.path(path_id)) {
            fail("paths do not match on path " + to_string(path_id), graph);
        }
    }
    if (tableau.src_id < graph.node_size()) {
        if (decompressed_tableau.src_id != trans[tableau.src_id] ||
            decompressed_tableau.snk_id != trans[tableau.snk_id] ||
            decompressed_tableau.src_sentinel != tableau.src_sentinel ||
            decompressed_tableau.snk_sentinel != tableau.snk_sentinel) {
            fail("tableaus do not match", graph);
        }
    }
}

void check_identical(const BaseGraph& graph1, const BaseGraph& graph2) {
    if (graph1.node_size() != graph2.node_size() || graph1.path_size() != graph2.path_size()) {
        fail("graphs are not identical sizes", graph1);
    }
    for (uint64_t node_id = 0; node_id < graph1.node_size(); ++node_id) {
        if (graph1.label(node_id) != graph2.label(node_id) ||
            graph1.next(node_id) != graph2.next(node_id) ||
            graph1.previous(node_id) != graph2.previous(node_id)) {
            fail("graphs are not identical at node " + to_string(node_id), graph1);
        }
    }
    for (uint64_t path_id = 0; path_id < graph1.path_size(); ++path_id) {
        if (graph1.path_name(path_id) != graph2.path_name(path_id) ||
            graph1.path(path_id) != graph2.path(path_id)) {
            fail("graphs are not identical at path " + to_string(path_id), graph1);
        }
    }
}

void do_test(const BaseGraph& graph, const SentinelTableau& tableau) {

    vector<uint64_t> trans;
    PackedSequenceGraph packed(graph, tableau, &trans);
    SentinelTableau decompressed_tableau = tableau;
    BaseGraph decompressed = packed.decompress(decompressed_tableau);
    check_translated(graph, tableau, decompressed, decompressed_tableau, trans);

    // the decompressed representation should be a fixed point
    vector<uint64_t> retrans;
    PackedSequenceGraph repacked(decompressed, decompressed_tableau, &retrans);
    SentinelTableau redecompressed_tableau = decompressed_tableau;
    BaseGraph redecompressed = repacked.decompress(redecompressed_tableau);
    for (uint64_t node_id = 0; node_id < retrans.size(); ++node_id) {
        if (retrans[node_id] != node_id) {
            fail("decompressed graph is not a fixed point of compression", graph);
        }
    }
    check_identical(decompressed, redecompressed);

    // and the wrapper should give the same result
    CompressibleGraph compressible(BaseGraph(graph), tableau);
    compressible.compress();
    compressible.compress();
    if (!compressible.is_compressed() || compressible.node_size() != graph.node_size()) {
        fail("wrapper did not compress", graph);
    }
    CompressibleGraph copied(compressible);
    CompressibleGraph moved(move(copied));
    check_identical(moved.get_base_graph(), decompressed);
    check_identical(compressible.get_base_graph(), decompressed);
    if (compressible.is_compressed() || compressible.get_tableau().src_id != decompressed_tableau.src_id) {
        fail("wrapper did not decompress", graph);
    }
}

template<class Generator>
void add_random_walks(BaseGraph& graph, size_t num_walks, size_t max_length, Generator& gen) {
    for (size_t i = 0; i < num_walks && graph.node_size() != 0; ++i) {
        uint64_t path_id = graph.add_path("walk" + to_string(graph.path_size()));
        uint64_t node_id = gen() % graph.node_size();
        graph.extend_path(path_id, node_id);
        size_t length = gen() % max_length;
        for (size_t j = 0; j < length && graph.next_size(node_id) != 0; ++j) {
            node_id = graph.next(node_id)[gen() % graph.next_size(node_id)];
            graph.extend_path(path_id, node_id);
        }
    }
}

int main(int argc, char* argv[]) {

    mt19937 gen(9824);

    // empty graph
    do_test(BaseGraph(), SentinelTableau());

    // single sequences with sentinels
    for (size_t len : {1, 2, 10, 1000}) {
        BaseGraph graph = make_base_graph("seq", random_sequence(len, gen));
        auto tableau = add_sentinels(graph, 5, 6);
        do_test(graph, tableau);
        reassign_sentinels(graph, tableau, 7, 8);
        do_test(graph, tableau);
    }

    // non-branching cycles, a self loop, and duplicate edges
    {
        BaseGraph graph;
        for (size_t i = 0; i < 10; ++i) {
            graph.add_node('A' + i);
        }
        graph.add_edge(0, 1);
        graph.add_edge(1, 2);
        graph.add_edge(2, 0);
        graph.add_edge(3, 3);
        graph.add_edge(4, 5);
        graph.add_edge(4, 5);
        graph.add_edge(9, 8);
        graph.add_edge(8, 7);
        graph.add_edge(7, 6);
        graph.add_edge(6, 9);
        uint64_t path_id = graph.add_path("cycle");
        for (uint64_t node_id : {1, 2, 0, 1, 2, 0}) {
            graph.extend_path(path_id, node_id);
        }
        path_id = graph.add_path("loop");
        for (uint64_t node_id : {3, 3, 3}) {
            graph.extend_path(path_id, node_id);
        }
        path_id = graph.add_path("empty");
        do_test(graph, SentinelTableau());

        // a path that doesn't follow edges can't be compressed
        path_id = graph.add_path("jumps");
        for (uint64_t node_id : {4, 5, 9}) {
            graph.extend_path(path_id, node_id);
        }
        bool threw = false;
        try {
            PackedSequenceGraph packed(graph, SentinelTableau());
        }
        catch (runtime_error& err) {
            threw = true;
        }
        if (!threw) {
            fail("compressed a path that doesn't follow edges", graph);
        }
    }

    // random graphs
    for (size_t rep = 0; rep < 200; ++rep) {
        size_t num_nodes = 1 + gen() % 60;
        size_t num_edges = gen() % (2 * num_nodes);
        bool acyclic = (rep % 2 == 0);
        BaseGraph graph = random_graph(num_nodes, num_edges, acyclic, gen);
        SentinelTableau tableau;
        if (acyclic) {
            add_random_path_cover(graph, gen);
            if (rep % 4 == 0) {
                tableau = add_sentinels(graph, 5, 6);
            }
        }
        add_random_walks(graph, 3, 20, gen);
        do_test(graph, tableau);
    }

    // challenge graphs with many bubbles
    for (size_t rep = 0; rep < 50; ++rep) {
        BaseGraph graph = random_challenge_graph(1 + gen() % 100, gen);
        add_random_path_cover(graph, gen);
        auto tableau = add_sentinels(graph, 5, 6);
        do_test(graph, tableau);
    }

    cerr << "passed all tests!" << endl;
    return 0;
}
