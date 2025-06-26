#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <random>

#include "centrolign/anchorer.hpp"
#include "centrolign/chain_merge.hpp"
#include "centrolign/path_merge.hpp"
#include "centrolign/modify_graph.hpp"
#include "centrolign/graph.hpp"
#include "centrolign/utility.hpp"
#include "centrolign/test_util.hpp"
#include "centrolign/match_finder.hpp"
#include "centrolign/packed_forward_edges.hpp"
#include "centrolign/forward_edges.hpp"

using namespace std;
using namespace centrolign;


class TestAnchorer : public Anchorer {
public:
    TestAnchorer(const ScoreFunction& score_function) : Anchorer(score_function) {}
    
    using Anchorer::AnchorGraph;
    using Anchorer::exhaustive_chain_dp;
    using Anchorer::sparse_chain_dp;
    using Anchorer::anchor_weight;
    using Anchorer::sparse_affine_chain_dp;
    using Anchorer::gap_open;
    using Anchorer::gap_extend;
    using Anchorer::edge_weight;
    using Anchorer::extract_graphs_between;
    using Anchorer::project_paths;
    using Anchorer::divvy_matches;
    using Anchorer::score_function;
};

// DFS to walk out all fixed-length paths
unordered_map<string, vector<vector<uint64_t>>> k_mer_walks(const BaseGraph& graph, size_t k) {
    
    assert(k != 0);
    
    unordered_map<string, vector<vector<uint64_t>>> to_return;
    
    for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
        
        // records of (node id, next edge)
        vector<pair<uint64_t, size_t>> stack;
        stack.emplace_back(node_id, 0);
        
        while (!stack.empty()) {
            
            if (stack.size() == k) {
                
                string seq;
                vector<uint64_t> walk;
                for (auto& rec : stack) {
                    seq.push_back(graph.label(rec.first));
                    walk.push_back(rec.first);
                }
                to_return[seq].push_back(walk);
                
                stack.pop_back();
            }
            else if (stack.back().second == graph.next_size(stack.back().first)) {
                stack.pop_back();
            }
            else {
                uint64_t next = graph.next(stack.back().first)[stack.back().second++];
                stack.emplace_back(next, 0);
            }
        }
        
    }
    
    return to_return;
}

std::vector<match_set_t> generate_anchor_set(const BaseGraph& graph1,
                                             const BaseGraph& graph2,
                                             size_t k) {
    
    auto kmers1 = k_mer_walks(graph1, k);
    auto kmers2 = k_mer_walks(graph2, k);
    
    std::vector<match_set_t> anchors;
    for (pair<const string, vector<vector<uint64_t>>>& entry : kmers1) {
        if (kmers2.count(entry.first)) {
            anchors.emplace_back();
            anchors.back().walks1 = std::move(entry.second);
            anchors.back().walks2 = std::move(kmers2[entry.first]);
            anchors.back().count1 = anchors.back().walks1.size();
            anchors.back().count2 = anchors.back().walks2.size();
            anchors.back().full_length = k;
        }
    }
    
    return anchors;
}

void print_anchor_set(const vector<match_set_t>& anchors, size_t i) {
    cerr << i << "\twalks on graph1:\n";
    for (const auto& walk : anchors[i].walks1) {
        cerr << "\t\t";
        for (size_t j = 0; j < walk.size(); ++j) {
            if (j) {
                cerr << ',';
            }
            cerr << walk[j];
        }
        cerr << '\n';
    }
    cerr << "\twalks on graph2:\n";
    for (const auto& walk : anchors[i].walks2) {
        cerr << "\t\t";
        for (size_t j = 0; j < walk.size(); ++j) {
            if (j) {
                cerr << ',';
            }
            cerr << walk[j];
        }
        cerr << '\n';
    }
}

void print_chain(const TestAnchorer& anchorer, const vector<anchor_t>& chain,
                 const vector<double>& gap_scores) {
    
    // non-affine, local, global
    assert(chain.size() == 0 || chain.size() == gap_scores.size() + 1 || gap_scores.size() == chain.size() + 1);
    
    if (gap_scores.size() > chain.size()) {
        std::cerr << "(initial gap " <<  gap_scores[0] << ")\n";
    }
    for (size_t i = 0; i < chain.size(); ++i) {
        auto& link = chain[i];
        cerr << i << " (count1 " << link.count1 << ", count2 " << link.count2 << ", score " << anchorer.score_function->anchor_weight(link.count1, link.count2, link.walk1.size(), link.full_length) << ", set " << link.match_set << ", idx1 " << link.idx1 << ", idx2 " << link.idx2 << "):\n";
        for (auto walk : {link.walk1, link.walk2}) {
            cerr << '\t';
            for (size_t j = 0; j < walk.size(); ++j) {
                if (j) {
                    cerr << ',';
                }
                cerr << walk[j];
            }
            cerr << '\n';
        }
        if (!gap_scores.empty() && i + 1 != chain.size()) {
            cerr << "(edge weight " << gap_scores[chain.size() > gap_scores.size() ? i : i + 1] << ")\n";
        }
    }
    if (gap_scores.size() > chain.size() && gap_scores.size() > 1) {
        std::cerr << "(final gap " <<  gap_scores.back() << ")\n";
    }
}

void test_sparse_dynamic_programming(const BaseGraph& graph1,
                                     const BaseGraph& graph2,
                                     const std::vector<match_set_t>& anchors,
                                     uint64_t src1, uint64_t snk1,
                                     uint64_t src2, uint64_t snk2,
                                     bool affine, bool global, bool packed) {
    
    PathMerge<> chain_merge1(graph1);
    PathMerge<> chain_merge2(graph2);
    
    ScoreFunction score_function;
    score_function.anchor_score_function = ScoreFunction::InverseCount;
    TestAnchorer anchorer(score_function);
    anchorer.gap_open[0] = 0.25;
    anchorer.gap_extend[0] = 0.25;
    anchorer.gap_open[1] = .4;
    anchorer.gap_extend[1] = 0.15;
    
    vector<uint64_t> sources1(1, src1);
    vector<uint64_t> sources2(1, src2);
    vector<uint64_t> sinks1(1, snk1);
    vector<uint64_t> sinks2(1, snk2);
    
    std::vector<anchor_t> exhaustive_chain, sparse_chain;
    {
        auto anchors_copy = anchors;
        if (affine) {
            if (global) {
                if (packed) {
                    sparse_chain = anchorer.sparse_affine_chain_dp<size_t, size_t, size_t, int64_t, size_t, float,
                                                                   VectorPair<SignedPackedVector, PackedVector, int64_t, PackedMatchBank<size_t>::match_id_t>,
                                                                   VectorPair<PackedVector, PackedVector, size_t, PackedMatchBank<size_t>::match_id_t>,
                                                                   PackedVector, PackedVector, PackedMatchBank<size_t>, PackedForwardEdges>
                                                                  (anchors_copy,
                                                                   graph1,
                                                                   graph2,
                                                                   chain_merge1,
                                                                   chain_merge2,
                                                                   anchorer.gap_open,
                                                                   anchorer.gap_extend, 1.0,
                                                                   anchors_copy.size(), true,
                                                                   &sources1, &sources2, &sinks1, &sinks2);
                }
                else {
                    sparse_chain = anchorer.sparse_affine_chain_dp<size_t, size_t, size_t, int64_t, size_t, float,
                                                                   std::vector<std::pair<int64_t, MatchBank<size_t, size_t>::match_id_t>>,
                                                                   std::vector<std::pair<size_t, MatchBank<size_t, size_t>::match_id_t>>,
                                                                   std::vector<size_t>, std::vector<size_t>, MatchBank<size_t, size_t>, ForwardEdges<uint64_t, uint64_t>>
                                                                  (anchors_copy,
                                                                   graph1,
                                                                   graph2,
                                                                   chain_merge1,
                                                                   chain_merge2,
                                                                   anchorer.gap_open,
                                                                   anchorer.gap_extend, 1.0,
                                                                   anchors_copy.size(), true,
                                                                   &sources1, &sources2, &sinks1, &sinks2);
                }
            }
            else {
                if (packed) {
                    
                    sparse_chain = anchorer.sparse_affine_chain_dp<size_t, size_t, size_t, int64_t, size_t, float,
                                                                  VectorPair<SignedPackedVector, PackedVector, int64_t, PackedMatchBank<size_t>::match_id_t>,
                                                                  VectorPair<PackedVector, PackedVector, size_t, PackedMatchBank<size_t>::match_id_t>,
                                                                  PackedVector, PackedVector, PackedMatchBank<size_t>, PackedForwardEdges>
                                                                 (anchors_copy,
                                                                  graph1,
                                                                  graph2,
                                                                  chain_merge1,
                                                                  chain_merge2,
                                                                  anchorer.gap_open,
                                                                  anchorer.gap_extend, 1.0,
                                                                  anchors_copy.size(), true);
                }
                else {
                    
                    sparse_chain = anchorer.sparse_affine_chain_dp<size_t, size_t, size_t, int64_t, size_t, float,
                                                                   std::vector<std::pair<int64_t, MatchBank<size_t, size_t>::match_id_t>>,
                                                                   std::vector<std::pair<size_t, MatchBank<size_t, size_t>::match_id_t>>,
                                                                   std::vector<size_t>, std::vector<size_t>, MatchBank<size_t, size_t>, ForwardEdges<uint64_t, uint64_t>>
                                                                  (anchors_copy,
                                                                   graph1,
                                                                   graph2,
                                                                   chain_merge1,
                                                                   chain_merge2,
                                                                   anchorer.gap_open,
                                                                   anchorer.gap_extend, 1.0,
                                                                   anchors_copy.size(), true);
                }
            }
        }
        else {
            if (global) {
                if (packed) {
                    sparse_chain = anchorer.sparse_chain_dp<size_t, size_t, size_t, size_t, float,
                                                            VectorPair<PackedVector, PackedVector, size_t, PackedMatchBank<size_t>::match_id_t>,
                                                            PackedVector, PackedMatchBank<size_t>, PackedForwardEdges>
                                                           (anchors_copy,
                                                            graph1,
                                                            chain_merge1,
                                                            chain_merge2, anchors_copy.size(), true,
                                                            &sources1, &sources2, &sinks1, &sinks2);
                }
                else {
                    
                    sparse_chain = anchorer.sparse_chain_dp<size_t, size_t, size_t, size_t, float,
                                                            std::vector<std::pair<size_t, MatchBank<size_t, size_t>::match_id_t>>,
                                                            std::vector<size_t>, MatchBank<size_t, size_t>, ForwardEdges<uint64_t, uint64_t>>
                                                           (anchors_copy,
                                                            graph1,
                                                            chain_merge1,
                                                            chain_merge2, anchors_copy.size(), true,
                                                            &sources1, &sources2, &sinks1, &sinks2);
                }
            }
            else {
                if (packed) {
                    sparse_chain = anchorer.sparse_chain_dp<size_t, size_t, size_t, size_t, float,
                                                            VectorPair<PackedVector, PackedVector, size_t, PackedMatchBank<size_t>::match_id_t>,
                                                            PackedVector, PackedMatchBank<size_t>, PackedForwardEdges>
                                                           (anchors_copy,
                                                            graph1,
                                                            chain_merge1,
                                                            chain_merge2, anchors_copy.size(), true);
                }
                else {
                    
                    sparse_chain = anchorer.sparse_chain_dp<size_t, size_t, size_t, size_t, float,
                                                            std::vector<std::pair<size_t, MatchBank<size_t, size_t>::match_id_t>>,
                                                            std::vector<size_t>, MatchBank<size_t, size_t>, ForwardEdges<uint64_t, uint64_t>>
                                                           (anchors_copy,
                                                            graph1,
                                                            chain_merge1,
                                                            chain_merge2, anchors_copy.size(), true);
                }
            }
        }
    }
    {
        auto anchors_copy = anchors;
        if (global) {
            exhaustive_chain = anchorer.exhaustive_chain_dp<size_t, size_t>(anchors_copy,
                                                            graph1, graph2,
                                                            chain_merge1,
                                                            chain_merge2,
                                                            affine, 1.0, anchors_copy.size(),
                                                            &sources1, &sources2, &sinks1, &sinks2);
        }
        else {
            exhaustive_chain = anchorer.exhaustive_chain_dp<size_t, size_t>(anchors_copy,
                                                            graph1, graph2,
                                                            chain_merge1,
                                                            chain_merge2,
                                                            affine, 1.0, anchors_copy.size());
        }
    }
    
    // score the anchors
    double exhaustive_score = 0.0, sparse_score = 0.0;
    for (auto& link : exhaustive_chain) {
        exhaustive_score += anchorer.score_function->anchor_weight(link.count1, link.count2, link.walk1.size(), link.full_length);
    }
    for (auto& link : sparse_chain) {
        sparse_score += anchorer.score_function->anchor_weight(link.count1, link.count2, link.walk1.size(), link.full_length);
    }
    
    std::vector<double> gap_costs_sparse, gap_costs_exhaustive;
    
    PostSwitchDistances<vector<size_t>> switch_dists1, switch_dists2;
    if (affine) {
        switch_dists1 = std::move(PostSwitchDistances<vector<size_t>>(graph1, chain_merge1));
        switch_dists2 = std::move(PostSwitchDistances<vector<size_t>>(graph2, chain_merge2));
        
        if (global) {
            // score up the first/last edge
            if (!exhaustive_chain.empty()) {
                gap_costs_exhaustive.push_back(anchorer.edge_weight(src1, exhaustive_chain.front().walk1.front(),
                                                                    src2, exhaustive_chain.front().walk2.front(), 1.0,
                                                                    chain_merge1, chain_merge2, switch_dists1, switch_dists2));
            }
            else {
                gap_costs_exhaustive.push_back(anchorer.edge_weight(src1, snk1,
                                                                    src2, snk2, 1.0,
                                                                    chain_merge1, chain_merge2, switch_dists1, switch_dists2));
            }
            
            if (!sparse_chain.empty()) {
                gap_costs_sparse.push_back(anchorer.edge_weight(src1, sparse_chain.front().walk1.front(),
                                                                src2, sparse_chain.front().walk2.front(), 1.0,
                                                                chain_merge1, chain_merge2, switch_dists1, switch_dists2));
            }
            else {
                gap_costs_sparse.push_back(anchorer.edge_weight(src1, snk1,
                                                                src2, snk2, 1.0,
                                                                chain_merge1, chain_merge2, switch_dists1, switch_dists2));
            }
        }
        
        // score up the edges
        for (size_t i = 1; i < exhaustive_chain.size(); ++i) {
            gap_costs_exhaustive.push_back(anchorer.edge_weight(exhaustive_chain[i-1].walk1.back(), exhaustive_chain[i].walk1.front(),
                                                                exhaustive_chain[i-1].walk2.back(), exhaustive_chain[i].walk2.front(), 1.0,
                                                                chain_merge1, chain_merge2, switch_dists1, switch_dists2));
        }
        for (size_t i = 1; i < sparse_chain.size(); ++i) {
            gap_costs_sparse.push_back(anchorer.edge_weight(sparse_chain[i-1].walk1.back(), sparse_chain[i].walk1.front(),
                                                            sparse_chain[i-1].walk2.back(), sparse_chain[i].walk2.front(), 1.0,
                                                            chain_merge1, chain_merge2, switch_dists1, switch_dists2));
        }
        
        if (global) {
            // score up the first/last edge
            if (!exhaustive_chain.empty()) {
                gap_costs_exhaustive.push_back(anchorer.edge_weight(exhaustive_chain.back().walk1.back(), snk1,
                                                                    exhaustive_chain.back().walk2.back(), snk2, 1.0,
                                                                    chain_merge1, chain_merge2, switch_dists1, switch_dists2));
            }

            if (!sparse_chain.empty()) {
                gap_costs_sparse.push_back(anchorer.edge_weight(sparse_chain.back().walk1.back(), snk1,
                                                                sparse_chain.back().walk2.back(), snk2, 1.0,
                                                                chain_merge1, chain_merge2, switch_dists1, switch_dists2));
            }
        }
    }
    
    for (auto c : gap_costs_exhaustive) {
        exhaustive_score += c;
    }
    for (auto c : gap_costs_sparse) {
        sparse_score += c;
    }
    
    if (abs(exhaustive_score - sparse_score) > 1e-6) {
        cerr << "did not find equivalent chains with sparse and exhaustive DP, affine? " << affine << ", global? " << global << ", packed? " << packed  << "\n";
        cerr << "boundaries: " << src1 << ":" << snk1 << ", " << src2 << ":" << snk2 << '\n';
        cerr << "anchor sets:\n";
        for (size_t i = 0; i < anchors.size(); ++i) {
            print_anchor_set(anchors, i);
        }
        cerr << "graphs:\n";
        cerr << cpp_representation(graph1, "graph1") << '\n';
        cerr << cpp_representation(graph2, "graph2") << '\n';
        cerr << "exhaustive chain (score " << exhaustive_score << "):\n";
        print_chain(anchorer, exhaustive_chain, gap_costs_exhaustive);
        cerr << "sparse chain (score " << sparse_score << "):\n";
        print_chain(anchorer, sparse_chain, gap_costs_sparse);
        exit(1);
    }
}


template<class Generator>
pair<uint64_t, uint64_t> generate_source_sink(const BaseGraph& graph, Generator& gen) {
    
    std::uniform_int_distribution<uint64_t> node_distr(0, graph.node_size() - 1);
    
    uint64_t src, snk;
    do {
        src = node_distr(gen);
        snk = node_distr(gen);
    } while (!is_reachable(graph, src, snk));
    
    return make_pair(src, snk);
}


int main(int argc, char* argv[]) {
     
    random_device rd;
    default_random_engine gen(rd());
    
    // heaviest weight path
    {
        TestAnchorer::AnchorGraph graph;
        uint64_t n0 = graph.add_node(0, 0, 0, 1.0);
        uint64_t n1 = graph.add_node(0, 0, 0, 0.5);
        uint64_t n2 = graph.add_node(0, 0, 0, 0.0);
        uint64_t n3 = graph.add_node(0, 0, 0, 0.5);
        uint64_t n4 = graph.add_node(0, 0, 0, 1.0);
        uint64_t n5 = graph.add_node(0, 0, 0, 1.0);
        uint64_t n6 = graph.add_node(0, 0, 0, 0.5);
        uint64_t n7 = graph.add_node(0, 0, 0, 1.0);

        graph.add_edge(n0, n2);
        graph.add_edge(n1, n2);
        graph.add_edge(n2, n3);
        graph.add_edge(n2, n4);
        graph.add_edge(n3, n5);
        graph.add_edge(n4, n5);
        graph.add_edge(n5, n6);
        graph.add_edge(n5, n7);

        vector<uint64_t> expected{n0, n2, n4, n5, n7};
        assert(graph.heaviest_weight_path() == expected);
    }
    // heaviest path with edge weights
    {
        TestAnchorer::AnchorGraph graph;
        uint64_t n0 = graph.add_node(0, 0, 0, 1.0);
        uint64_t n1 = graph.add_node(0, 0, 0, 1.5);
        uint64_t n2 = graph.add_node(0, 0, 0, 1.0);

        graph.add_edge(n0, n1, -1.0);
        graph.add_edge(n0, n2, 0.0);
        graph.add_edge(n1, n2, -1.0);

        vector<uint64_t> expected{n0, n2};
        assert(graph.heaviest_weight_path() == expected);
    }
    
    {
        BaseGraph graph1;
        for (auto c : std::string("ATGCATGTAA")) {
            graph1.add_node(c);
        }
        
        std::vector<std::pair<int, int>> graph1_edges{
            {0, 3},
            {0, 8},
            {0, 7},
            {0, 1},
            {1, 8},
            {2, 8},
            {2, 9},
            {2, 3},
            {2, 4},
            {3, 5},
            {3, 7},
            {3, 8},
            {4, 9},
            {4, 6},
            {4, 8},
            {5, 7},
            {6, 9},
            {8, 9}
        };
        
        std::vector<std::vector<int>> graph1_paths{
            {2, 3, 5, 7},
            {2, 4, 6, 9},
            {0, 8, 9},
            {0, 1, 8, 9}
        };
        
        for (auto e : graph1_edges) {
            graph1.add_edge(e.first, e.second);
        }
        
        for (size_t i = 0; i < graph1_paths.size(); ++i) {
            auto p = graph1.add_path(std::to_string(i));
            for (auto n : graph1_paths[i]) {
                graph1.extend_path(p, n);
            }
        }
        
        BaseGraph graph2;
        for (auto c : std::string("CACCACTTAT")) {
            graph2.add_node(c);
        }
        
        std::vector<std::pair<int, int>> graph2_edges{
            {0, 9},
            {0, 5},
            {0, 7},
            {1, 9},
            {1, 8},
            {1, 3},
            {1, 4},
            {2, 6},
            {2, 8},
            {2, 9},
            {3, 7},
            {4, 5},
            {4, 9},
            {5, 9},
            {6, 8},
            {6, 9},
            {7, 8},
            {7, 9}
        };
        
        std::vector<std::vector<int>> graph2_paths{
            {2, 6, 8},
            {0, 5, 9},
            {1, 4, 9},
            {1, 3, 7, 9}
        };
        
        for (auto e : graph2_edges) {
            graph2.add_edge(e.first, e.second);
        }
        
        for (size_t i = 0; i < graph2_paths.size(); ++i) {
            auto p = graph2.add_path(std::to_string(i));
            for (auto n : graph2_paths[i]) {
                graph2.extend_path(p, n);
            }
        }

        auto anchors = generate_anchor_set(graph1, graph2, 2);
        test_sparse_dynamic_programming(graph1, graph2, anchors, 2, 8, 3, 9, true, false, false);
        test_sparse_dynamic_programming(graph1, graph2, anchors, 2, 8, 3, 9, true, false, true);
    }

    {
        BaseGraph graph1;
        for (auto c : std::string("ATGGACCCTTGTGGCA")) {
            graph1.add_node(c);
        }
        
        std::vector<std::pair<int, int>> graph1_edges{
            {0, 15},
            {0, 3},
            {0, 9},
            {0, 11},
            {0, 2},
            {1, 8},
            {1, 6},
            {1, 9},
            {2, 6},
            {2, 11},
            {2, 4},
            {2, 15},
            {2, 3},
            {3, 10},
            {3, 13},
            {3, 7},
            {3, 8},
            {4, 11},
            {4, 6},
            {5, 10},
            {5, 6},
            {5, 14},
            {6, 9},
            {6, 11},
            {7, 13},
            {7, 10},
            {8, 11},
            {8, 13},
            {9, 13},
            {9, 10},
            {9, 12},
            {10, 11},
            {10, 13},
            {12, 15},
            {12, 14}
        };
        
        std::vector<std::vector<int>> graph1_paths{
            {1, 8, 11},
            {0, 2, 3, 13},
            {0, 2, 4, 6, 9, 12, 14},
            {5, 10, 13},
            {0, 2, 15},
            {0, 3, 7, 10, 13}
        };
        
        for (auto e : graph1_edges) {
            graph1.add_edge(e.first, e.second);
        }
        
        for (size_t i = 0; i < graph1_paths.size(); ++i) {
            auto p = graph1.add_path(std::to_string(i));
            for (auto n : graph1_paths[i]) {
                graph1.extend_path(p, n);
            }
        }
        
        BaseGraph graph2;
        for (auto c : std::string("CGGGGAACGCACCTCA")) {
            graph2.add_node(c);
        }
        
        std::vector<std::pair<int, int>> graph2_edges{
            {0, 2},
            {0, 13},
            {0, 14},
            {1, 15},
            {1, 10},
            {1, 12},
            {1, 4},
            {2, 13},
            {2, 7},
            {2, 6},
            {2, 3},
            {2, 4},
            {2, 5},
            {2, 12},
            {2, 9},
            {3, 6},
            {3, 15},
            {3, 13},
            {3, 12},
            {4, 15},
            {4, 6},
            {4, 11},
            {4, 10},
            {5, 8},
            {5, 6},
            {5, 14},
            {6, 10},
            {6, 8},
            {7, 11},
            {7, 10},
            {8, 12},
            {9, 12},
            {9, 10},
            {10, 15},
            {10, 11}
        };
        
        std::vector<std::vector<int>> graph2_paths{
            {1, 4, 15},
            {0, 2, 13},
            {0, 2, 7, 11},
            {0, 2, 5, 8, 12},
            {0, 14},
            {0, 2, 9, 10, 15},
            {0, 2, 3, 6, 10, 15}
        };
        
        for (auto e : graph2_edges) {
            graph2.add_edge(e.first, e.second);
        }
        
        for (size_t i = 0; i < graph2_paths.size(); ++i) {
            auto p = graph2.add_path(std::to_string(i));
            for (auto n : graph2_paths[i]) {
                graph2.extend_path(p, n);
            }
        }
        
        auto anchors = generate_anchor_set(graph1, graph2, 2);
        test_sparse_dynamic_programming(graph1, graph2, anchors, 4, 15, 3, 13, true, false, false);
    }

    {
        BaseGraph graph1;
        for (auto c : std::string("CAACCCAATCCAACCCAACCCCACCAACAG")) {
            graph1.add_node(c);
        }

        std::vector<std::pair<int, int>> graph1_edges{
            {0, 1},
            {1, 2},
            {2, 3},
            {2, 4},
            {2, 8},
            {3, 4},
            {4, 5},
            {5, 6},
            {6, 7},
            {6, 25},
            {7, 8},
            {8, 9},
            {8, 10},
            {9, 10},
            {10, 11},
            {11, 12},
            {12, 13},
            {12, 18},
            {13, 14},
            {14, 15},
            {14, 26},
            {15, 16},
            {16, 17},
            {16, 28},
            {16, 19},
            {17, 18},
            {18, 19},
            {19, 20},
            {19, 29},
            {20, 21},
            {21, 22},
            {22, 23},
            {22, 27},
            {23, 24},
            {25, 8},
            {26, 16},
            {27, 24},
            {28, 18},
            {29, 21}
        };

        std::vector<std::vector<int>> graph1_paths{
            {0, 1, 2, 4, 5, 6, 25, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 29, 21, 22, 23, 24},
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 26, 16, 28, 18, 19, 20, 21, 22, 27, 24}
        };

        for (auto e : graph1_edges) {
            graph1.add_edge(e.first, e.second);
        }

        for (size_t i = 0; i < graph1_paths.size(); ++i) {
            auto p = graph1.add_path(std::to_string(i));
            for (auto n : graph1_paths[i]) {
                graph1.extend_path(p, n);
            }
        }

        BaseGraph graph2;
        for (auto c : std::string("TTAAAAAAAAAAAAAAAAAAAAAAAATAAA")) {
            graph2.add_node(c);
        }

        std::vector<std::pair<int, int>> graph2_edges{
            {0, 1},
            {1, 2},
            {1, 26},
            {1, 3},
            {2, 3},
            {3, 4},
            {3, 27},
            {4, 5},
            {4, 25},
            {5, 6},
            {5, 12},
            {6, 7},
            {7, 8},
            {8, 9},
            {9, 10},
            {10, 11},
            {11, 12},
            {12, 13},
            {13, 14},
            {14, 15},
            {15, 16},
            {16, 17},
            {16, 29},
            {16, 18},
            {16, 22},
            {17, 18},
            {18, 19},
            {19, 20},
            {20, 21},
            {20, 22},
            {20, 23},
            {21, 22},
            {21, 23},
            {22, 23},
            {23, 24},
            {23, 28},
            {25, 6},
            {26, 3},
            {27, 5},
            {27, 25},
            {27, 8},
            {29, 18}
        };

        std::vector<std::vector<int>> graph2_paths{
            {0, 1, 26, 3, 4, 25, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24},
            {0, 1, 2, 3, 27, 5, 12, 13, 14, 15, 16, 29, 18, 19, 20, 21, 23, 28}
        };

        for (auto e : graph2_edges) {
            graph2.add_edge(e.first, e.second);
        }

        for (size_t i = 0; i < graph2_paths.size(); ++i) {
            auto p = graph2.add_path(std::to_string(i));
            for (auto n : graph2_paths[i]) {
                graph2.extend_path(p, n);
            }
        }

        auto anchors = generate_anchor_set(graph1, graph2, 2);
        test_sparse_dynamic_programming(graph1, graph2, anchors, 29, 23, 17, 19, true, true, false);
    }

    {
        BaseGraph graph1;
        for (auto c : std::string("AAAAAGCACA")) {
            graph1.add_node(c);
        }

        std::vector<std::pair<int, int>> graph1_edges{
            {0, 1},
            {1, 2},
            {1, 7},
            {1, 3},
            {1, 4},
            {2, 3},
            {3, 4},
            {3, 9},
            {3, 9},
            {5, 1},
            {6, 1},
            {7, 3},
            {8, 1}
        };

        std::vector<std::vector<int>> graph1_paths{
            {6, 1, 2, 3, 4},
            {8, 1, 7, 3, 9},
            {5, 1, 4},
            {0, 1, 4}
        };

        for (auto e : graph1_edges) {
            graph1.add_edge(e.first, e.second);
        }

        for (size_t i = 0; i < graph1_paths.size(); ++i) {
            auto p = graph1.add_path(std::to_string(i));
            for (auto n : graph1_paths[i]) {
                graph1.extend_path(p, n);
            }
        }

        BaseGraph graph2;
        for (auto c : std::string("AAACAATAAC")) {
            graph2.add_node(c);
        }

        std::vector<std::pair<int, int>> graph2_edges{
            {0, 1},
            {0, 9},
            {0, 1},
            {1, 2},
            {1, 7},
            {1, 8},
            {1, 3},
            {1, 6},
            {1, 4},
            {2, 3},
            {2, 6},
            {2, 4},
            {2, 3},
            {3, 4},
            {5, 1},
            {5, 9},
            {6, 4},
            {7, 3},
            {7, 6},
            {7, 4},
            {7, 4},
            {8, 3},
            {8, 6},
            {8, 4},
            {9, 2},
            {9, 7},
            {9, 8},
            {9, 3},
            {9, 6},
            {9, 4}
        };

        std::vector<std::vector<int>> graph2_paths{
            {0, 9, 4},
            {5, 1, 2, 3, 4},
            {0, 9, 7, 6, 4},
            {5, 9, 8, 6, 4}
        };

        for (auto e : graph2_edges) {
            graph2.add_edge(e.first, e.second);
        }

        for (size_t i = 0; i < graph2_paths.size(); ++i) {
            auto p = graph2.add_path(std::to_string(i));
            for (auto n : graph2_paths[i]) {
                graph2.extend_path(p, n);
            }
        }

        auto anchors = generate_anchor_set(graph1, graph2, 2);
        test_sparse_dynamic_programming(graph1, graph2, anchors, 0, 2, 5, 2, true, true, false);
    }

    {
        BaseGraph graph1;
        for (auto c : std::string("ACCACCACCACCACCCACCG")) {
            graph1.add_node(c);
        }

        std::vector<std::pair<int, int>> graph1_edges{
            {0, 1},
            {1, 2},
            {1, 15},
            {1, 3},
            {2, 3},
            {3, 4},
            {4, 5},
            {5, 6},
            {5, 19},
            {6, 7},
            {7, 8},
            {8, 9},
            {8, 18},
            {9, 10},
            {9, 17},
            {9, 11},
            {9, 14},
            {9, 12},
            {10, 11},
            {10, 11},
            {11, 12},
            {12, 13},
            {12, 18},
            {12, 14},
            {12, 16},
            {13, 14},
            {13, 16},
            {15, 3},
            {17, 11},
            {18, 14},
            {18, 16},
            {19, 7}
        };

        std::vector<std::vector<int>> graph1_paths{
            {0, 1, 15, 3, 4, 5, 6, 7, 8, 9, 14},
            {0, 1, 2, 3, 4, 5, 19, 7, 8, 9, 10, 11, 12, 18, 16},
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 17, 11, 12, 13, 14}
        };

        for (auto e : graph1_edges) {
            graph1.add_edge(e.first, e.second);
        }

        for (size_t i = 0; i < graph1_paths.size(); ++i) {
            auto p = graph1.add_path(std::to_string(i));
            for (auto n : graph1_paths[i]) {
                graph1.extend_path(p, n);
            }
        }

        BaseGraph graph2;
        for (auto c : std::string("CCCCCCCCCCCCCCCGGGCC")) {
            graph2.add_node(c);
        }

        std::vector<std::pair<int, int>> graph2_edges{
            {0, 1},
            {0, 2},
            {0, 18},
            {0, 5},
            {1, 2},
            {1, 18},
            {2, 3},
            {3, 4},
            {3, 5},
            {4, 5},
            {5, 6},
            {5, 7},
            {5, 19},
            {5, 19},
            {6, 7},
            {7, 8},
            {8, 9},
            {8, 19},
            {9, 10},
            {9, 15},
            {10, 11},
            {11, 12},
            {12, 13},
            {13, 14},
            {13, 16},
            {15, 11},
            {15, 13},
            {17, 1},
            {17, 2},
            {17, 18},
            {18, 3},
            {19, 10},
            {19, 15}
        };

        std::vector<std::vector<int>> graph2_paths{
            {0, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14},
            {17, 1, 2, 3, 4, 5, 19, 15, 13, 16},
            {17, 1, 18, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 16}
        };

        for (auto e : graph2_edges) {
            graph2.add_edge(e.first, e.second);
        }

        for (size_t i = 0; i < graph2_paths.size(); ++i) {
            auto p = graph2.add_path(std::to_string(i));
            for (auto n : graph2_paths[i]) {
                graph2.extend_path(p, n);
            }
        }

        auto anchors = generate_anchor_set(graph1, graph2, 2);
        test_sparse_dynamic_programming(graph1, graph2, anchors, 15, 6, 4, 13, false, true, false);
    }

    {
        BaseGraph graph1;
        for (auto c : std::string("CAACCCAATCCAACCCAACCCCACCAACAG")) {
            graph1.add_node(c);
        }

        std::vector<std::pair<int, int>> graph1_edges{
            {0, 1},
            {1, 2},
            {2, 3},
            {2, 4},
            {2, 8},
            {3, 4},
            {4, 5},
            {5, 6},
            {6, 7},
            {6, 25},
            {7, 8},
            {8, 9},
            {8, 10},
            {9, 10},
            {10, 11},
            {11, 12},
            {12, 13},
            {12, 18},
            {13, 14},
            {14, 15},
            {14, 26},
            {15, 16},
            {16, 17},
            {16, 28},
            {16, 19},
            {17, 18},
            {18, 19},
            {19, 20},
            {19, 29},
            {20, 21},
            {21, 22},
            {22, 23},
            {22, 27},
            {23, 24},
            {25, 8},
            {26, 16},
            {27, 24},
            {28, 18},
            {29, 21}
        };

        std::vector<std::vector<int>> graph1_paths{
            {0, 1, 2, 4, 5, 6, 25, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 29, 21, 22, 23, 24},
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 26, 16, 28, 18, 19, 20, 21, 22, 27, 24}
        };

        for (auto e : graph1_edges) {
            graph1.add_edge(e.first, e.second);
        }

        for (size_t i = 0; i < graph1_paths.size(); ++i) {
            auto p = graph1.add_path(std::to_string(i));
            for (auto n : graph1_paths[i]) {
                graph1.extend_path(p, n);
            }
        }

        BaseGraph graph2;
        for (auto c : std::string("TTAAAAAAAAAAAAAAAAAAAAAAAATAAA")) {
            graph2.add_node(c);
        }

        std::vector<std::pair<int, int>> graph2_edges{
            {0, 1},
            {1, 2},
            {1, 26},
            {1, 3},
            {2, 3},
            {3, 4},
            {3, 27},
            {4, 5},
            {4, 25},
            {5, 6},
            {5, 12},
            {6, 7},
            {7, 8},
            {8, 9},
            {9, 10},
            {10, 11},
            {11, 12},
            {12, 13},
            {13, 14},
            {14, 15},
            {15, 16},
            {16, 17},
            {16, 29},
            {16, 18},
            {16, 22},
            {17, 18},
            {18, 19},
            {19, 20},
            {20, 21},
            {20, 22},
            {20, 23},
            {21, 22},
            {21, 23},
            {22, 23},
            {23, 24},
            {23, 28},
            {25, 6},
            {26, 3},
            {27, 5},
            {27, 25},
            {27, 8},
            {29, 18}
        };

        std::vector<std::vector<int>> graph2_paths{
            {0, 1, 26, 3, 4, 25, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24},
            {0, 1, 2, 3, 27, 5, 12, 13, 14, 15, 16, 29, 18, 19, 20, 21, 23, 28}
        };

        for (auto e : graph2_edges) {
            graph2.add_edge(e.first, e.second);
        }

        for (size_t i = 0; i < graph2_paths.size(); ++i) {
            auto p = graph2.add_path(std::to_string(i));
            for (auto n : graph2_paths[i]) {
                graph2.extend_path(p, n);
            }
        }

        auto anchors = generate_anchor_set(graph1, graph2, 2);
        test_sparse_dynamic_programming(graph1, graph2, anchors, 29, 23, 17, 19, true, false, false);
    }

    {
        BaseGraph graph1;
        for (auto c : std::string("CCCCCCCCCCTCCCCACCCCCCCGCAATAT")) {
            graph1.add_node(c);
        }

        std::vector<std::pair<int, int>> graph1_edges{
            {0, 1},
            {1, 2},
            {1, 3},
            {2, 3},
            {3, 4},
            {4, 5},
            {5, 6},
            {6, 7},
            {7, 8},
            {8, 9},
            {9, 10},
            {9, 26},
            {9, 16},
            {10, 11},
            {11, 12},
            {12, 13},
            {12, 29},
            {13, 14},
            {14, 15},
            {15, 16},
            {15, 19},
            {16, 17},
            {16, 23},
            {17, 18},
            {18, 19},
            {19, 20},
            {19, 28},
            {20, 21},
            {20, 25},
            {21, 22},
            {22, 23},
            {22, 24},
            {23, 24},
            {25, 22},
            {26, 11},
            {27, 1},
            {28, 21},
            {28, 25},
            {29, 14}
        };

        std::vector<std::vector<int>> graph1_paths{
            {27, 1, 3, 4, 5, 6, 7, 8, 9, 26, 11, 12, 29, 14, 15, 19, 28, 25, 22, 23, 24},
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24}
        };

        for (auto e : graph1_edges) {
            graph1.add_edge(e.first, e.second);
        }

        for (size_t i = 0; i < graph1_paths.size(); ++i) {
            auto p = graph1.add_path(std::to_string(i));
            for (auto n : graph1_paths[i]) {
                graph1.extend_path(p, n);
            }
        }

        BaseGraph graph2;
        for (auto c : std::string("CCCACCCAACCCCACCCGACCCCACGACGC")) {
            graph2.add_node(c);
        }

        std::vector<std::pair<int, int>> graph2_edges{
            {0, 1},
            {1, 2},
            {2, 3},
            {3, 4},
            {4, 5},
            {5, 6},
            {5, 25},
            {6, 7},
            {7, 8},
            {7, 10},
            {8, 9},
            {8, 25},
            {8, 10},
            {9, 10},
            {10, 11},
            {10, 28},
            {11, 12},
            {12, 13},
            {12, 27},
            {13, 14},
            {13, 26},
            {14, 15},
            {14, 28},
            {15, 16},
            {16, 17},
            {16, 28},
            {17, 18},
            {18, 19},
            {19, 20},
            {20, 21},
            {21, 22},
            {22, 23},
            {22, 29},
            {23, 24},
            {25, 10},
            {26, 15},
            {27, 14},
            {27, 26},
            {28, 18},
            {29, 24}
        };

        std::vector<std::vector<int>> graph2_paths{
            {0, 1, 2, 3, 4, 5, 6, 7, 10, 11, 12, 27, 14, 28, 18, 19, 20, 21, 22, 23, 24},
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 26, 15, 16, 17, 18, 19, 20, 21, 22, 29, 24},
            {0, 1, 2, 3, 4, 5, 25, 10, 28, 18, 19, 20, 21, 22, 23, 24}
        };

        for (auto e : graph2_edges) {
            graph2.add_edge(e.first, e.second);
        }

        for (size_t i = 0; i < graph2_paths.size(); ++i) {
            auto p = graph2.add_path(std::to_string(i));
            for (auto n : graph2_paths[i]) {
                graph2.extend_path(p, n);
            }
        }

        auto anchors = generate_anchor_set(graph1, graph2, 2);
        test_sparse_dynamic_programming(graph1, graph2, anchors, 2, 3, 6, 10, true, true, false);
    }

    {
        BaseGraph graph1;
        for (auto c : std::string("TACCTCTGTTCACTACGTAC")) {
            graph1.add_node(c);
        }

        std::vector<std::pair<int, int>> graph1_edges{
            {0, 9},
            {0, 13},
            {0, 3},
            {0, 2},
            {0, 17},
            {0, 14},
            {0, 10},
            {1, 9},
            {1, 18},
            {1, 7},
            {1, 11},
            {1, 5},
            {1, 6},
            {1, 12},
            {1, 16},
            {2, 5},
            {3, 8},
            {3, 7},
            {3, 6},
            {3, 17},
            {3, 4},
            {3, 18},
            {3, 13},
            {4, 19},
            {4, 5},
            {4, 8},
            {4, 18},
            {4, 6},
            {5, 9},
            {5, 6},
            {5, 14},
            {5, 7},
            {5, 18},
            {5, 8},
            {5, 12},
            {5, 10},
            {5, 19},
            {6, 8},
            {6, 11},
            {6, 16},
            {6, 15},
            {6, 7},
            {7, 16},
            {7, 14},
            {7, 18},
            {7, 11},
            {7, 12},
            {7, 9},
            {8, 14},
            {8, 9},
            {8, 10},
            {9, 11},
            {9, 18},
            {9, 10},
            {9, 16},
            {9, 13},
            {9, 17},
            {9, 15},
            {10, 15},
            {10, 16},
            {10, 12},
            {10, 11},
            {11, 18},
            {11, 14},
            {11, 13},
            {11, 17},
            {11, 12},
            {12, 19},
            {13, 14},
            {13, 15},
            {13, 16},
            {13, 17},
            {13, 19},
            {13, 18},
            {14, 16},
            {15, 18},
            {16, 18},
            {16, 17},
            {17, 19},
            {18, 19}
        };

        std::vector<std::vector<int>> graph1_paths{
            {0, 3, 7, 11, 18, 19},
            {1, 5, 14, 16, 17, 19},
            {0, 3, 4, 6, 8, 9, 15, 18, 19},
            {0, 3, 4, 8, 9, 10, 12, 19},
            {0, 2, 5, 12, 19},
            {0, 3, 13, 19}
        };

        for (auto e : graph1_edges) {
            graph1.add_edge(e.first, e.second);
        }

        for (size_t i = 0; i < graph1_paths.size(); ++i) {
            auto p = graph1.add_path(std::to_string(i));
            for (auto n : graph1_paths[i]) {
                graph1.extend_path(p, n);
            }
        }

        BaseGraph graph2;
        for (auto c : std::string("TATTTCCACTTAAACACCTA")) {
            graph2.add_node(c);
        }

        std::vector<std::pair<int, int>> graph2_edges{
            {0, 3},
            {0, 19},
            {0, 10},
            {0, 16},
            {0, 17},
            {1, 13},
            {1, 15},
            {1, 17},
            {1, 4},
            {1, 18},
            {1, 12},
            {1, 16},
            {2, 3},
            {2, 14},
            {2, 17},
            {2, 9},
            {2, 8},
            {2, 4},
            {2, 18},
            {3, 12},
            {3, 17},
            {3, 7},
            {3, 15},
            {3, 8},
            {3, 9},
            {3, 5},
            {3, 13},
            {4, 8},
            {4, 5},
            {4, 10},
            {4, 6},
            {4, 19},
            {4, 18},
            {4, 12},
            {4, 16},
            {5, 13},
            {5, 10},
            {5, 18},
            {5, 9},
            {6, 8},
            {6, 7},
            {6, 18},
            {6, 13},
            {6, 17},
            {6, 9},
            {6, 19},
            {7, 11},
            {7, 19},
            {7, 18},
            {7, 17},
            {7, 16},
            {8, 13},
            {8, 10},
            {8, 12},
            {8, 15},
            {8, 11},
            {9, 12},
            {9, 10},
            {9, 19},
            {10, 13},
            {10, 18},
            {10, 14},
            {11, 15},
            {11, 16},
            {11, 19},
            {11, 12},
            {12, 13},
            {12, 15},
            {13, 14},
            {13, 19},
            {13, 18},
            {13, 17},
            {14, 16},
            {14, 15},
            {14, 17},
            {14, 18},
            {15, 16},
            {15, 19},
            {16, 18},
            {16, 19}
        };

        std::vector<std::vector<int>> graph2_paths{
            {1, 4, 16, 18},
            {2, 3, 8, 12, 13, 19},
            {0, 17},
            {1, 4, 6, 7, 11, 15, 19},
            {2, 3, 5, 9, 10, 14, 18}
        };

        for (auto e : graph2_edges) {
            graph2.add_edge(e.first, e.second);
        }

        for (size_t i = 0; i < graph2_paths.size(); ++i) {
            auto p = graph2.add_path(std::to_string(i));
            for (auto n : graph2_paths[i]) {
                graph2.extend_path(p, n);
            }
        }

        auto anchors = generate_anchor_set(graph1, graph2, 2);
        test_sparse_dynamic_programming(graph1, graph2, anchors, 3, 18, 2, 19, true, true, false);
    }

    {
        BaseGraph graph1;
        for (auto c : std::string("ACACTACGATATGAATCAGG")) {
            graph1.add_node(c);
        }

        std::vector<std::pair<int, int>> graph1_edges{
            {0, 3},
            {0, 17},
            {0, 6},
            {0, 12},
            {0, 9},
            {0, 7},
            {0, 2},
            {0, 13},
            {1, 14},
            {1, 12},
            {1, 15},
            {1, 9},
            {1, 5},
            {1, 19},
            {2, 11},
            {2, 18},
            {2, 7},
            {2, 17},
            {2, 15},
            {2, 19},
            {2, 8},
            {3, 8},
            {3, 9},
            {3, 14},
            {3, 5},
            {3, 19},
            {3, 7},
            {4, 13},
            {4, 6},
            {4, 19},
            {4, 18},
            {5, 14},
            {5, 13},
            {5, 9},
            {5, 6},
            {5, 8},
            {5, 7},
            {5, 17},
            {6, 13},
            {6, 14},
            {6, 12},
            {6, 9},
            {6, 10},
            {6, 18},
            {7, 15},
            {7, 17},
            {8, 12},
            {8, 16},
            {8, 18},
            {8, 10},
            {8, 11},
            {8, 14},
            {9, 13},
            {9, 14},
            {9, 11},
            {9, 19},
            {9, 18},
            {9, 10},
            {10, 16},
            {10, 18},
            {10, 13},
            {10, 19},
            {11, 12},
            {11, 18},
            {11, 14},
            {11, 17},
            {11, 19},
            {11, 15},
            {12, 13},
            {12, 15},
            {12, 17},
            {13, 18},
            {13, 19},
            {14, 19},
            {15, 19},
            {15, 18},
            {16, 19},
            {16, 17},
            {17, 19},
            {18, 19}
        };

        std::vector<std::vector<int>> graph1_paths{
            {0, 3, 5, 6, 10, 13, 19},
            {0, 2, 8, 12, 17, 19},
            {1, 5, 7, 15, 18, 19},
            {0, 3, 8, 16, 17, 19},
            {4, 13, 18, 19},
            {1, 9, 14, 19},
            {0, 3, 8, 11, 12, 15, 19}
        };

        for (auto e : graph1_edges) {
            graph1.add_edge(e.first, e.second);
        }

        for (size_t i = 0; i < graph1_paths.size(); ++i) {
            auto p = graph1.add_path(std::to_string(i));
            for (auto n : graph1_paths[i]) {
                graph1.extend_path(p, n);
            }
        }

        BaseGraph graph2;
        for (auto c : std::string("CTGATGGCTTAGAGACAAGC")) {
            graph2.add_node(c);
        }

        std::vector<std::pair<int, int>> graph2_edges{
            {0, 11},
            {0, 17},
            {0, 3},
            {0, 14},
            {0, 13},
            {0, 19},
            {0, 2},
            {0, 7},
            {1, 7},
            {1, 10},
            {1, 8},
            {1, 18},
            {1, 5},
            {1, 12},
            {1, 3},
            {1, 6},
            {1, 19},
            {2, 9},
            {2, 17},
            {2, 10},
            {2, 4},
            {2, 18},
            {2, 5},
            {2, 6},
            {2, 13},
            {2, 3},
            {2, 12},
            {2, 7},
            {3, 15},
            {3, 19},
            {3, 18},
            {3, 10},
            {3, 11},
            {3, 16},
            {4, 5},
            {4, 9},
            {4, 8},
            {5, 13},
            {5, 7},
            {5, 8},
            {5, 14},
            {5, 12},
            {6, 16},
            {6, 15},
            {6, 18},
            {6, 10},
            {6, 7},
            {7, 14},
            {7, 19},
            {8, 12},
            {8, 14},
            {8, 11},
            {8, 13},
            {8, 10},
            {8, 19},
            {9, 12},
            {9, 18},
            {9, 14},
            {9, 16},
            {10, 11},
            {10, 15},
            {10, 13},
            {10, 14},
            {11, 14},
            {11, 16},
            {11, 12},
            {11, 17},
            {12, 18},
            {12, 16},
            {12, 19},
            {12, 14},
            {13, 14},
            {13, 17},
            {13, 15},
            {14, 19},
            {14, 16},
            {14, 17},
            {15, 19},
            {17, 18},
            {18, 19}
        };

        std::vector<std::vector<int>> graph2_paths{
            {0, 2, 6, 15, 19},
            {1, 12, 14, 17, 18, 19},
            {0, 2, 4, 5, 7, 14, 16},
            {0, 2, 4, 8, 10, 13, 14, 19},
            {1, 3, 11, 16},
            {0, 2, 9, 18, 19}
        };

        for (auto e : graph2_edges) {
            graph2.add_edge(e.first, e.second);
        }

        for (size_t i = 0; i < graph2_paths.size(); ++i) {
            auto p = graph2.add_path(std::to_string(i));
            for (auto n : graph2_paths[i]) {
                graph2.extend_path(p, n);
            }
        }

        auto anchors = generate_anchor_set(graph1, graph2, 3);
        test_sparse_dynamic_programming(graph1, graph2, anchors, 3, 16, 1, 17, true, true, false);
    }

    {
        BaseGraph graph1;
        for (auto c : std::string("GCATCTACTG")) {
            graph1.add_node(c);
        }

        std::vector<std::pair<int, int>> graph1_edges{
            {0, 9},
            {0, 1},
            {0, 5},
            {1, 5},
            {1, 3},
            {1, 6},
            {2, 8},
            {3, 8},
            {3, 7},
            {4, 5},
            {4, 7},
            {5, 6},
            {5, 8},
            {5, 9},
            {5, 7},
            {6, 8},
            {6, 9},
            {7, 8}
        };

        std::vector<std::vector<int>> graph1_paths{
            {2, 8},
            {0, 1, 3, 7, 8},
            {4, 5, 6, 9}
        };

        for (auto e : graph1_edges) {
            graph1.add_edge(e.first, e.second);
        }

        for (size_t i = 0; i < graph1_paths.size(); ++i) {
            auto p = graph1.add_path(std::to_string(i));
            for (auto n : graph1_paths[i]) {
                graph1.extend_path(p, n);
            }
        }

        BaseGraph graph2;
        for (auto c : std::string("ACCTATGACC")) {
            graph2.add_node(c);
        }

        std::vector<std::pair<int, int>> graph2_edges{
            {0, 6},
            {0, 2},
            {0, 4},
            {1, 8},
            {1, 9},
            {1, 7},
            {1, 6},
            {2, 5},
            {2, 6},
            {2, 3},
            {3, 7},
            {3, 9},
            {3, 5},
            {3, 8},
            {4, 6},
            {4, 8},
            {5, 6},
            {5, 8}
        };

        std::vector<std::vector<int>> graph2_paths{
            {0, 4, 8},
            {0, 2, 3, 7},
            {0, 2, 3, 5, 6},
            {1, 9}
        };

        for (auto e : graph2_edges) {
            graph2.add_edge(e.first, e.second);
        }

        for (size_t i = 0; i < graph2_paths.size(); ++i) {
            auto p = graph2.add_path(std::to_string(i));
            for (auto n : graph2_paths[i]) {
                graph2.extend_path(p, n);
            }
        }

        auto anchors = generate_anchor_set(graph1, graph2, 2);
        test_sparse_dynamic_programming(graph1, graph2, anchors, 5, 7, 3, 7, true, true, false);
    }

    {
        BaseGraph graph1;
        for (auto c : std::string("TCGTCGA")) {
            graph1.add_node(c);
        }

        std::vector<std::pair<int, int>> graph1_edges{
            {0, 5},
            {0, 1},
            {0, 2},
            {0, 3},
            {1, 4},
            {2, 3},
            {2, 4},
            {3, 6},
            {4, 5},
            {4, 6}
        };

        std::vector<std::vector<int>> graph1_paths{
            {0, 3, 6},
            {0, 2, 4, 5},
            {0, 1, 4, 5}
        };

        for (auto e : graph1_edges) {
            graph1.add_edge(e.first, e.second);
        }

        for (size_t i = 0; i < graph1_paths.size(); ++i) {
            auto p = graph1.add_path(std::to_string(i));
            for (auto n : graph1_paths[i]) {
                graph1.extend_path(p, n);
            }
        }

        BaseGraph graph2;
        for (auto c : std::string("ACGAAGC")) {
            graph2.add_node(c);
        }

        std::vector<std::pair<int, int>> graph2_edges{
            {0, 6},
            {0, 3},
            {0, 4},
            {0, 1},
            {1, 5},
            {1, 6},
            {1, 2},
            {2, 5},
            {4, 5},
            {5, 6}
        };

        std::vector<std::vector<int>> graph2_paths{
            {0, 1, 2, 5, 6},
            {0, 3},
            {0, 4, 5, 6}
        };

        for (auto e : graph2_edges) {
            graph2.add_edge(e.first, e.second);
        }

        for (size_t i = 0; i < graph2_paths.size(); ++i) {
            auto p = graph2.add_path(std::to_string(i));
            for (auto n : graph2_paths[i]) {
                graph2.extend_path(p, n);
            }
        }

        auto anchors = generate_anchor_set(graph1, graph2, 2);
        test_sparse_dynamic_programming(graph1, graph2, anchors, 4, 6, 0, 6, false, true, false);
    }

    {
        BaseGraph graph1;
        for (auto c : std::string("GACTTTT")) {
            graph1.add_node(c);
        }

        std::vector<std::pair<int, int>> graph1_edges{
            {0, 3},
            {1, 5},
            {1, 3},
            {1, 6},
            {1, 2},
            {2, 5},
            {2, 3},
            {3, 4},
            {4, 6},
            {5, 6}
        };

        std::vector<std::vector<int>> graph1_paths{
            {1, 3, 4, 6},
            {0, 3, 4, 6},
            {1, 2, 5, 6}
        };

        for (auto e : graph1_edges) {
            graph1.add_edge(e.first, e.second);
        }

        for (size_t i = 0; i < graph1_paths.size(); ++i) {
            auto p = graph1.add_path(std::to_string(i));
            for (auto n : graph1_paths[i]) {
                graph1.extend_path(p, n);
            }
        }

        BaseGraph graph2;
        for (auto c : std::string("GGATCTA")) {
            graph2.add_node(c);
        }

        std::vector<std::pair<int, int>> graph2_edges{
            {0, 3},
            {1, 5},
            {2, 3},
            {2, 6},
            {2, 4},
            {3, 5},
            {3, 6},
            {3, 4},
            {4, 5},
            {4, 6}
        };

        std::vector<std::vector<int>> graph2_paths{
            {2, 6},
            {0, 3, 4, 5},
            {1, 5}
        };

        for (auto e : graph2_edges) {
            graph2.add_edge(e.first, e.second);
        }

        for (size_t i = 0; i < graph2_paths.size(); ++i) {
            auto p = graph2.add_path(std::to_string(i));
            for (auto n : graph2_paths[i]) {
                graph2.extend_path(p, n);
            }
        }

        auto anchors = generate_anchor_set(graph1, graph2, 2);
        test_sparse_dynamic_programming(graph1, graph2, anchors, 3, 6, 0, 6, false, true, false);
    }

    {
        BaseGraph graph1;
        for (auto c : std::string("CGAACGT")) {
            graph1.add_node(c);
        }

        std::vector<std::pair<int, int>> graph1_edges{
            {0, 4},
            {0, 3},
            {1, 2},
            {1, 5},
            {2, 5},
            {2, 3},
            {3, 4},
            {3, 5},
            {4, 5},
            {4, 6}
        };

        std::vector<std::vector<int>> graph1_paths{
            {1, 2, 3, 5},
            {0, 4, 6}
        };

        for (auto e : graph1_edges) {
            graph1.add_edge(e.first, e.second);
        }

        for (size_t i = 0; i < graph1_paths.size(); ++i) {
            auto p = graph1.add_path(std::to_string(i));
            for (auto n : graph1_paths[i]) {
                graph1.extend_path(p, n);
            }
        }

        BaseGraph graph2;
        for (auto c : std::string("TGTGAGC")) {
            graph2.add_node(c);
        }

        std::vector<std::pair<int, int>> graph2_edges{
            {0, 4},
            {1, 4},
            {1, 2},
            {2, 4},
            {2, 3},
            {2, 5},
            {3, 4},
            {3, 6},
            {3, 5},
            {4, 5}
        };

        std::vector<std::vector<int>> graph2_paths{
            {1, 2, 3, 6},
            {0, 4, 5}
        };

        for (auto e : graph2_edges) {
            graph2.add_edge(e.first, e.second);
        }

        for (size_t i = 0; i < graph2_paths.size(); ++i) {
            auto p = graph2.add_path(std::to_string(i));
            for (auto n : graph2_paths[i]) {
                graph2.extend_path(p, n);
            }
        }

        auto anchors = generate_anchor_set(graph1, graph2, 2);
        test_sparse_dynamic_programming(graph1, graph2, anchors, 2, 6, 1, 4, false, true, false);
    }

    {
        BaseGraph graph1;
        for (auto c : std::string("TCCTTCC")) {
            graph1.add_node(c);
        }

        std::vector<std::pair<int, int>> graph1_edges{
            {0, 4},
            {0, 2},
            {0, 3},
            {1, 3},
            {2, 3},
            {2, 4},
            {3, 5},
            {3, 4},
            {3, 6},
            {4, 5}
        };

        std::vector<std::vector<int>> graph1_paths{
            {1, 3, 6},
            {0, 2, 4, 5}
        };

        for (auto e : graph1_edges) {
            graph1.add_edge(e.first, e.second);
        }

        for (size_t i = 0; i < graph1_paths.size(); ++i) {
            auto p = graph1.add_path(std::to_string(i));
            for (auto n : graph1_paths[i]) {
                graph1.extend_path(p, n);
            }
        }

        BaseGraph graph2;
        for (auto c : std::string("AGGTTCG")) {
            graph2.add_node(c);
        }

        std::vector<std::pair<int, int>> graph2_edges{
            {0, 3},
            {0, 1},
            {0, 6},
            {0, 2},
            {0, 5},
            {1, 3},
            {1, 4},
            {3, 5},
            {3, 6},
            {4, 6}
        };

        std::vector<std::vector<int>> graph2_paths{
            {0, 5},
            {0, 1, 4, 6},
            {0, 3, 5},
            {0, 2}
        };

        for (auto e : graph2_edges) {
            graph2.add_edge(e.first, e.second);
        }

        for (size_t i = 0; i < graph2_paths.size(); ++i) {
            auto p = graph2.add_path(std::to_string(i));
            for (auto n : graph2_paths[i]) {
                graph2.extend_path(p, n);
            }
        }

        auto anchors = generate_anchor_set(graph1, graph2, 2);
        test_sparse_dynamic_programming(graph1, graph2, anchors, 0, 3, 4, 6, false, true, false);
    }

    {
        BaseGraph graph1;
        for (auto c : std::string("CCAAACCAACCGAACCCAACCCAACCCACT")) {
            graph1.add_node(c);
        }

        std::vector<std::pair<int, int>> graph1_edges{
            {0, 1},
            {0, 27},
            {1, 2},
            {2, 3},
            {3, 4},
            {4, 5},
            {5, 6},
            {6, 7},
            {7, 8},
            {7, 11},
            {8, 9},
            {8, 25},
            {8, 28},
            {8, 14},
            {9, 10},
            {9, 11},
            {9, 26},
            {10, 11},
            {10, 26},
            {11, 12},
            {12, 13},
            {13, 14},
            {14, 15},
            {14, 19},
            {15, 16},
            {15, 17},
            {16, 17},
            {17, 18},
            {18, 19},
            {19, 20},
            {20, 21},
            {20, 29},
            {21, 22},
            {22, 23},
            {22, 24},
            {23, 24},
            {25, 10},
            {25, 11},
            {25, 26},
            {26, 12},
            {27, 2},
            {27, 4},
            {28, 10},
            {28, 11},
            {28, 26},
            {29, 22}
        };

        std::vector<std::vector<int>> graph1_paths{
            {0, 27, 2, 3, 4, 5, 6, 7, 8, 9, 10, 26, 12, 13, 14, 15, 17, 18, 19, 20, 29, 22, 23, 24},
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 28, 11, 12, 13, 14, 19, 20, 21, 22, 23, 24},
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 25, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 24}
        };

        for (auto e : graph1_edges) {
            graph1.add_edge(e.first, e.second);
        }

        for (size_t i = 0; i < graph1_paths.size(); ++i) {
            auto p = graph1.add_path(std::to_string(i));
            for (auto n : graph1_paths[i]) {
                graph1.extend_path(p, n);
            }
        }

        BaseGraph graph2;
        for (auto c : std::string("AAACAAAACAAAACAAAGCAAAACAATACC")) {
            graph2.add_node(c);
        }

        std::vector<std::pair<int, int>> graph2_edges{
            {0, 1},
            {0, 28},
            {1, 2},
            {2, 3},
            {3, 4},
            {4, 5},
            {5, 6},
            {6, 7},
            {7, 8},
            {7, 29},
            {8, 9},
            {9, 10},
            {10, 11},
            {11, 12},
            {11, 15},
            {12, 13},
            {13, 14},
            {14, 15},
            {14, 25},
            {15, 16},
            {16, 17},
            {16, 18},
            {16, 26},
            {17, 18},
            {17, 26},
            {18, 19},
            {18, 27},
            {19, 20},
            {20, 21},
            {21, 22},
            {22, 23},
            {23, 24},
            {25, 16},
            {25, 21},
            {26, 19},
            {26, 27},
            {27, 20},
            {28, 2},
            {29, 9}
        };

        std::vector<std::vector<int>> graph2_paths{
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 15, 16, 17, 26, 27, 20, 21, 22, 23, 24},
            {0, 28, 2, 3, 4, 5, 6, 7, 29, 9, 10, 11, 12, 13, 14, 25, 16, 18, 19, 20, 21, 22, 23, 24}
        };

        for (auto e : graph2_edges) {
            graph2.add_edge(e.first, e.second);
        }

        for (size_t i = 0; i < graph2_paths.size(); ++i) {
            auto p = graph2.add_path(std::to_string(i));
            for (auto n : graph2_paths[i]) {
                graph2.extend_path(p, n);
            }
        }

        auto anchors = generate_anchor_set(graph1, graph2, 2);
        test_sparse_dynamic_programming(graph1, graph2, anchors, -1, -1, -1, -1, true, false, false);
    }

    {
        BaseGraph graph1;
        for (auto c : std::string("CATTGTCCTAAGGAAT")) {
            graph1.add_node(c);
        }

        std::vector<std::pair<int, int>> graph1_edges{
            {0, 15},
            {0, 6},
            {0, 14},
            {0, 11},
            {1, 7},
            {1, 10},
            {1, 6},
            {1, 8},
            {1, 9},
            {2, 3},
            {2, 14},
            {2, 15},
            {2, 11},
            {3, 13},
            {3, 9},
            {4, 11},
            {4, 6},
            {5, 11},
            {5, 6},
            {5, 15},
            {5, 10},
            {6, 15},
            {6, 13},
            {7, 15},
            {7, 11},
            {8, 10},
            {8, 12},
            {9, 15},
            {9, 12},
            {9, 11},
            {10, 13},
            {10, 14},
            {11, 13},
            {13, 14},
            {14, 15}
        };

        std::vector<std::vector<int>> graph1_paths{
            {0, 11, 13, 14, 15},
            {1, 7, 11, 13, 14, 15},
            {2, 3, 9, 12},
            {1, 8, 10, 14, 15},
            {4, 6, 15},
            {5, 10, 14, 15}
        };

        for (auto e : graph1_edges) {
            graph1.add_edge(e.first, e.second);
        }

        for (size_t i = 0; i < graph1_paths.size(); ++i) {
            auto p = graph1.add_path(std::to_string(i));
            for (auto n : graph1_paths[i]) {
                graph1.extend_path(p, n);
            }
        }

        BaseGraph graph2;
        for (auto c : std::string("TGGCTGTGGCACGGAG")) {
            graph2.add_node(c);
        }

        std::vector<std::pair<int, int>> graph2_edges{
            {0, 1},
            {1, 11},
            {1, 13},
            {1, 2},
            {2, 14},
            {2, 10},
            {2, 13},
            {2, 5},
            {3, 12},
            {4, 6},
            {4, 12},
            {5, 15},
            {5, 7},
            {5, 8},
            {5, 10},
            {6, 9},
            {6, 8},
            {6, 15},
            {6, 12},
            {7, 9},
            {7, 12},
            {8, 15},
            {8, 9},
            {8, 14},
            {9, 10},
            {9, 13},
            {9, 12},
            {9, 11},
            {10, 14},
            {10, 11},
            {11, 14},
            {11, 13},
            {11, 15},
            {12, 13},
            {13, 15}
        };

        std::vector<std::vector<int>> graph2_paths{
            {0, 1, 11, 15},
            {4, 12, 13, 15},
            {4, 6, 9, 10, 14},
            {3, 12, 13, 15},
            {0, 1, 2, 5, 8, 14},
            {0, 1, 2, 5, 7, 9, 10, 11, 14}
        };

        for (auto e : graph2_edges) {
            graph2.add_edge(e.first, e.second);
        }

        for (size_t i = 0; i < graph2_paths.size(); ++i) {
            auto p = graph2.add_path(std::to_string(i));
            for (auto n : graph2_paths[i]) {
                graph2.extend_path(p, n);
            }
        }

        auto anchors = generate_anchor_set(graph1, graph2, 2);
        test_sparse_dynamic_programming(graph1, graph2, anchors, -1, -1, -1, -1, true, false, false);
    }

    {
        BaseGraph graph1;
        for (auto c : std::string("AAGAAGCACCTCAACT")) {
            graph1.add_node(c);
        }

        std::vector<std::pair<int, int>> graph1_edges{
            {0, 6},
            {0, 11},
            {0, 8},
            {0, 3},
            {1, 2},
            {1, 5},
            {1, 9},
            {1, 13},
            {1, 11},
            {2, 4},
            {2, 11},
            {2, 8},
            {3, 6},
            {3, 8},
            {3, 5},
            {4, 13},
            {4, 9},
            {5, 8},
            {5, 15},
            {5, 10},
            {5, 6},
            {6, 10},
            {6, 15},
            {6, 14},
            {6, 13},
            {6, 11},
            {7, 15},
            {7, 13},
            {8, 10},
            {8, 12},
            {9, 15},
            {10, 12},
            {11, 14},
            {11, 13},
            {13, 14}
        };

        std::vector<std::vector<int>> graph1_paths{
            {1, 5, 10, 12},
            {1, 2, 4, 9, 15},
            {0, 3, 8, 12},
            {0, 6, 11, 13, 14},
            {7, 13, 14}
        };

        for (auto e : graph1_edges) {
            graph1.add_edge(e.first, e.second);
        }

        for (size_t i = 0; i < graph1_paths.size(); ++i) {
            auto p = graph1.add_path(std::to_string(i));
            for (auto n : graph1_paths[i]) {
                graph1.extend_path(p, n);
            }
        }

        BaseGraph graph2;
        for (auto c : std::string("TAGGAGCGCTCGTGGA")) {
            graph2.add_node(c);
        }

        std::vector<std::pair<int, int>> graph2_edges{
            {0, 4},
            {0, 10},
            {0, 14},
            {0, 7},
            {0, 12},
            {0, 15},
            {1, 5},
            {1, 4},
            {1, 10},
            {1, 14},
            {1, 13},
            {1, 7},
            {3, 7},
            {3, 6},
            {4, 9},
            {4, 14},
            {4, 11},
            {4, 15},
            {4, 6},
            {5, 6},
            {5, 10},
            {5, 14},
            {5, 8},
            {5, 12},
            {6, 10},
            {6, 14},
            {7, 11},
            {7, 15},
            {9, 13},
            {9, 15},
            {9, 12},
            {10, 13},
            {10, 14},
            {11, 12},
            {12, 15}
        };

        std::vector<std::vector<int>> graph2_paths{
            {1, 5, 14},
            {0, 4, 6, 10, 13},
            {3, 7, 15},
            {2},
            {1, 4, 11, 12, 15},
            {0, 4, 9, 12, 15},
            {1, 5, 8}
        };

        for (auto e : graph2_edges) {
            graph2.add_edge(e.first, e.second);
        }

        for (size_t i = 0; i < graph2_paths.size(); ++i) {
            auto p = graph2.add_path(std::to_string(i));
            for (auto n : graph2_paths[i]) {
                graph2.extend_path(p, n);
            }
        }

        auto anchors = generate_anchor_set(graph1, graph2, 2);
        test_sparse_dynamic_programming(graph1, graph2, anchors, -1, -1, -1, -1, true, false, false);
    }

    {
        BaseGraph graph1;
        for (auto c : std::string("ATTCGGGCAC")) {
            graph1.add_node(c);
        }

        std::vector<std::pair<int, int>> graph1_edges{
            {0, 6},
            {0, 1},
            {0, 3},
            {0, 8},
            {0, 4},
            {1, 5},
            {1, 2},
            {1, 7},
            {1, 8},
            {1, 4},
            {2, 9},
            {3, 6},
            {3, 9},
            {4, 9},
            {5, 6},
            {6, 7},
            {6, 8},
            {8, 9}
        };

        std::vector<std::vector<int>> graph1_paths{
            {0, 6, 7},
            {0, 1, 8, 9},
            {0, 3, 9},
            {0, 1, 5, 6, 8, 9},
            {0, 1, 2, 9},
            {0, 1, 4, 9}
        };

        for (auto e : graph1_edges) {
            graph1.add_edge(e.first, e.second);
        }

        for (size_t i = 0; i < graph1_paths.size(); ++i) {
            auto p = graph1.add_path(std::to_string(i));
            for (auto n : graph1_paths[i]) {
                graph1.extend_path(p, n);
            }
        }

        BaseGraph graph2;
        for (auto c : std::string("AGCTCAAACA")) {
            graph2.add_node(c);
        }

        std::vector<std::pair<int, int>> graph2_edges{
            {0, 1},
            {0, 2},
            {1, 4},
            {1, 2},
            {1, 5},
            {2, 3},
            {2, 5},
            {3, 7},
            {3, 9},
            {4, 5},
            {4, 8},
            {5, 7},
            {5, 8},
            {5, 9},
            {6, 8},
            {6, 9},
            {7, 9},
            {8, 9}
        };

        std::vector<std::vector<int>> graph2_paths{
            {0, 2, 3, 7, 9},
            {0, 1, 5, 8, 9},
            {6, 9},
            {0, 1, 4, 5, 8, 9}
        };

        for (auto e : graph2_edges) {
            graph2.add_edge(e.first, e.second);
        }

        for (size_t i = 0; i < graph2_paths.size(); ++i) {
            auto p = graph2.add_path(std::to_string(i));
            for (auto n : graph2_paths[i]) {
                graph2.extend_path(p, n);
            }
        }

        auto anchors = generate_anchor_set(graph1, graph2, 2);
        test_sparse_dynamic_programming(graph1, graph2, anchors, -1, -1, -1, -1, true, false, false);
    }

    {
        BaseGraph graph1;
        for (auto c : std::string("GGCAGTTAGCCATCGC")) {
            graph1.add_node(c);
        }

        std::vector<std::pair<int, int>> graph1_edges{
            {0, 2},
            {0, 8},
            {0, 6},
            {0, 4},
            {1, 10},
            {1, 9},
            {1, 13},
            {1, 7},
            {1, 4},
            {2, 8},
            {2, 7},
            {2, 6},
            {2, 13},
            {3, 8},
            {3, 4},
            {3, 11},
            {3, 9},
            {4, 6},
            {4, 9},
            {4, 11},
            {5, 9},
            {5, 7},
            {6, 14},
            {6, 10},
            {7, 8},
            {7, 10},
            {8, 15},
            {8, 13},
            {8, 14},
            {9, 10},
            {9, 11},
            {10, 11},
            {12, 15},
            {12, 13},
            {13, 15}
        };

        std::vector<std::vector<int>> graph1_paths{
            {1, 7, 10, 11},
            {0, 4, 9, 11},
            {0, 2, 13, 15},
            {3, 8, 14},
            {1, 4, 6, 10, 11},
            {12, 13, 15},
            {5, 9, 10, 11}
        };

        for (auto e : graph1_edges) {
            graph1.add_edge(e.first, e.second);
        }

        for (size_t i = 0; i < graph1_paths.size(); ++i) {
            auto p = graph1.add_path(std::to_string(i));
            for (auto n : graph1_paths[i]) {
                graph1.extend_path(p, n);
            }
        }
        BaseGraph graph2;
        for (auto c : std::string("ACTACATCGGCCGCCG")) {
            graph2.add_node(c);
        }

        std::vector<std::pair<int, int>> graph2_edges{
            {0, 13},
            {0, 15},
            {0, 11},
            {0, 3},
            {1, 8},
            {1, 12},
            {1, 3},
            {1, 7},
            {1, 4},
            {1, 13},
            {2, 3},
            {2, 6},
            {2, 10},
            {3, 15},
            {3, 8},
            {3, 12},
            {4, 13},
            {4, 14},
            {4, 8},
            {5, 6},
            {5, 12},
            {5, 15},
            {5, 11},
            {6, 12},
            {6, 8},
            {6, 7},
            {6, 14},
            {7, 15},
            {8, 14},
            {9, 10},
            {10, 12},
            {11, 14},
            {11, 12},
            {12, 13},
            {14, 15}
        };

        std::vector<std::vector<int>> graph2_paths{
            {5, 12, 13},
            {2, 6, 7, 15},
            {0, 3, 8, 14, 15},
            {9, 10, 12, 13},
            {1, 4, 14, 15},
            {0, 11, 14, 15}
        };

        for (auto e : graph2_edges) {
            graph2.add_edge(e.first, e.second);
        }

        for (size_t i = 0; i < graph2_paths.size(); ++i) {
            auto p = graph2.add_path(std::to_string(i));
            for (auto n : graph2_paths[i]) {
                graph2.extend_path(p, n);
            }
        }

        auto anchors = generate_anchor_set(graph1, graph2, 2);
        test_sparse_dynamic_programming(graph1, graph2, anchors, -1, -1, -1, -1, true, false, false);
    }

    // sparse chaining tests
    {
        // baby 2-bubble graphs

        BaseGraph graph1;
        uint64_t n10 = graph1.add_node('A');
        uint64_t n11 = graph1.add_node('A');
        uint64_t n12 = graph1.add_node('G');
        uint64_t n13 = graph1.add_node('T');
        uint64_t n14 = graph1.add_node('G');
        uint64_t n15 = graph1.add_node('C');
        uint64_t n16 = graph1.add_node('C');

        graph1.add_edge(n10, n11);
        graph1.add_edge(n10, n12);
        graph1.add_edge(n11, n13);
        graph1.add_edge(n12, n13);
        graph1.add_edge(n13, n14);
        graph1.add_edge(n13, n15);
        graph1.add_edge(n14, n16);
        graph1.add_edge(n15, n16);

        uint64_t p10 = graph1.add_path("1");
        graph1.extend_path(p10, n10);
        graph1.extend_path(p10, n11);
        graph1.extend_path(p10, n13);
        graph1.extend_path(p10, n14);
        graph1.extend_path(p10, n16);
        uint64_t p11 = graph1.add_path("2");
        graph1.extend_path(p11, n10);
        graph1.extend_path(p11, n12);
        graph1.extend_path(p11, n13);
        graph1.extend_path(p11, n15);
        graph1.extend_path(p11, n16);

        BaseGraph graph2;
        uint64_t n20 = graph2.add_node('G');
        uint64_t n21 = graph2.add_node('A');
        uint64_t n22 = graph2.add_node('G');
        uint64_t n23 = graph2.add_node('A');
        uint64_t n24 = graph2.add_node('T');
        uint64_t n25 = graph2.add_node('C');
        uint64_t n26 = graph2.add_node('C');

        graph2.add_edge(n20, n21);
        graph2.add_edge(n20, n22);
        graph2.add_edge(n21, n23);
        graph2.add_edge(n22, n23);
        graph2.add_edge(n23, n24);
        graph2.add_edge(n23, n25);
        graph2.add_edge(n24, n26);
        graph2.add_edge(n25, n26);

        uint64_t p20 = graph2.add_path("3");
        graph2.extend_path(p20, n20);
        graph2.extend_path(p20, n21);
        graph2.extend_path(p20, n23);
        graph2.extend_path(p20, n25);
        graph2.extend_path(p20, n26);
        uint64_t p21 = graph2.add_path("4");
        graph2.extend_path(p21, n20);
        graph2.extend_path(p21, n22);
        graph2.extend_path(p21, n23);
        graph2.extend_path(p21, n24);
        graph2.extend_path(p21, n26);

        {
            vector<match_set_t> anchors(2);
            anchors[0].walks1.emplace_back(vector<uint64_t>{n10, n11});
            anchors[0].walks2.emplace_back(vector<uint64_t>{n21, n23});
            anchors[1].walks1.emplace_back(vector<uint64_t>{n15, n16});
            anchors[1].walks2.emplace_back(vector<uint64_t>{n25, n26});


            anchors[0].count1 = anchors[0].walks1.size();
            anchors[0].count2 = anchors[0].walks2.size();
            anchors[0].full_length = anchors[0].walks1.front().size();
            anchors[1].count1 = anchors[1].walks1.size();
            anchors[1].count2 = anchors[1].walks2.size();
            anchors[1].full_length = anchors[1].walks1.front().size();

            for (auto affine : {true, false}) {
                for (auto packed : {false, true}) {
                    test_sparse_dynamic_programming(graph1, graph2, anchors, -1, -1, -1, -1, affine, false, packed);
                }
            }
        }

        // these aren't real matches anymore, but the algorithm doesn't pay
        // attention to the sequence, so whatever
        {
            // put them out of order on graph2
            vector<match_set_t> anchors(2);
            anchors[0].walks1.emplace_back(vector<uint64_t>{n10, n11});
            anchors[0].walks2.emplace_back(vector<uint64_t>{n25, n26});
            anchors[1].walks1.emplace_back(vector<uint64_t>{n15, n16});
            anchors[1].walks2.emplace_back(vector<uint64_t>{n21, n23});

            anchors[0].count1 = anchors[0].walks1.size();
            anchors[0].count2 = anchors[0].walks2.size();
            anchors[0].full_length = anchors[0].walks1.front().size();
            anchors[1].count1 = anchors[1].walks1.size();
            anchors[1].count2 = anchors[1].walks2.size();
            anchors[1].full_length = anchors[1].walks1.front().size();

            for (auto affine : {true, false}) {
                for (auto packed : {true, false}) {
                    test_sparse_dynamic_programming(graph1, graph2, anchors, -1, -1, -1, -1, affine, false, packed);
                }
            }
        }
        {
            // make it have to choose a better score
            vector<match_set_t> anchors(3);
            anchors[0].walks1.emplace_back(vector<uint64_t>{n10, n11});
            anchors[0].walks2.emplace_back(vector<uint64_t>{n20, n22});
            anchors[1].walks1.emplace_back(vector<uint64_t>{n14});
            anchors[1].walks2.emplace_back(vector<uint64_t>{n24});
            anchors[2].walks1.emplace_back(vector<uint64_t>{n15});
            anchors[2].walks2.emplace_back(vector<uint64_t>{n25});
            anchors[2].walks2.emplace_back(vector<uint64_t>{n26});

            anchors[0].count1 = anchors[0].walks1.size();
            anchors[0].count2 = anchors[0].walks2.size();
            anchors[0].full_length = anchors[0].walks1.front().size();
            anchors[1].count1 = anchors[1].walks1.size();
            anchors[1].count2 = anchors[1].walks2.size();
            anchors[1].full_length = anchors[1].walks1.front().size();
            anchors[2].count1 = anchors[2].walks1.size();
            anchors[2].count2 = anchors[2].walks2.size();
            anchors[2].full_length = anchors[2].walks1.front().size();

            for (auto affine : {true, false}) {
                for (auto packed : {true, false}) {
                    test_sparse_dynamic_programming(graph1, graph2, anchors, -1, -1, -1, -1, affine, false, packed);
                }
            }
        }
    }

    // sparse affine dynamic programming on simple linear graphs
    {
        BaseGraph graph1 = make_base_graph("1", "AACCGTT");
        BaseGraph graph2 = make_base_graph("2", "AAGCCTT");

        vector<match_set_t> anchors(3);
        anchors[0].walks1.emplace_back(vector<uint64_t>{0, 1});
        anchors[0].walks2.emplace_back(vector<uint64_t>{0, 1});
        anchors[0].count1 = 1;
        anchors[0].count2 = 1;
        anchors[0].full_length = 2;
        anchors[1].walks1.emplace_back(vector<uint64_t>{2, 3});
        anchors[1].walks2.emplace_back(vector<uint64_t>{3, 4});
        anchors[1].count1 = 1;
        anchors[1].count2 = 1;
        anchors[1].full_length = 2;
        anchors[2].walks1.emplace_back(vector<uint64_t>{5, 6});
        anchors[2].walks2.emplace_back(vector<uint64_t>{5, 6});
        anchors[2].count1 = 1;
        anchors[2].count2 = 1;
        anchors[2].full_length = 2;

        ScoreFunction score_function;
        score_function.anchor_score_function = ScoreFunction::InverseCount;
        TestAnchorer anchorer(score_function);
        anchorer.gap_open[0] = 1.0;
        anchorer.gap_extend[0] = 1.0;
        anchorer.gap_open[1] = 3.0;
        anchorer.gap_extend[1] = 0.5;
        anchorer.gap_open[2] = 5.0;
        anchorer.gap_extend[2] = 5.0;

        PathMerge<> chain_merge1(graph1);
        PathMerge<> chain_merge2(graph2);

        auto chain = anchorer.sparse_affine_chain_dp<size_t, size_t, size_t, int64_t, size_t, float, std::vector<std::pair<int64_t, MatchBank<size_t, size_t>::match_id_t>>, std::vector<std::pair<size_t, MatchBank<size_t, size_t>::match_id_t>>, std::vector<size_t>, std::vector<size_t>, MatchBank<size_t, size_t>, ForwardEdges<uint64_t, uint64_t>>
                                                    (anchors, graph1, graph2, chain_merge1, chain_merge2,anchorer.gap_open,
                                                     anchorer.gap_extend, 1.0,anchors.size(), true);
        
        bool correct = (chain.size() == 2);
        correct &= (chain[0].walk1 == vector<uint64_t>{0, 1});
        correct &= (chain[0].walk2 == vector<uint64_t>{0, 1});
        correct &= (chain[1].walk1 == vector<uint64_t>{5, 6});
        correct &= (chain[1].walk2 == vector<uint64_t>{5, 6});
        assert(correct);
    }
    
    // test the fill in anchoring algorithm
    {
        BaseGraph graph1;
        for (auto c : std::string("GATCGAT")) {
            graph1.add_node(c);
        }
        
        std::vector<std::pair<int, int>> graph1_edges{
            {0, 1},
            {1, 2},
            {2, 3},
            {3, 4},
            {4, 5},
            {5, 6}
        };
        
        std::vector<std::vector<int>> graph1_paths{
            {0, 1, 2, 3, 4, 5, 6}
        };
        
        for (auto e : graph1_edges) {
            graph1.add_edge(e.first, e.second);
        }
        
        for (size_t i = 0; i < graph1_paths.size(); ++i) {
            auto p = graph1.add_path(std::to_string(i));
            for (auto n : graph1_paths[i]) {
                graph1.extend_path(p, n);
            }
        }
        
        BaseGraph graph2;
        for (auto c : std::string("TAGCTAG")) {
            graph2.add_node(c);
        }
        
        std::vector<std::pair<int, int>> graph2_edges{
            {0, 1},
            {1, 2},
            {2, 3},
            {3, 4},
            {4, 5},
            {5, 6}
        };
        
        std::vector<std::vector<int>> graph2_paths{
            {0, 1, 2, 3, 4, 5, 6}
        };
        
        for (auto e : graph2_edges) {
            graph2.add_edge(e.first, e.second);
        }
        
        for (size_t i = 0; i < graph2_paths.size(); ++i) {
            auto p = graph2.add_path(std::to_string(i));
            for (auto n : graph2_paths[i]) {
                graph2.extend_path(p, n);
            }
        }
        
        BaseGraph sub1;
        for (auto c : std::string("GAT")) {
            sub1.add_node(c);
        }
        
        std::vector<std::pair<int, int>> sub1_edges{
            {0, 1},
            {1, 2}
        };
        
        std::vector<std::vector<int>> sub1_paths{
            {0, 1, 2}
        };
        
        for (auto e : sub1_edges) {
            sub1.add_edge(e.first, e.second);
        }
        
        for (size_t i = 0; i < sub1_paths.size(); ++i) {
            auto p = sub1.add_path(std::to_string(i));
            for (auto n : sub1_paths[i]) {
                sub1.extend_path(p, n);
            }
        }
        
        BaseGraph sub2;
        for (auto c : std::string("TAG")) {
            sub2.add_node(c);
        }
        
        std::vector<std::pair<int, int>> sub2_edges{
            {0, 1},
            {1, 2}
        };
        
        std::vector<std::vector<int>> sub2_paths{
            {0, 1, 2}
        };
        
        for (auto e : sub2_edges) {
            sub2.add_edge(e.first, e.second);
        }
        
        for (size_t i = 0; i < sub2_paths.size(); ++i) {
            auto p = sub2.add_path(std::to_string(i));
            for (auto n : sub2_paths[i]) {
                sub2.extend_path(p, n);
            }
        }
        
        auto tableau1 = add_sentinels(graph1, '^', '$');
        auto tableau2 = add_sentinels(graph2, '^', '$');
        
        std::vector<match_set_t> match_sets(4);
        // A
        match_sets[0].walks1.push_back({1});
        match_sets[0].walks1.push_back({5});
        match_sets[0].walks2.push_back({1});
        match_sets[0].walks2.push_back({5});
        // C
        match_sets[1].walks1.push_back({3});
        match_sets[1].walks2.push_back({3});
        // G
        match_sets[2].walks1.push_back({0});
        match_sets[2].walks1.push_back({4});
        match_sets[2].walks2.push_back({2});
        match_sets[2].walks2.push_back({6});
        // T
        match_sets[3].walks1.push_back({2});
        match_sets[3].walks1.push_back({6});
        match_sets[3].walks2.push_back({0});
        match_sets[3].walks2.push_back({4});
        
        for (auto& m : match_sets) {
            m.count1 = m.walks1.size();
            m.count2 = m.walks2.size();
            m.full_length = 1;
        }
        
        std::vector<anchor_t> chain(1);
        chain.front().walk1 = match_sets[1].walks1.front();
        chain.front().walk2 = match_sets[1].walks2.front();
        chain.front().count1 = 1;
        chain.front().count2 = 1;
        
        ChainMerge merge1(graph1, tableau1);
        ChainMerge merge2(graph2, tableau2);
        
        ScoreFunction score_function;
        score_function.anchor_score_function = ScoreFunction::InverseCount;
        TestAnchorer anchorer(score_function);
        
        auto stitch_graphs = anchorer.extract_graphs_between(chain, graph1, graph2,
                                                             tableau1, tableau2,
                                                             merge1, merge2);
        
        anchorer.project_paths(graph1, graph2, stitch_graphs);
        
        std::vector<std::vector<std::pair<size_t, std::pair<std::vector<size_t>, std::vector<size_t>>>>> origins;
        auto divvied = anchorer.divvy_matches(match_sets, graph1, graph2, stitch_graphs, origins);
        
        assert(stitch_graphs.size() == 2);
        assert(divvied.size() == 2);
        
        for (size_t i = 0; i < stitch_graphs.size(); ++i) {
            
            auto& stitch_graph = stitch_graphs[i];
            
            assert(stitch_graph.first.back_translation.size() == stitch_graph.first.subgraph.node_size());
            assert(stitch_graph.second.back_translation.size() == stitch_graph.second.subgraph.node_size());
            
            if (!possibly_isomorphic(sub1, stitch_graph.first.subgraph)) {
                cerr << "isomorphism 1 failure\n";
                exit(1);
            }
            if (!possibly_isomorphic(sub2, stitch_graph.second.subgraph)) {
                cerr << "isomorphism 2 failure\n";
                exit(1);
            }
            
            assert(divvied[i].size() == 3);
            for (auto& set : divvied[i]) {
                assert(set.walks1.size() == 1);
                assert(set.walks2.size() == 1);
            }
        }
    }
    

    vector<pair<int, int>> graph_sizes{{7, 10}, {10, 18}, {16, 25}, {16, 35}, {20, 80}};
    for (auto size : graph_sizes) {
        for (int rep = 0; rep < 5; ++rep) {
            BaseGraph graph1 = random_graph(size.first, size.second, true, gen);
            BaseGraph graph2 = random_graph(size.first, size.second, true, gen);
            add_random_path_cover(graph1, gen);
            add_random_path_cover(graph2, gen);
//            cerr << "graph1:\n";
//            cerr << cpp_representation(graph1, "graph1");
//            cerr << "graph2:\n";
//            cerr << cpp_representation(graph2, "graph2");
            for (int k : {2, 3}) {
                auto anchors = generate_anchor_set(graph1, graph2, k);
                for (auto affine : {true, false}) {
                    for (auto global : {true, false}) {
                        for (auto packed : {true, false}) {
                            auto src_snk1 = generate_source_sink(graph1, gen);
                            auto src_snk2 = generate_source_sink(graph2, gen);
                            test_sparse_dynamic_programming(graph1, graph2, anchors,
                                                            src_snk1.first, src_snk1.second,
                                                            src_snk2.first, src_snk2.second,
                                                            affine, global, packed);
                        }
                    }
                }
            }
        }
    }

    vector<int> challenge_graph_sizes{10, 20, 30};
    for (auto size : challenge_graph_sizes) {
        for (int rep = 0; rep < 5; ++rep) {
            BaseGraph graph1 = random_challenge_graph(size, gen);
            BaseGraph graph2 = random_challenge_graph(size, gen);
            add_random_path_cover(graph1, gen);
            add_random_path_cover(graph2, gen);
            for (int k : {2, 3, 4}) {
                auto anchors = generate_anchor_set(graph1, graph2, k);
                for (auto affine : {true, false}) {
                    for (auto global : {true, false}) {
                        for (auto packed : {true, false}) {
                            auto src_snk1 = generate_source_sink(graph1, gen);
                            auto src_snk2 = generate_source_sink(graph2, gen);
                            test_sparse_dynamic_programming(graph1, graph2, anchors,
                                                            src_snk1.first, src_snk1.second,
                                                            src_snk2.first, src_snk2.second,
                                                            affine, global, packed);
                        }
                    }
                }
            }
        }
    }
    
    cerr << "passed all tests!" << endl;
}
