#include "centrolign/modify_graph.hpp"

#include <unordered_set>
#include <vector>
#include <cassert>
#include <map>

#include "centrolign/utility.hpp"
#include "centrolign/logging.hpp"
#include "centrolign/step_index.hpp"


namespace centrolign {

using namespace std;

SequenceGraph make_sequence_graph(const std::string& name,
                                  const std::string& sequence) {
    assert(!sequence.empty());
    assert(!name.empty());
    SequenceGraph graph;
    uint64_t node_id = graph.add_node(encode_seq(sequence));
    uint64_t path_id = graph.add_path(name);
    graph.extend_path(path_id, node_id);
    return graph;
}

BaseGraph make_base_graph(const std::string& name,
                          const std::string& sequence) {
    assert(!sequence.empty());
    assert(!name.empty());
    BaseGraph graph;
    uint64_t prev_id = graph.add_node(encode_base(sequence[0]));
    uint64_t path_id = graph.add_path(name);
    graph.extend_path(path_id, prev_id);
    for (size_t i = 1; i < sequence.size(); ++i) {
        uint64_t node_id = graph.add_node(encode_base(sequence[i]));
        graph.add_edge(prev_id, node_id);
        graph.extend_path(path_id, node_id);
        prev_id = node_id;
    }
    return graph;
}

SentinelTableau add_sentinels(BaseGraph& graph,
                              char src_sentinel, char snk_sentinel) {
    
    std::unordered_set<uint64_t> path_endings, path_beginnings;
    for (uint64_t path_id = 0; path_id < graph.path_size(); ++path_id) {
        path_beginnings.insert(graph.path(path_id).front());
        path_endings.insert(graph.path(path_id).back());
    }
    
    SentinelTableau tableau;
    tableau.src_id = graph.add_node(src_sentinel);
    tableau.snk_id = graph.add_node(snk_sentinel);
    tableau.src_sentinel = src_sentinel;
    tableau.snk_sentinel = snk_sentinel;
    
    if (graph.node_size() == 2) {
        graph.add_edge(tableau.src_id, tableau.snk_id);
    }
    else {
        for (uint64_t node_id = 0; node_id < tableau.src_id; ++node_id) {
            if (graph.next_size(node_id) == 0 || path_endings.count(node_id)) {
                graph.add_edge(node_id, tableau.snk_id);
            }
            if (graph.previous_size(node_id) == 0 || path_beginnings.count(node_id)) {
                graph.add_edge(tableau.src_id, node_id);
            }
        }
    }
    
    return tableau;
}

void reassign_sentinels(BaseGraph& graph, SentinelTableau& tableau,
                        char src_sentinel, char snk_sentinel) {
    
    tableau.src_sentinel = src_sentinel;
    tableau.snk_sentinel = snk_sentinel;
    graph.relabel(tableau.src_id, src_sentinel);
    graph.relabel(tableau.snk_id, snk_sentinel);
}


void purge_uncovered_nodes(BaseGraph& graph, SentinelTableau& tableau) {
    
    static const bool debug = false;
    
    vector<bool> covered(graph.node_size(), false);
    covered[tableau.src_id] = true;
    covered[tableau.snk_id] = true;
    
    for (uint64_t path_id = 0; path_id < graph.path_size(); ++path_id) {
        for (auto node_id : graph.path(path_id)) {
            covered[node_id] = true;
        }
    }
    
    bool all_covered = true;
    for (size_t i = 0; i < covered.size() && all_covered; ++i) {
        all_covered = covered[i];
    }
    
    if (debug) {
        cerr << "covered vector:\n";
        for (size_t i = 0; i < covered.size(); ++i) {
            cerr << '\t' << i << ": " << covered[i] << '\n';
        }
    }
    
    if (!all_covered) {
        
        logging::log(logging::Debug, "Not all nodes are covered by paths, purging");
        
        BaseGraph purged;
        vector<uint64_t> removed_before(covered.size() + 1, 0);
        for (size_t i = 0; i < covered.size(); ++i) {
            if (covered[i]) {
                purged.add_node(graph.label(i));
                removed_before[i + 1] = removed_before[i];
            }
            else {
                removed_before[i + 1] = removed_before[i] + 1;
            }
        }
        
        if (debug) {
            cerr << "removed_before vector:\n";
            for (size_t i = 0; i < removed_before.size(); ++i) {
                cerr << '\t' << i << ": " << removed_before[i] << '\n';
            }
        }
        
        for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
            if (covered[node_id]) {
                for (auto next_id : graph.next(node_id)) {
                    if (covered[next_id]) {
                        purged.add_edge(node_id - removed_before[node_id],
                                        next_id - removed_before[next_id]);
                    }
                }
            }
        }
        
        for (uint64_t path_id = 0; path_id < graph.path_size(); ++path_id) {
            uint64_t new_path_id = purged.add_path(graph.path_name(path_id));
            for (auto node_id : graph.path(path_id)) {
                purged.extend_path(new_path_id, node_id - removed_before[node_id]);
            }
        }
        
        tableau.src_id -= removed_before[tableau.src_id];
        tableau.snk_id -= removed_before[tableau.snk_id];
        
        logging::log(logging::Debug, "Removed " + to_string(graph.node_size() - purged.node_size()) + " uncovered nodes");
        
        graph = move(purged);
    }
}

void simplify_bubbles(BaseGraph& graph, SentinelTableau& tableau) {
    
    static const bool debug = false;
    
    StepIndex step_index(graph);
    
    bool did_simplify = false;
    
    for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
        
        if (graph.next_size(node_id) <= 1) {
            // not a bubble source
            continue;
        }
        
        if (debug) {
            std::cerr << "checking for bubble at source node " << node_id << '\n';
        }
        
        // rigid test for a simple, possibly multialleleic single base bubble
        bool is_bubble = true;
        uint64_t sink_id = -1;
        for (auto next_id : graph.next(node_id)) {
            if (next_id == sink_id) {
                // deletion allele
                is_bubble = false;
                if (debug) {
                    std::cerr << "fails deletion test on node " << next_id << '\n';
                }
                break;
            }
            // all alleles should be nonbranching
            if (graph.next_size(next_id) != 1) {
                is_bubble = false;
                if (debug) {
                    std::cerr << "fails nonbranching test on test to node " << next_id << '\n';
                }
                break;
            }
            // all alleles should hit the same sink
            if (sink_id == -1) {
                sink_id = graph.next(next_id).front();
                if (debug) {
                    std::cerr << "sink candidate is " << sink_id << '\n';
                }
            }
            else if (graph.next(next_id).front() != sink_id) {
                is_bubble = false;
                if (debug) {
                    std::cerr << "fails single sink test on test to node " << next_id << '\n';
                }
                break;
            }
        }
        if (debug) {
            if (is_bubble && sink_id == -1) {
                std::cerr << "failed to find sink\n";
            }
        }
        // did we find a sink?
        is_bubble = is_bubble && sink_id != -1;
        if (is_bubble) {
            // does the sink have any other edges incoming?
            is_bubble = graph.previous_size(sink_id) == graph.next_size(node_id);
            if (debug) {
                if (!is_bubble) {
                    std::cerr << "failed equal degree test\n";
                }
            }
        }
        
        if (is_bubble) {
            if (debug) {
                std::cerr << "verified bubble starting at " << node_id << '\n';
            }
            // aggregate alleles with the same label
            // note: ordered map for machine independence
            std::map<char, std::vector<uint64_t>> label_sets;
            for (auto next_id : graph.next(node_id)) {
                label_sets[graph.label(next_id)].push_back(next_id);
            }
            
            for (const auto& label_set : label_sets) {
                if (label_set.second.size() > 1) {
                    // there are multiple alleles with this label, reassign the paths to one node
                    if (debug) {
                        std::cerr << "label set for " << label_set.first << " can be collapsed" << '\n';
                    }
                    for (size_t j = 1; j < label_set.second.size(); ++j) {
                        for (const auto& step : step_index.path_steps(label_set.second[j])) {
                            graph.reassign_subpath(step.first, step.second, {label_set.second.front()});
                        }
                    }
                    did_simplify = true;
                }
            }
        }
    }
    
    if (did_simplify) {
        // we removed paths from some nodes
        purge_uncovered_nodes(graph, tableau);
    }
}

}
