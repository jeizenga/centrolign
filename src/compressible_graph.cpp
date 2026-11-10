#include "centrolign/compressible_graph.hpp"

#include <stdexcept>
#include <chrono>
#include <algorithm>

#include "centrolign/logging.hpp"
#include "centrolign/utility.hpp"

namespace centrolign {

using namespace std;

PackedSequenceGraph::PackedSequenceGraph(const BaseGraph& graph, const SentinelTableau& tableau,
                                         std::vector<uint64_t>* translation) {

    has_sentinels = (tableau.src_id < graph.node_size() && tableau.snk_id < graph.node_size() && tableau.src_id != tableau.snk_id);

    auto is_sentinel = [&](uint64_t node_id) {
        return node_id == tableau.src_id || node_id == tableau.snk_id;
    };

    // identify the nodes that cannot be merged into their predecessor
    vector<bool> is_head(graph.node_size(), true);
    for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
        if (graph.previous_size(node_id) == 1 && !is_sentinel(node_id)) {
            uint64_t prev_id = graph.previous(node_id).front();
            if (prev_id != node_id && graph.next_size(prev_id) == 1 && !is_sentinel(prev_id)) {
                is_head[node_id] = false;
            }
        }
    }
    // paths must start at the beginning of a sequence node and end at the end of one
    for (uint64_t path_id = 0; path_id < graph.path_size(); ++path_id) {
        const auto& path = graph.path(path_id);
        if (!path.empty()) {
            is_head[path.front()] = true;
            if (graph.next_size(path.back()) == 1) {
                is_head[graph.next(path.back()).front()] = true;
            }
        }
    }

    // a node that is not the end of a sequence node has exactly one successor, which follows it
    auto continues = [&](uint64_t node_id) {
        return graph.next_size(node_id) == 1 && !is_head[graph.next(node_id).front()];
    };

    // assign decompressed node IDs so that each sequence node is a consecutive range
    vector<uint64_t> new_ids(graph.node_size(), -1);
    // the index of the sequence node that contains each original node
    vector<uint64_t> sequence_node(graph.node_size());
    vector<uint64_t> tails;
    uint64_t next_new_id = 0;
    auto add_sequence_node = [&](uint64_t node_id) {
        new_ids[node_id] = next_new_id++;
        sequence_node[node_id] = tails.size();
        while (continues(node_id)) {
            node_id = graph.next(node_id).front();
            new_ids[node_id] = next_new_id++;
            sequence_node[node_id] = tails.size();
        }
        tails.push_back(node_id);
    };

    for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
        if (is_head[node_id] && !is_sentinel(node_id)) {
            add_sequence_node(node_id);
        }
    }
    for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
        if (new_ids[node_id] == -1 && !is_sentinel(node_id)) {
            // this is part of a non-branching cycle, which has no natural head
            is_head[node_id] = true;
            add_sequence_node(node_id);
        }
    }
    if (has_sentinels) {
        add_sequence_node(tableau.src_id);
        add_sequence_node(tableau.snk_id);
        src_label = graph.label(tableau.src_id);
        snk_label = graph.label(tableau.snk_id);
    }

    // encode the labels
    uint64_t max_label = 0;
    for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
        if (!is_sentinel(node_id)) {
            max_label = max<uint64_t>(max_label, uint8_t(graph.label(node_id)));
        }
    }
    labels = SizedVectorFactory<PackedVector>::make(graph.node_size(), max_label);
    for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
        if (!is_sentinel(node_id)) {
            labels[new_ids[node_id]] = uint8_t(graph.label(node_id));
        }
    }

    // encode the sequence node boundaries
    node_starts = SizedVectorFactory<PackedVector>::make(tails.size() + 1, graph.node_size());
    for (size_t i = 0; i + 1 < tails.size(); ++i) {
        node_starts[i + 1] = new_ids[tails[i]] + 1;
    }
    node_starts[tails.size()] = graph.node_size();

    // encode the edges between sequence nodes
    size_t num_edges = 0;
    size_t max_degree = 0;
    for (auto tail_id : tails) {
        num_edges += graph.next_size(tail_id);
        max_degree = max(max_degree, graph.next_size(tail_id));
    }
    out_degrees = SizedVectorFactory<PackedVector>::make(tails.size(), max_degree);
    edge_targets = SizedVectorFactory<PackedVector>::make(num_edges, tails.empty() ? 0 : tails.size() - 1);
    for (size_t i = 0, e = 0; i < tails.size(); ++i) {
        const auto& next_ids = graph.next(tails[i]);
        out_degrees[i] = next_ids.size();
        for (auto next_id : next_ids) {
            edge_targets[e++] = sequence_node[next_id];
        }
    }

    // encode the paths
    path_names.reserve(graph.path_size());
    paths.resize(graph.path_size());
    for (uint64_t path_id = 0; path_id < graph.path_size(); ++path_id) {

        path_names.emplace_back(graph.path_name(path_id));

        const auto& path = graph.path(path_id);
        auto& packed_path = paths[path_id];
        if (path.empty()) {
            continue;
        }
        packed_path.start = sequence_node[path.front()];

        // record which edge is taken whenever the path leaves a sequence node
        vector<uint64_t> choices;
        uint64_t max_choice = 0;
        for (size_t i = 1; i < path.size(); ++i) {
            const auto& next_ids = graph.next(path[i - 1]);
            auto it = find(next_ids.begin(), next_ids.end(), path[i]);
            if (it == next_ids.end()) {
                throw runtime_error("Path " + graph.path_name(path_id) + " does not follow an edge from node " + to_string(path[i - 1]) + " to node " + to_string(path[i]) + ", so it cannot be compressed");
            }
            if (!continues(path[i - 1])) {
                choices.push_back(it - next_ids.begin());
                max_choice = max<uint64_t>(max_choice, choices.back());
            }
        }

        packed_path.steps = SizedVectorFactory<PackedVector>::make(choices.size(), max_choice);
        for (size_t i = 0; i < choices.size(); ++i) {
            packed_path.steps[i] = choices[i];
        }
    }

    if (translation) {
        *translation = move(new_ids);
    }
}

BaseGraph PackedSequenceGraph::decompress(SentinelTableau& tableau) const {

    BaseGraph graph;

    for (uint64_t node_id = 0; node_id < labels.size(); ++node_id) {
        graph.add_node(char(labels.at(node_id)));
    }
    if (has_sentinels) {
        tableau.src_id = labels.size() - 2;
        tableau.snk_id = labels.size() - 1;
        graph.relabel(tableau.src_id, src_label);
        graph.relabel(tableau.snk_id, snk_label);
    }

    // the offset of each sequence node's edges in the edge targets
    vector<uint64_t> edge_offsets(out_degrees.size() + 1, 0);
    for (size_t i = 0, e = 0; i < out_degrees.size(); ++i) {
        uint64_t begin = node_starts.at(i), end = node_starts.at(i + 1);
        for (uint64_t node_id = begin + 1; node_id < end; ++node_id) {
            graph.add_edge(node_id - 1, node_id);
        }
        for (size_t j = 0, n = out_degrees.at(i); j < n; ++j) {
            graph.add_edge(end - 1, node_starts.at(edge_targets.at(e++)));
        }
        edge_offsets[i + 1] = e;
    }

    for (size_t i = 0; i < paths.size(); ++i) {
        uint64_t path_id = graph.add_path(path_names[i]);
        const auto& packed_path = paths[i];
        if (packed_path.start != -1) {
            // the path covers each sequence node it visits entirely
            uint64_t seq_node = packed_path.start;
            for (size_t c = 0; ; ++c) {
                for (uint64_t node_id = node_starts.at(seq_node), end = node_starts.at(seq_node + 1); node_id < end; ++node_id) {
                    graph.extend_path(path_id, node_id);
                }
                if (c == packed_path.steps.size()) {
                    break;
                }
                seq_node = edge_targets.at(edge_offsets[seq_node] + packed_path.steps.at(c));
            }
        }
    }

    return graph;
}

size_t PackedSequenceGraph::node_size() const {
    return labels.size();
}

size_t PackedSequenceGraph::memory_size() const {
    size_t size = sizeof(*this);
    size += labels.memory_size() + node_starts.memory_size() + out_degrees.memory_size() + edge_targets.memory_size();
    size += path_names.capacity() * sizeof(decltype(path_names)::value_type);
    for (const auto& name : path_names) {
        size += name.capacity();
    }
    size += paths.capacity() * sizeof(decltype(paths)::value_type);
    for (const auto& path : paths) {
        size += path.steps.memory_size();
    }
    return size;
}

CompressibleGraph::CompressibleGraph(BaseGraph&& graph, const SentinelTableau& tableau) noexcept
    : graph(move(graph)), tableau(tableau) {
    cached_memory_size.store(this->graph.memory_size());
}

CompressibleGraph::CompressibleGraph(const CompressibleGraph& other) noexcept {
    *this = other;
}

CompressibleGraph::CompressibleGraph(CompressibleGraph&& other) noexcept {
    *this = move(other);
}

CompressibleGraph& CompressibleGraph::operator=(const CompressibleGraph& other) noexcept {
    if (this != &other) {
        graph = other.graph;
        tableau = other.tableau;
        packed_graph = other.packed_graph;
        compressed = other.compressed;
        cached_memory_size.store(other.cached_memory_size.load());
    }
    return *this;
}

CompressibleGraph& CompressibleGraph::operator=(CompressibleGraph&& other) noexcept {
    if (this != &other) {
        graph = move(other.graph);
        tableau = other.tableau;
        packed_graph = move(other.packed_graph);
        compressed = other.compressed;
        cached_memory_size.store(other.cached_memory_size.load());

        other.graph = BaseGraph();
        other.packed_graph = PackedSequenceGraph();
        other.compressed = false;
        other.cached_memory_size.store(0);
    }
    return *this;
}

void CompressibleGraph::compress() {
    if (compressed) {
        return;
    }
    auto start = chrono::steady_clock::now();

    packed_graph = PackedSequenceGraph(graph, tableau);
    // release the memory
    graph = BaseGraph();
    compressed = true;

    size_t size_before = cached_memory_size.load();
    cached_memory_size.store(packed_graph.memory_size());

    if (logging::level >= logging::Debug) {
        double secs = chrono::duration<double>(chrono::steady_clock::now() - start).count();
        logging::log(logging::Debug, "Compressed graph with " + to_string(packed_graph.node_size()) + " nodes from " + format_memory_usage(size_before) + " to " + format_memory_usage(cached_memory_size.load()) + " in " + to_string(secs) + " s.");
    }
}

void CompressibleGraph::decompress() {
    if (!compressed) {
        return;
    }
    auto start = chrono::steady_clock::now();

    graph = packed_graph.decompress(tableau);
    // release the memory
    packed_graph = PackedSequenceGraph();
    compressed = false;

    cached_memory_size.store(graph.memory_size());

    if (logging::level >= logging::Debug) {
        double secs = chrono::duration<double>(chrono::steady_clock::now() - start).count();
        logging::log(logging::Debug, "Decompressed graph with " + to_string(graph.node_size()) + " nodes to " + format_memory_usage(cached_memory_size.load()) + " in " + to_string(secs) + " s.");
    }
}

BaseGraph& CompressibleGraph::get_base_graph() {
    decompress();
    return graph;
}

SentinelTableau& CompressibleGraph::get_tableau() {
    decompress();
    return tableau;
}

const BaseGraph& CompressibleGraph::get_base_graph() const {
    if (compressed) {
        throw runtime_error("Attempted to access a compressed graph without decompressing it");
    }
    return graph;
}

const SentinelTableau& CompressibleGraph::get_tableau() const {
    if (compressed) {
        throw runtime_error("Attempted to access a compressed graph's tableau without decompressing it");
    }
    return tableau;
}

bool CompressibleGraph::is_compressed() const {
    return compressed;
}

size_t CompressibleGraph::node_size() const {
    return compressed ? packed_graph.node_size() : graph.node_size();
}

size_t CompressibleGraph::memory_size() const {
    return cached_memory_size.load();
}

}
