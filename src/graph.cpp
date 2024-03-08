#include "centrolign/graph.hpp"

namespace centrolign {

using namespace std;

uint64_t SequenceGraph::add_node(const string& sequence) {
    nodes.emplace_back(sequence);
    return nodes.size() - 1;
}

void SequenceGraph::add_edge(uint64_t node_id_from, uint64_t node_id_to) {
    nodes[node_id_from].next.push_back(node_id_to);
    nodes[node_id_to].prev.push_back(node_id_from);
}

void SequenceGraph::remove_edge(uint64_t node_id_from, uint64_t node_id_to) {
    auto& edges_from = nodes[node_id_from].next;
    auto& edges_to = nodes[node_id_to].prev;
    for (size_t i = 0; i < edges_from.size(); ++i) {
        if (edges_from[i] == node_id_to) {
            edges_from[i] = edges_from.back();
            edges_from.pop_back();
            break;
        }
    }
    for (size_t i = 0; i < edges_to.size(); ++i) {
        if (edges_to[i] == node_id_from) {
            edges_to[i] = edges_to.back();
            edges_to.pop_back();
            break;
        }
    }
}

void SequenceGraph::relabel(uint64_t node_id, const std::string& sequence) {
    nodes[node_id].seq = sequence;
}

string SequenceGraph::label(uint64_t node_id) const {
    return nodes[node_id].seq;
}

size_t SequenceGraph::label_size(uint64_t node_id) const {
    return nodes[node_id].seq.size();
}

size_t SequenceGraph::node_size() const {
    return nodes.size();
}

const vector<uint64_t>& SequenceGraph::next(uint64_t node_id) const {
    return nodes[node_id].next;
}

const vector<uint64_t>& SequenceGraph::previous(uint64_t node_id) const {
    return nodes[node_id].prev;
}

size_t SequenceGraph::next_size(uint64_t node_id) const {
    return nodes[node_id].next.size();
}

size_t SequenceGraph::previous_size(uint64_t node_id) const {
    return nodes[node_id].prev.size();
}

size_t SequenceGraph::path_size() const {
    return paths.size();
}

const std::string& SequenceGraph::path_name(uint64_t path_id) const {
    return paths[path_id].first;
}

const std::vector<uint64_t>& SequenceGraph::path(uint64_t path_id) const {
    return paths[path_id].second;
}

uint64_t SequenceGraph::path_id(const std::string& name) const {
    return name_to_id.at(name);
}

uint64_t SequenceGraph::add_path(const std::string& name) {
    name_to_id[name] = paths.size();
    paths.emplace_back(name, std::vector<uint64_t>());
    return paths.size() - 1;
}

void SequenceGraph::extend_path(uint64_t path_id, uint64_t node_id) {
    paths[path_id].second.push_back(node_id);
}

void SequenceGraph::pre_extend_path(uint64_t path_id, uint64_t node_id) {
    auto& path = paths[path_id].second;
    path.insert(path.begin(), node_id);
}

BaseGraphOverlay::BaseGraphOverlay(const SequenceGraph* sequence_graph) noexcept : seq_graph(sequence_graph) {
    
    cumul_len.reserve(sequence_graph->node_size() + 1);
    origin.reserve(sequence_graph->node_size());
    cumul_len.push_back(0);
    for (size_t i = 0, n = sequence_graph->node_size(); i < n; ++i) {
        size_t node_len = sequence_graph->label(i).size();
        cumul_len.push_back(cumul_len.back() + node_len);
        for (size_t j = 0; j < node_len; ++j) {
            origin.push_back(i);
        }
    }
}

size_t BaseGraphOverlay::node_size() const {
    return origin.size();
}

char BaseGraphOverlay::label(uint64_t node_id) const {
    uint64_t node_origin = origin[node_id];
    return seq_graph->label(node_origin)[node_id - cumul_len[node_origin]];
}

size_t BaseGraphOverlay::label_size(uint64_t node_id) const {
    return 1;
}

vector<uint64_t> BaseGraphOverlay::next(uint64_t node_id) const {
    vector<uint64_t> return_val;
    uint64_t node_origin = origin[node_id];
    if (node_id + 1 < cumul_len[node_origin + 1]) {
        return_val.push_back(node_id + 1);
    }
    else {
        for (auto next_id : seq_graph->next(node_origin)) {
            return_val.push_back(cumul_len[next_id]);
        }
    }
    return return_val;
}

vector<uint64_t> BaseGraphOverlay::previous(uint64_t node_id) const {
    vector<uint64_t> return_val;
    uint64_t node_origin = origin[node_id];
    if (node_id > cumul_len[node_origin]) {
        return_val.push_back(node_id - 1);
    }
    else {
        for (auto prev_id : seq_graph->previous(node_origin)) {
            return_val.push_back(cumul_len[prev_id + 1] - 1);
        }
    }
    return return_val;
}

size_t BaseGraphOverlay::next_size(uint64_t node_id) const {
    uint64_t node_origin = origin[node_id];
    if (node_id + 1 < cumul_len[node_origin + 1]) {
        return 1;
    }
    else {
        return seq_graph->next_size(node_origin);
    }
}

size_t BaseGraphOverlay::previous_size(uint64_t node_id) const {
    uint64_t node_origin = origin[node_id];
    if (node_id > cumul_len[node_origin]) {
        return 1;
    }
    else {
        return seq_graph->previous_size(node_origin);
    }
}

size_t BaseGraphOverlay::path_size() const {
    return seq_graph->path_size();
}

const string& BaseGraphOverlay::path_name(uint64_t path_id) const {
    return seq_graph->path_name(path_id);
}

vector<uint64_t> BaseGraphOverlay::path(uint64_t path_id) const {
    vector<uint64_t> path;
    for (uint64_t node_id : seq_graph->path(path_id)) {
        for (uint64_t base_id = cumul_len[node_id], n = cumul_len[node_id + 1]; base_id < n; ++base_id) {
            path.push_back(base_id);
        }
    }
    return path;
}

uint64_t BaseGraphOverlay::path_id(const std::string& name) const {
    return seq_graph->path_id(name);
}

uint64_t BaseGraph::add_node(char base) {
    nodes.emplace_back(base);
    return nodes.size() - 1;
}

void BaseGraph::add_edge(uint64_t node_id_from, uint64_t node_id_to) {
    nodes[node_id_from].next.push_back(node_id_to);
    nodes[node_id_to].prev.push_back(node_id_from);
}

void BaseGraph::remove_edge(uint64_t node_id_from, uint64_t node_id_to) {
    auto& edges_from = nodes[node_id_from].next;
    auto& edges_to = nodes[node_id_to].prev;
    for (size_t i = 0; i < edges_from.size(); ++i) {
        if (edges_from[i] == node_id_to) {
            edges_from[i] = edges_from.back();
            edges_from.pop_back();
            break;
        }
    }
    for (size_t i = 0; i < edges_to.size(); ++i) {
        if (edges_to[i] == node_id_from) {
            edges_to[i] = edges_to.back();
            edges_to.pop_back();
            break;
        }
    }
}

void BaseGraph::relabel(uint64_t node_id, char base) {
    nodes[node_id].base = base;
}

char BaseGraph::label(uint64_t node_id) const {
    return nodes[node_id].base;
}

size_t BaseGraph::label_size(uint64_t node_id) const {
    return 1;
}

size_t BaseGraph::node_size() const {
    return nodes.size();
}

const vector<uint64_t>& BaseGraph::next(uint64_t node_id) const {
    return nodes[node_id].next;
}

const vector<uint64_t>& BaseGraph::previous(uint64_t node_id) const {
    return nodes[node_id].prev;
}

size_t BaseGraph::next_size(uint64_t node_id) const {
    return nodes[node_id].next.size();
}

size_t BaseGraph::previous_size(uint64_t node_id) const {
    return nodes[node_id].prev.size();
}

// TODO: path interface is entirely repeated from SequenceGraph
size_t BaseGraph::path_size() const {
    return paths.size();
}

const std::string& BaseGraph::path_name(uint64_t path_id) const {
    return paths[path_id].first;
}

const std::vector<uint64_t>& BaseGraph::path(uint64_t path_id) const {
    return paths[path_id].second;
}

uint64_t BaseGraph::path_id(const std::string& name) const {
    return name_to_id.at(name);
}

uint64_t BaseGraph::add_path(const std::string& name) {
    name_to_id[name] = paths.size();
    paths.emplace_back(name, std::vector<uint64_t>());
    return paths.size() - 1;
}

void BaseGraph::extend_path(uint64_t path_id, uint64_t node_id) {
    paths[path_id].second.push_back(node_id);
}

void BaseGraph::pre_extend_path(uint64_t path_id, uint64_t node_id) {
    auto& path = paths[path_id].second;
    path.insert(path.begin(), node_id);
}

size_t BaseGraph::memory_size() const {
    size_t size = 0;
    size += nodes.capacity() * sizeof(BaseGraphNode);
    for (const auto& node : nodes) {
        size += (node.next.size() + node.prev.size()) * sizeof(uint64_t);
    }
    // TODO: i'm not fully accounting for the memory usage here, but it's probably okay
    // estimating as if each bucket is a singly-linked list
    size += name_to_id.bucket_count() * sizeof(void*);
    size += name_to_id.size() * (sizeof(uint64_t) + sizeof(string) + sizeof(void*));
    for (const auto& rec : name_to_id) {
        size += rec.first.capacity();
    }
    size += paths.capacity() * sizeof(std::pair<std::string, std::vector<uint64_t>>);
    for (const auto& path_rec : paths) {
        size += path_rec.first.capacity();
        size += path_rec.second.capacity() * sizeof(uint64_t);
    }
    return size;
}

}
