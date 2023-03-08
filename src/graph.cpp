#include "graph.hpp"

namespace centrolign {

using namespace std;

uint64_t SequenceGraph::add_node(const string& sequence) {
    nodes.emplace(sequence);
    return nodes.size() - 1;
}

void SequenceGraph::add_edge(uint64_t node_id_from, uint64_t node_id_to) {
    nodes[node_id_from].next.push_back(node_id_to);
    nodes[node_id_to].prev.push_back(node_id_from);
}

string SequenceGraph::sequence(uint64_t node_id) const {
    return nodes[node_id].seq;
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

BaseGraphOverlay::BaseGraphOverlay(const SequenceGraph* sequence_graph) : seq_graph(sequence_graph) {
    
    cumul_len.reserve(sequence_graph->node_size() + 1);
    origin.reserve(sequence_graph->node_size());
    cumul_len.push_back(0);
    for (size_t i = 0, n = sequence_graph->node_size(); i < n; ++i) {
        size_t node_len = sequence_graph->sequence(i).size();
        cumul_len.push_back(cumul_len.back() + node_len);
        for (size_t j = 0; j < node_len; ++j) {
            origin.push_back(i);
        }
    }
}

size_t BaseGraphOverlay::node_size() const {
    return origin.size();
}

char BaseGraphOverlay::base(uint64_t node_id) const {
    uint64_t node_origin = origin[node_id];
    return seq_graph->sequence(node_origin)[node_id - cumul_len[node_origin]];
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
        for (auto prev_id : seq_graph->prev(node_origin)) {
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

uint64_t BaseGraph::add_node(char base) {
    nodes.emplace(base);
    return nodes.size() - 1;
}

void BaseGraph::add_edge(uint64_t node_id_from, uint64_t node_id_to) {
    nodes[node_id_from].next.push_back(node_id_to);
    nodes[node_id_to].prev.push_back(node_id_from);
}

char BaseGraph::base(uint64_t node_id) const {
    return nodes[node_id].base;
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

}
