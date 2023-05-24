#include "centrolign/compacted_graph.hpp"


namespace centrolign {

using namespace std;

size_t CompactedGraph::node_size() const {
    return nodes.size();
}

size_t CompactedGraph::label_size(uint64_t node_id) const {
    return nodes[node_id].size;
}

uint64_t CompactedGraph::front(uint64_t node_id) const {
    return nodes[node_id].front;
}

uint64_t CompactedGraph::back(uint64_t node_id) const {
    return nodes[node_id].back;
}

const std::vector<uint64_t>& CompactedGraph::next(uint64_t node_id) const {
    return nodes[node_id].next;
}

const std::vector<uint64_t>& CompactedGraph::previous(uint64_t node_id) const {
    return nodes[node_id].prev;
}

size_t CompactedGraph::next_size(uint64_t node_id) const {
    return nodes[node_id].next.size();
}

size_t CompactedGraph::previous_size(uint64_t node_id) const {
    return nodes[node_id].prev.size();
}



}
