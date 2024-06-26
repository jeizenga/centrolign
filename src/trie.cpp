#include "centrolign/trie.hpp"

#include <cctype>
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <stdexcept>

namespace centrolign {

using namespace std;

Trie::Trie() {
    // init the root
    nodes.emplace_back();
}

uint64_t Trie::insert_sequence(const std::string& name,
                               const std::vector<uint64_t>& sequence) {
    
    paths.emplace_back();
    paths.back().first = name;
    auto& path = paths.back().second;
    
    uint64_t here = get_root();
    for (size_t i = 0; i < sequence.size(); ++i) {
        uint64_t next = follow(here, sequence[i]);
        if (next != -1) {
            path.push_back(next);
            here = next;
        }
        else {
            nodes[here].children[sequence[i]] = nodes.size();
            path.push_back(nodes.size());
            
            nodes.emplace_back();
            auto& node = nodes.back();
            node.parent = here;
            node.label = sequence[i];
            
            here = nodes.size() - 1;
        }
    }
    
    return paths.size() - 1;
}

void Trie::clear() {
    nodes.clear();
    paths.clear();
}

size_t Trie::node_size() const {
    return nodes.size();
}

uint64_t Trie::get_root() const {
    return 0;
}

uint64_t Trie::get_parent(uint64_t node_id) const {
    return nodes[node_id].parent;
}

std::vector<uint64_t> Trie::next(uint64_t node_id) const {
    std::vector<uint64_t> children;
    for (const auto& edge : nodes[node_id].children) {
        children.push_back(edge.second);
    }
    return children;
}

size_t Trie::next_size(uint64_t node_id) const {
    return nodes[node_id].children.size();
}

std::vector<uint64_t> Trie::previous(uint64_t node_id) const {
    std::vector<uint64_t> prev;
    if (node_id != get_root()) {
        prev.push_back(nodes[node_id].parent);
    }
    return prev;
}

size_t Trie::previous_size(uint64_t node_id) const {
    return node_id == get_root() ? 0 : 1;
}

uint64_t Trie::label(uint64_t node_id) const {
    return nodes[node_id].label;
}

uint64_t Trie::follow(uint64_t node_id, uint64_t label) const {
    auto it = nodes[node_id].children.find(label);
    if (it != nodes[node_id].children.end()) {
        return it->second;
    }
    else {
        return -1;
    }
}

size_t Trie::path_size() const {
    return paths.size();
}

const std::string& Trie::path_name(uint64_t path_id) const {
    return paths[path_id].first;
}

const std::vector<uint64_t>& Trie::path(uint64_t path_id) const {
    return paths[path_id].second;
}



}
