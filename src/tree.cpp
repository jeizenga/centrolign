#include "centrolign/tree.hpp"

#include <cctype>
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <stdexcept>

namespace centrolign {

using namespace std;

const bool Tree::debug = false;

Tree::Tree(const string& newick) {
    
    vector<uint64_t> stack;
    
    // the characters that could indicate the boundary of a node
    const static string terminators = ",();";
    
    auto cursor = newick.begin();
    uint64_t ascending_node = -1;;
    while (cursor < newick.end()) {
        
        auto next = find_first_of(cursor, newick.end(),
                                  terminators.begin(), terminators.end());
        
        if (*next == ';') {
            if (ascending_node != -1) {
                parse_label(ascending_node, cursor, next);
            }
        }
        else if (*next == '(') {
            // start a new node
            
            uint64_t node_id = 0;
            if (stack.empty()) {
                if (root != -1 || !nodes.empty()) {
                    throw runtime_error("Newick string encodes a disconnected tree");
                }
                nodes.emplace_back();
                root = node_id;
            }
            else {
                node_id = add_child(stack.back());
            }
            // this node has children
            stack.push_back(node_id);
            ascending_node = -1;
        }
        else if (*next == ',') {
            if (ascending_node == -1) {
                // leaf node
                uint64_t node_id = add_child(stack.back());
                parse_label(node_id, cursor, next);
            }
            else {
                parse_label(ascending_node, cursor, next);
            }
            ascending_node = -1;
        }
        else if (*next == ')') {
            if (ascending_node == -1) {
                // final leaf node among children
                uint64_t node_id = add_child(stack.back());
                parse_label(node_id, cursor, next);
            }
            else {
                parse_label(ascending_node, cursor, next);
            }
            // we will no longer be adding children to this node
            ascending_node = stack.back();
            stack.pop_back();
        }
        
        cursor = next + 1;
    }
    
    
    if (debug) {
        debug_print(cerr);
    }
}

uint64_t Tree::add_child(uint64_t parent_id) {
    
    uint64_t node_id = nodes.size();
    nodes.emplace_back();
    auto& node = nodes.back();
    
    auto& parent = nodes[parent_id];
    parent.children.push_back(node_id);
    node.parent = parent_id;
    
    return node_id;
}

void Tree::parse_label(uint64_t node_id,
                       std::string::const_iterator begin,
                       std::string::const_iterator end) {
    
    auto& node = nodes[node_id];
    
    // find divider for label and distance (or go to the end if no disatnce)
    auto div = find(begin, end, ':');
    auto b = begin;
    // skip leading whitespace
    while (b < div && isspace(*b)) {
        ++b;
    }
    // skip trailing whitespace
    auto e = div;
    while (e - 1 > b && isspace(*(e - 1))) {
        --e;
    }
    
    node.label = string(b, e);
    
    if (div != end) {
        auto b = div + 1;
        // skip leading whitespace
        while (b < end && isspace(*b)) {
            ++b;
        }
        
        if (b == end) {
            throw runtime_error("Newick string has ':' without a distance following it");
        }
        
        // parse the distance
        node.distance = strtod(&(*b), nullptr);
    }
}

void Tree::debug_print(ostream& out) const {
    if (root != -1) {
        debug_print_internal(root, 0, out);
    }
}

void Tree::debug_print_internal(uint64_t node, size_t level, ostream& out) const {
    const auto& n = nodes[node];
    for (size_t i = 0; i < level; ++i) {
        out << '\t';
    }
    out << "-> " << node;
    if (!n.label.empty()) {
        out << ' ' << n.label;
    }
    if (n.distance != numeric_limits<double>::max()) {
        out << ' ' << n.distance;
    }
    out << '\n';
    for (auto child : n.children) {
        debug_print_internal(child, level + 1, out);
    }
}
}
