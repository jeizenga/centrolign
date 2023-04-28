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
    
    // set up the label query interface
    for (uint64_t node_id = 0; node_id < nodes.size(); ++node_id) {
        const auto& node = nodes[node_id];
        if (find(node.label.begin(), node.label.end(), '#') != node.label.end()) {
            throw runtime_error("Tree labels may not include '#': " + node.label + "\n");
        }
        
        if (!node.label.empty()) {
            if (label_map.count(node.label)) {
                throw runtime_error("Duplicate label " + node.label + " in Newick tree");
            }
            label_map[node.label] = node_id;
        }
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

uint64_t Tree::get_id(const std::string& label) const {
    return label_map.at(label);
}

uint64_t Tree::get_root() const {
    return root;
}

uint64_t Tree::get_parent(uint64_t node_id) const {
    return nodes[node_id].parent;
}

const std::vector<uint64_t>& Tree::get_children(uint64_t node_id) const {
    return nodes[node_id].children;
}

const std::string& Tree::label(uint64_t node_id) const {
    return nodes[node_id].label;
}

void Tree::binarize() {
    for (uint64_t node_id = 0, end = nodes.size(); node_id < end; ++node_id) {
        if (nodes[node_id].children.size() > 2) {
            
            
            string label = nodes[node_id].label;
            int label_num = 0;
            nodes[node_id].label = label + "#" + to_string(label_num++);
            
            // steal the children away so we can modify them
            auto children = move(nodes[node_id].children);
            
            // add the initial left child
            nodes[node_id].children.clear();
            nodes[node_id].children.push_back(children.front());
            
            uint64_t prev_node_id = node_id;
            
            for (size_t i = 2; i < children.size(); ++i) {
                
                // make a new node
                uint64_t new_node_id = nodes.size();
                nodes.emplace_back();
                auto& node = nodes[new_node_id];
                node.label = label + "#" + to_string(label_num++);
                node.distance = nodes[prev_node_id].distance;
               
                // set it as the right child of the previous
                node.parent = prev_node_id;
                nodes[prev_node_id].children.push_back(new_node_id);
                
                // give it a left children from the original children list
                node.children.push_back(children[i - 1]);
                nodes[node.children.front()].parent = nodes.size() - 1;
                
                prev_node_id = new_node_id;
            }
            
            // add the final right child
            nodes.back().children.push_back(children.back());
            nodes[children.back()].parent = prev_node_id;
        }
    }
}

std::vector<uint64_t> Tree::postorder() const {
    
    // bool indicates whether its edges have been followed
    vector<pair<uint64_t, bool>> stack;
    
    vector<uint64_t> order;
    order.reserve(nodes.size());
    
    if (root != -1) {
        stack.emplace_back(root, false);
    }
    
    while (!stack.empty()) {
        auto& top = stack.back();
        if (top.second) {
            // we've cleared the subtrees below this node
            order.push_back(top.first);
            stack.pop_back();
        }
        else {
            // queue up this nodes edges
            top.second = true;
            for (auto child : get_children(top.first)) {
                stack.emplace_back(child, false);
            }
        }
    }
    
    return order;
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
