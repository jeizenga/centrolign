#include "centrolign/tree.hpp"

#include <cctype>
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <stdexcept>
#include <sstream>
#include <functional>

namespace centrolign {

using namespace std;

const bool Tree::debug = false;

Tree::Tree(const string& newick) {
    
    {
        // check basic formatting validity
        auto it = find(newick.begin(), newick.end(), ';');
        if (it == newick.end()) {
            throw runtime_error("Newick string is missing a terminating ';'");
        }
        ++it;
        for (; it != newick.end(); ++it) {
            if (!isspace(*it)) {
                throw runtime_error("Newick string includes characters after the terminating ';'");
            }
        }
        
        uint64_t num_quotes = count(newick.begin(), newick.end(), '"');
        if (num_quotes % 2 == 1) {
            throw runtime_error("Newick string has an odd number of quotation marks");
        }
        if (find(newick.begin(), newick.end(), '\'') != newick.end()) {
            throw runtime_error("Newick string parser does not support single quotes (')");
        }
    }
    
    vector<uint64_t> stack;
    auto cursor = newick.begin();
    uint64_t ascending_node = -1;
    
    while (cursor < newick.end()) {
        
        // the characters that could indicate the boundary of a node
        auto next = find_skipping_quotes(cursor, newick.end(), ",();");
        
        if (*next == ';') {
            if (ascending_node != -1) {
                parse_label(ascending_node, cursor, next);
            }
            // we've already checked that there's nothing else in the string
            break;
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
        else {
            // this should never happen
            assert(false);
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
                throw runtime_error("Duplicate label " + node.label + " in guide tree");
            }
            label_map[node.label] = node_id;
        }
    }
}

std::string::const_iterator Tree::find_skipping_quotes(std::string::const_iterator begin,
                                                       std::string::const_iterator end,
                                                       const std::string& values) const {
        
    bool in_quote = false;
    for (auto it = begin; it != end; ++it) {
        if (*it == '"') {
            // entering or leaving quotation-delimited substring
            in_quote = !in_quote;
        }
        else if (!in_quote) {
            for (auto v_it = values.begin(); v_it != values.end(); ++v_it) {
                if (*it == *v_it) {
                    return it;
                }
            }
        }
    }
    return end;
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
    auto div = find_skipping_quotes(begin, end, ":");
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
    for (auto c = b + 1; c + 1 < e; ++c) {
        if (*c == '"') {
            throw runtime_error("Newick string label has internal quotation mark: " + string(b, e));
        }
    }
    // handle quotation marks
    if (b < e && *b == '"') {
        if (b + 1 == e) {
            throw runtime_error("Newick string label consists of only one quotation mark: " + string(b, e));
        }
        if (*(e - 1) != '"') {
            throw runtime_error("Newick string label has unmatched quotation mark: " + string(b, e));
        }
        // move past the quotation marks
        ++b;
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

size_t Tree::node_size() const {
    return nodes.size();
}

bool Tree::has_label(const std::string& label) const {
    return label_map.count(label);
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

double Tree::distance(uint64_t node_id) const {
    return nodes[node_id].distance;
}

bool Tree::is_leaf(uint64_t node_id) const {
    return nodes[node_id].children.empty();
}

void Tree::binarize() {
    for (uint64_t node_id = 0, end = nodes.size(); node_id < end; ++node_id) {
        if (nodes[node_id].children.size() > 2) {

            string label = std::move(nodes[node_id].label);
            int label_num = 0;
            if (!label.empty()) {
                nodes[node_id].label = label + "#" + to_string(label_num++);
            }
            
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
                if (!label.empty()) {
                    node.label = label + "#" + to_string(label_num++);
                }
                // this node is just a stand-in for the original node, so it occurs
                // at no additional distance
                node.distance = 0.0;
               
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

void Tree::prune(const std::vector<uint64_t>& node_ids) {
    
    // mark the nodes that we're going to keep
    std::vector<bool> keep(nodes.size(), false);
    for (auto node_id : node_ids) {
        auto here = node_id;
        while (here != -1 && !keep[here]) {
            keep[here] = true;
            here = get_parent(here);
        }
    }
    
    auto keep_children = [&](uint64_t node_id) {
        std::vector<uint64_t> kc;
        for (auto c : nodes[node_id].children) {
            if (keep[c]) {
                kc.push_back(c);
            }
        }
        return kc;
    };
    
    // walk down to find the LCA of the keep nodes
    uint64_t here = root;
    while (here != -1 && keep[here] && keep_children(here).size() == 1) {
        keep[here] = false;
        here = keep_children(here).front();
    }
    
    // handle the case of a single node, where the LCA algorithm breaks down
    if (!node_ids.empty()) {
        keep[node_ids.front()] = true;
    }
    
    filter(keep);
    // make sure the root doesn't have a distance from parent
    if (!nodes.empty()) {
        nodes[root].distance = std::numeric_limits<double>::max();
    }
}

void Tree::compact() {
    
    std::vector<bool> keep(nodes.size(), true);
    
    for (uint64_t node_id = 0; node_id < nodes.size(); ++node_id) {
        auto& node = nodes[node_id];
        if (node.children.size() == 1) {
            // this is a non-branching path
            keep[node_id] = false;
            // rewire the edges around this node
            if (node_id == root) {
                root = node.children.front();
                nodes[node.children.front()].parent = -1;
            }
            else {
                nodes[node.parent].children.push_back(node.children.front());
                nodes[node.children.front()].parent = node.parent;
            }
        }
    }
    
    // add the distances onto the bottom of each compacted path
    for (uint64_t node_id = 0; node_id < nodes.size(); ++node_id) {
        if (keep[node_id]) {
            
            auto& node = nodes[node_id];
            
            auto here = node.parent;
            while (here != -1 && !keep[here] &&
                   node.distance != std::numeric_limits<double>::max()) {
                if (nodes[here].distance != std::numeric_limits<double>::max()) {
                    node.distance += nodes[here].distance;
                }
                else {
                    node.distance = std::numeric_limits<double>::max();
                }
            }
        }
    }
    
    filter(keep);
}

void Tree::filter(const std::vector<bool>& keep) {
    
    // move keep nodes into a prefix
    std::vector<size_t> removed_so_far(nodes.size() + 1, 0);
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (!keep[i]) {
            if (!nodes[i].label.empty()) {
                label_map.erase(nodes[i].label);
            }
            removed_so_far[i + 1] = removed_so_far[i] + 1;
        }
        else {
            if (removed_so_far[i]) {
                nodes[i - removed_so_far[i]] = move(nodes[i]);
                if (!nodes[i].label.empty()) {
                    label_map[nodes[i].label] = i - removed_so_far[i];
                }
            }
            removed_so_far[i + 1] = removed_so_far[i];
        }
    }
    
    if (removed_so_far.back()) {
        // remove the positions that we moved up
        nodes.resize(nodes.size() - removed_so_far.back());
        
        // update child indexes
        for (auto& node : nodes) {
            size_t children_removed = 0;
            for (size_t i = 0; i < node.children.size(); ++i) {
                auto child = node.children[i];
                if (!keep[child]) {
                    ++children_removed;
                }
                else {
                    node.children[i - children_removed] = child - removed_so_far[child];
                }
            }
            node.children.resize(node.children.size() - children_removed);
            // erase the parent pointer
            node.parent = -1;
        }
        
        // reset the parent pointers
        for (uint64_t node_id = 0; node_id < nodes.size(); ++node_id) {
            for (auto child_id : nodes[node_id].children) {
                nodes[child_id].parent = node_id;
            }
        }
        
        // find the new root
        for (uint64_t node_id = 0; node_id < nodes.size(); ++node_id) {
            if (nodes[node_id].parent == -1) {
                root = node_id;
                break;
            }
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

std::vector<uint64_t> Tree::preorder() const {
    
    vector<uint64_t> stack;
    
    vector<uint64_t> order;
    order.reserve(nodes.size());
    
    if (root != -1) {
        stack.emplace_back(root);
    }
    
    while (!stack.empty()) {
        auto top = stack.back();
        stack.pop_back();
        order.push_back(top);
        for (auto child : get_children(top)) {
            stack.emplace_back(child);
        }
    }
    
    return order;
}

std::string Tree::to_newick() const {
    
    std::stringstream strm;
    
    // recursive internal function
    function<void(uint64_t)> recurse = [&](uint64_t node_id) {
        const auto& node = nodes[node_id];
        if (!node.children.empty()) {
            strm << '(';
            for (size_t i = 0; i < node.children.size(); ++i) {
                if (i != 0) {
                    strm << ',';
                }
                recurse(node.children[i]);
            }
            strm << ')';
        }
        if (!node.label.empty()) {
            strm << '"' << node.label << '"';
        }
        if (node.distance != std::numeric_limits<double>::max()) {
            strm << ':' << node.distance;
        }
    };
    
    if (root != -1) {
        recurse(root);
    }
    
    strm << ';';
    
    return strm.str();
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
