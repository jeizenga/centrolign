#ifndef centrolign_max_search_tree_hpp
#define centrolign_max_search_tree_hpp

#include <vector>
#include <cstdint>
#include <iostream>

namespace centrolign {

/*
 * Static-topology search tree that supports range max queries and value updates
 */
template<typename K, typename V>
class MaxSearchTree {
public:
    
    // vector will be sorted
    MaxSearchTree(std::vector<std::pair<K, V>>& values);
    MaxSearchTree() = default;
    ~MaxSearchTree() = default;
    
    class iterator; // forward declaration
    
    // standard iteration interface
    iterator begin() const;
    iterator end() const;
    // the range of iterators that have this key
    std::pair<iterator, iterator> equal_range(const K& key) const;
    // an arbitrary iterator that has this key
    iterator find(const K& key) const;
    // change value at the key the iterator points to
    void update(const iterator& it, const V& value);
    // returns iterator to max value in key range [lo, hi)
    iterator range_max(const K& lo, const K& hi) const;
    
    class iterator {
    public:
        iterator() = default;
        ~iterator() = default;
        
        iterator& operator++();
        const std::pair<K, V>& operator*() const;
        const std::pair<K, V>* operator->() const;
        bool operator==(const iterator& other) const;
        bool operator!=(const iterator& other) const;
        
        
    private:
        friend class MaxSearchTree<K, V>;
        
        // internal constructor
        iterator(const MaxSearchTree<K, V>& iteratee, size_t i);
        // tree we're iterating over
        const MaxSearchTree<K, V>* iteratee = nullptr;
        // index of the node
        size_t i = 0;
    };
    
private:
    
    static const bool debug_mst = false;
    
    struct Node {
        Node() = default;
        Node(const std::pair<K, V>& key_value, size_t idx) : key_value(key_value), subtree_max(idx) { }
        std::pair<K, V> key_value;
        size_t subtree_max = 0;
    };
    
    void reidentify_subtree_max(size_t x);
    
    static inline size_t left(size_t x);
    static inline size_t right(size_t x);
    static inline size_t parent(size_t x);
    
    std::vector<Node> nodes;
    
    friend class iterator;
};

/*
 * Template implementations
 */

template<typename K, typename V>
MaxSearchTree<K, V>::MaxSearchTree(std::vector<std::pair<K, V>>& values) {
    
    // handle this as a special case so we can
    if (values.empty()) {
        return;
    }
    
    // comparing only on key, not value
    auto cmp = [](const std::pair<K, V>& a, const std::pair<K, V>& b) {
        return a.first < b.first;
    };
    if (!std::is_sorted(values.begin(), values.end(), cmp)) {
        std::stable_sort(values.begin(), values.end(), cmp);
    }
    
    // figure out the height of the tree
    size_t height = 0;
    while ((1 << height) - 1 < values.size()) { // max capacity of a tree with height
        ++height;
    }
    
    if (debug_mst) {
        std::cerr << "building tree of height " << height << " for " << values.size() << " key-value pairs\n";
    }
    
    // init the implicit tree
    nodes.resize(values.size());
    
    // the next item in the vector we will assign to a node
    size_t vec_idx = 0;
    // records of (node_idx, queued left)
    std::vector<std::pair<size_t, bool>> stack;
    stack.emplace_back(0, false);
    
    // in-order traversal of the tree through a stack
    while (!stack.empty()) {
        auto& top = stack.back();
        if (!top.second) {
            // we haven't traversed the left branch yet
            top.second = true;
            size_t l = left(top.first);
            if (l < nodes.size()) {
                stack.emplace_back(l, false);
            }
        }
        else {
            // the left is traversed, this is the node's position
            size_t node_idx = top.first;
            nodes[node_idx] = Node(values[vec_idx++], node_idx);
            stack.pop_back();
            // queue up the right branch as well
            size_t r = right(node_idx);
            if (r < nodes.size()) {
                stack.emplace_back(r, false);
            }
        }
    }
    
    if (debug_mst) {
        std::cerr << "identifying subtree maxima\n";
    }
    
    // all nodes are now initialized to point to themselves for subtree max, we have to
    // work upward through the tree to point at the correct location
    
    for (size_t i = nodes.size() - 1; i > 0; --i) {
        const auto& node = nodes[i];
        auto& par = nodes[parent(i)];
        if (nodes[node.subtree_max].key_value.second > nodes[par.subtree_max].key_value.second) {
            par.subtree_max = node.subtree_max;
        }
    }
    if (debug_mst) {
        std::cerr << "finished constructing tree:\n";
        for (size_t i = 0; i < nodes.size(); ++i) {
            auto& node = nodes[i];
//            std::cerr << i << ", key " << node.key_value.first;
//            std::cerr << ", val " << node.key_value.second.first << ',' << node.key_value.second.second;
            std::cerr << ", max " << node.subtree_max;
            if (left(i) < nodes.size()) {
                std::cerr << ", left " << left(i);
            }
            if (right(i) < nodes.size()) {
                std::cerr << ", right " << right(i);
            }
            std::cerr << '\n';
        }
    }
}

template<typename K, typename V>
inline size_t MaxSearchTree<K, V>::left(size_t x) {
    return 2 * x + 1;
}

template<typename K, typename V>
inline size_t MaxSearchTree<K, V>::right(size_t x) {
    return 2 * x + 2;
}

template<typename K, typename V>
inline size_t MaxSearchTree<K, V>::parent(size_t x) {
    return (x - 1) / 2;
}


template<typename K, typename V>
typename MaxSearchTree<K, V>::iterator MaxSearchTree<K, V>::begin() const {
    size_t i = 0;
    size_t l = left(i);
    while (l < nodes.size()) {
        i = l;
        l = left(i);
    }
    return iterator(*this, i);
}

template<typename K, typename V>
typename MaxSearchTree<K, V>::iterator MaxSearchTree<K, V>::end() const {
    return iterator(*this, nodes.size());
}


template<typename K, typename V>
typename MaxSearchTree<K, V>::iterator MaxSearchTree<K, V>::find(const K& key) const {
    size_t cursor = 0;
    while (cursor < nodes.size()) {
        if (nodes[cursor].key_value.first == key) {
            return iterator(*this, cursor);
        }
        else if (nodes[cursor].key_value.first > key) {
            cursor = left(cursor);
        }
        else {
            cursor = right(cursor);
        }
    }
    return end();
}

template<typename K, typename V>
std::pair<typename MaxSearchTree<K, V>::iterator, typename MaxSearchTree<K, V>::iterator>
MaxSearchTree<K, V>::equal_range(const K& key) const {
    
    if (debug_mst) {
        std::cerr << "finding equal range for key ";
//        std::cerr << key;
        std::cerr << '\n';
    }
    
    size_t lower = -1, upper = -1;
    
    // find the lowest node that is inside the range
    size_t cursor = 0;
    while (cursor < nodes.size()) {
        if (nodes[cursor].key_value.first == key) {
            if (debug_mst) {
                std::cerr << "found new lower bound " << cursor << '\n';
            }
            lower = cursor;
            cursor = left(cursor);
        }
        else if (nodes[cursor].key_value.first > key) {
            cursor = left(cursor);
        }
        else {
            cursor = right(cursor);
        }
        if (debug_mst) {
            std::cerr << "advance left cursor to " << cursor << '\n';
        }
    }
    
    // find the lowest node that is outside the range
    cursor = 0;
    while (cursor < nodes.size()) {
        if (nodes[cursor].key_value.first > key) {
            if (debug_mst) {
                std::cerr << "found new upper bound " << cursor << '\n';
            }
            upper = cursor;
            cursor = left(cursor);
        }
        else {
            cursor = right(cursor);
        }
        if (debug_mst) {
            std::cerr << "advance right cursor to " << cursor << '\n';
        }
    }
    
    if (lower != -1 && upper != -1) {
        // we found bounding nodes
        return std::make_pair(iterator(*this, lower), iterator(*this, upper));
    }
    else if (lower != -1) {
        // we found a lower bound node (so the range is not empty), but no upper bound
        // so it must go all the way through the end
        return std::make_pair(iterator(*this, lower), end());
    }
    else {
        return std::make_pair(end(), end());
    }
    
}

template<typename K, typename V>
void MaxSearchTree<K, V>::reidentify_subtree_max(size_t x) {
    size_t new_max = x;
    size_t l = left(x);
    if (l < nodes.size() && nodes[l].key_value.second > nodes[new_max].key_value.second) {
        new_max = l;
    }
    size_t r = right(x);
    if (r < nodes.size() && nodes[r].key_value.second > nodes[new_max].key_value.second) {
        new_max = r;
    }
    nodes[x].subtree_max = new_max;
}

template<typename K, typename V>
void MaxSearchTree<K, V>::update(const MaxSearchTree<K, V>::iterator& it, const V& value) {
    
    if (debug_mst) {
        std::cerr << "performing update\n";
    }
    
    // TODO: these could be simplified into only the second case, but the first one seems
    // more likely to me, and it is a simpler loop...
    
    auto& node = nodes[it.i];
    if (value > nodes[node.subtree_max].key_value.second) {
        // we're making a new max for the subtree rooted here
        node.subtree_max = it.i;
        
        // this node might become the new max on the path to the root
        size_t here = it.i;
        while (here != 0) {
            here = parent(here);
            if (value > nodes[nodes[here].subtree_max].key_value.second) {
                nodes[here].subtree_max = it.i;
            }
            else {
                break;
            }
        }
        node.key_value.second = value;
    }
    else {
        // this value is less than the max
        node.key_value.second = value;
        if (node.subtree_max == it.i) {
            // the max was formerly at this node, but it might not be now
            reidentify_subtree_max(it.i);
            
            // we may also need to update the max toward the root
            size_t here = it.i;
            while (here != 0) {
                here = parent(here);
                if (nodes[here].subtree_max != it.i) {
                    // the other subtree has a larger value anyway, we can stop early
                    break;
                }
                reidentify_subtree_max(here);
            }
        }
    }
}

template<typename K, typename V>
typename MaxSearchTree<K, V>::iterator MaxSearchTree<K, V>::range_max(const K& lo, const K& hi) const {
    
    if (debug_mst) {
        std::cerr << "searching for range max\n";
    }
    
    // traverse downward until finding the first value in the interval
    size_t cursor = 0;
    while (cursor < nodes.size() &&
           (nodes[cursor].key_value.first < lo || nodes[cursor].key_value.first >= hi)) {
        if (nodes[cursor].key_value.first >= lo) {
            cursor = left(cursor);
        }
        else {
            cursor = right(cursor);
        }
        if (debug_mst) {
            std::cerr << "advance initial search cursor to " << cursor << '\n';
        }
    }
    if (debug_mst) {
        std::cerr << "initial search ends at " << cursor << '\n';
    }
    
    if (cursor >= nodes.size()) {
        // we never found a value that was in the interval
        if (debug_mst) {
            std::cerr << "range is empty\n";
        }
        return end();
    }
    
    size_t max_idx = cursor;
    
    size_t right_cursor = right(cursor);
    size_t left_cursor = left(cursor);
    if (debug_mst) {
        std::cerr << "start left search at " << left_cursor << '\n';
    }
    // split off the leftward traversal where rightward off-path edges are included entirely
    while (left_cursor < nodes.size()) {
        if (nodes[left_cursor].key_value.first >= lo) {
            if (nodes[left_cursor].key_value.second > nodes[max_idx].key_value.second) {
                max_idx = left_cursor;
            }
            size_t r = right(left_cursor);
            if (r < nodes.size() && nodes[nodes[r].subtree_max].key_value.second > nodes[max_idx].key_value.second) {
                max_idx = nodes[r].subtree_max;
            }
            left_cursor = left(left_cursor);
        }
        else {
            left_cursor = right(left_cursor);
        }
        if (debug_mst) {
            std::cerr << "advance left search cursor to " << left_cursor << '\n';
        }
    }
    if (debug_mst) {
        std::cerr << "start right search at " << left_cursor << '\n';
    }
    // and next the rightward traversal where leftward off-path edges are included entirely
    while (right_cursor < nodes.size()) {
        if (nodes[right_cursor].key_value.first < hi) {
            if (nodes[right_cursor].key_value.second > nodes[max_idx].key_value.second) {
                max_idx = right_cursor;
            }
            size_t l = left(right_cursor);
            if (l < nodes.size() && nodes[nodes[l].subtree_max].key_value.second > nodes[max_idx].key_value.second) {
                max_idx = nodes[l].subtree_max;
            }
            right_cursor = right(right_cursor);
        }
        else {
            right_cursor = left(right_cursor);
        }
        if (debug_mst) {
            std::cerr << "advance right search cursor to " << right_cursor << '\n';
        }
    }
    if (debug_mst) {
        std::cerr << "returning max at node index " << max_idx << '\n';
    }
    return iterator(*this, max_idx);
}

template<typename K, typename V>
MaxSearchTree<K, V>::iterator::iterator(const MaxSearchTree<K, V>& iteratee, size_t i) :
    iteratee(&iteratee), i(i)
{
    // nothing to do
}

template<typename K, typename V>
typename MaxSearchTree<K, V>::iterator& MaxSearchTree<K, V>::iterator::operator++() {
    size_t r = MaxSearchTree<K, V>::right(i);
    if (r < iteratee->nodes.size()) {
        // go to leftmost node of right subtree
        i = r;
        while (MaxSearchTree<K, V>::left(i) < iteratee->nodes.size()) {
            i = MaxSearchTree<K, V>::left(i);
        }
    }
    else if (i == 0) {
        // we are the root, and there is no right branch, this is the end
        i = iteratee->nodes.size();
    }
    else {
        // walk until we cross an edge rightward
        while (true) {
            size_t p = MaxSearchTree<K, V>::parent(i);
            if (i == MaxSearchTree<K, V>::left(p)) {
                // edge is rightward
                i = p;
                break;
            }
            if (p == 0) {
                // we walked only leftward edges to the root, we must have been at the end
                i = iteratee->nodes.size();
                break;
            }
            i = p;
        }
    }
    return *this;
}

template<typename K, typename V>
const std::pair<K, V>& MaxSearchTree<K, V>::iterator::operator*() const {
    return iteratee->nodes[i].key_value;
}

template<typename K, typename V>
const std::pair<K, V>* MaxSearchTree<K, V>::iterator::operator->() const {
    return &(iteratee->nodes[i].key_value);
}

template<typename K, typename V>
bool MaxSearchTree<K, V>::iterator::operator==(const iterator& other) const {
    return other.iteratee == this->iteratee && other.i == this->i;
}

template<typename K, typename V>
bool MaxSearchTree<K, V>::iterator::operator!=(const iterator& other) const {
    return !(*this == other);
}
    
}

#endif /* centrolign_max_search_tree_hpp */
