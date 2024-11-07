#ifndef centrolign_max_search_tree_hpp
#define centrolign_max_search_tree_hpp

#include <vector>
#include <cstdint>
#include <iostream>

#include "centrolign/utility.hpp"

namespace centrolign {

/*
 * Static-topology search tree that supports range max queries and value updates
 */
template<typename K, typename V, class KeyVector = std::vector<K>, class ValueVector = std::vector<V>, class IndexVector = std::vector<size_t>>
class MaxSearchTree {
public:
    
    // vector will be sorted
    MaxSearchTree(std::vector<std::pair<K, V>>& data);
    MaxSearchTree() = default;
    MaxSearchTree(const MaxSearchTree<K, V, KeyVector, ValueVector, IndexVector>& other) noexcept = default;
    MaxSearchTree(MaxSearchTree<K, V, KeyVector, ValueVector, IndexVector>&& other) noexcept = default;
    ~MaxSearchTree() = default;
    MaxSearchTree& operator=(const MaxSearchTree<K, V, KeyVector, ValueVector, IndexVector>& other) noexcept = default;
    MaxSearchTree& operator=(MaxSearchTree<K, V, KeyVector, ValueVector, IndexVector>&& other) noexcept = default;
    
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

    size_t size() const;
    bool empty() const;
    
    class iterator {
    public:
        iterator() = default;
        ~iterator() = default;
        
        iterator& operator++();
        const std::pair<K, V> operator*() const;
        bool operator==(const iterator& other) const;
        bool operator!=(const iterator& other) const;
        
        
    private:
        friend class MaxSearchTree<K, V, KeyVector, ValueVector, IndexVector>;
        
        // internal constructor
        iterator(const MaxSearchTree<K, V, KeyVector, ValueVector, IndexVector>& iteratee, size_t i);
        // tree we're iterating over
        const MaxSearchTree<K, V, KeyVector, ValueVector, IndexVector>* iteratee = nullptr;
        // index of the node
        size_t i = 0;
    };
    
    size_t memory_size() const;
    
private:
    
    static const bool debug_mst = false;
    
    KeyVector key;
    ValueVector value;
    IndexVector subtree_max;
    
    void reidentify_subtree_max(size_t x);
    
    static inline size_t left(size_t x);
    static inline size_t right(size_t x);
    static inline size_t parent(size_t x);
        
    friend class iterator;
};

/*
 * Template implementations
 */

template<typename K, typename V, class KeyVector, class ValueVector, class IndexVector>
MaxSearchTree<K, V, KeyVector, ValueVector, IndexVector>::MaxSearchTree(std::vector<std::pair<K, V>>& data) : key(data.size()), value(data.size()), subtree_max(data.size()) {
    
    // handle this as a special case so we can
    if (data.empty()) {
        return;
    }
    
    // comparing only on key, not value
    auto cmp = [](const std::pair<K, V>& a, const std::pair<K, V>& b) {
        return a.first < b.first;
    };
    if (!std::is_sorted(data.begin(), data.end(), cmp)) {
        std::stable_sort(data.begin(), data.end(), cmp);
    }
    
    // figure out the height of the tree
    size_t height = 0;
    while ((1 << height) - 1 < data.size()) { // max capacity of a tree with height
        ++height;
    }
    
    if (debug_mst) {
        std::cerr << "building tree of height " << height << " for " << data.size() << " key-value pairs\n";
    }
    
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
            if (l < size()) {
                stack.emplace_back(l, false);
            }
        }
        else {
            // the left is traversed, this is the node's position
            key[top.first] = data[vec_idx].first;
            value[top.first] = data[vec_idx].second;
            subtree_max[top.first] = top.first;
            ++vec_idx;
            
            // queue up the right branch as well
            size_t r = right(top.first);
            stack.pop_back();
            if (r < size()) {
                stack.emplace_back(r, false);
            }
        }
    }
    
    if (debug_mst) {
        std::cerr << "identifying subtree maxima\n";
    }
    
    // all nodes are now initialized to point to themselves for subtree max, we have to
    // work upward through the tree to point at the correct location
    
    for (size_t i = size() - 1; i > 0; --i) {
        auto par = parent(i);
        if (value[subtree_max[i]] > value[subtree_max[par]]) {
            subtree_max[par] = subtree_max[i];
        }
    }
    if (debug_mst) {
        std::cerr << "finished constructing tree:\n";
        for (size_t i = 0; i < size(); ++i) {
//            std::cerr << i << ", key " << node.key_value.first;
//            std::cerr << ", val " << node.key_value.second.first << ',' << node.key_value.second.second;
            std::cerr << ", max " << subtree_max[i];
            if (left(i) < size()) {
                std::cerr << ", left " << left(i);
            }
            if (right(i) < size()) {
                std::cerr << ", right " << right(i);
            }
            std::cerr << '\n';
        }
    }
}

template<typename K, typename V, class KeyVector, class ValueVector, class IndexVector>
inline size_t MaxSearchTree<K, V, KeyVector, ValueVector, IndexVector>::left(size_t x) {
    return 2 * x + 1;
}

template<typename K, typename V, class KeyVector, class ValueVector, class IndexVector>
inline size_t MaxSearchTree<K, V, KeyVector, ValueVector, IndexVector>::right(size_t x) {
    return 2 * x + 2;
}

template<typename K, typename V, class KeyVector, class ValueVector, class IndexVector>
inline size_t MaxSearchTree<K, V, KeyVector, ValueVector, IndexVector>::parent(size_t x) {
    return (x - 1) / 2;
}


template<typename K, typename V, class KeyVector, class ValueVector, class IndexVector>
typename MaxSearchTree<K, V, KeyVector, ValueVector, IndexVector>::iterator MaxSearchTree<K, V, KeyVector, ValueVector, IndexVector>::begin() const {
    size_t i = 0;
    size_t l = left(i);
    while (l < size()) {
        i = l;
        l = left(i);
    }
    return iterator(*this, i);
}

template<typename K, typename V, class KeyVector, class ValueVector, class IndexVector>
typename MaxSearchTree<K, V, KeyVector, ValueVector, IndexVector>::iterator MaxSearchTree<K, V, KeyVector, ValueVector, IndexVector>::end() const {
    return iterator(*this, size());
}


template<typename K, typename V, class KeyVector, class ValueVector, class IndexVector>
typename MaxSearchTree<K, V, KeyVector, ValueVector, IndexVector>::iterator MaxSearchTree<K, V, KeyVector, ValueVector, IndexVector>::find(const K& search_key) const {
    size_t cursor = 0;
    while (cursor < size()) {
        if (key[cursor] == search_key) {
            return iterator(*this, cursor);
        }
        else if (key[cursor] > search_key) {
            cursor = left(cursor);
        }
        else {
            cursor = right(cursor);
        }
    }
    return end();
}

template<typename K, typename V, class KeyVector, class ValueVector, class IndexVector>
std::pair<typename MaxSearchTree<K, V, KeyVector, ValueVector, IndexVector>::iterator, typename MaxSearchTree<K, V, KeyVector, ValueVector, IndexVector>::iterator>
MaxSearchTree<K, V, KeyVector, ValueVector, IndexVector>::equal_range(const K& search_key) const {
    
    if (debug_mst) {
        std::cerr << "finding equal range for key ";
//        std::cerr << key;
        std::cerr << '\n';
    }
    
    size_t lower = -1, upper = -1;
    
    // find the lowest node that is inside the range
    size_t cursor = 0;
    while (cursor < size()) {
        if (key[cursor] == search_key) {
            if (debug_mst) {
                std::cerr << "found new lower bound " << cursor << '\n';
            }
            lower = cursor;
            cursor = left(cursor);
        }
        else if (key[cursor] > search_key) {
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
    while (cursor < size()) {
        if (key[cursor] > search_key) {
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

template<typename K, typename V, class KeyVector, class ValueVector, class IndexVector>
void MaxSearchTree<K, V, KeyVector, ValueVector, IndexVector>::reidentify_subtree_max(size_t x) {
    size_t new_max = x;
    size_t l = left(x);
    if (l < size() && value[subtree_max[l]] > value[new_max]) {
        new_max = subtree_max[l];
    }
    size_t r = right(x);
    if (r < size() && value[subtree_max[r]] > value[new_max]) {
        new_max = subtree_max[r];
    }
    subtree_max[x] = new_max;
}

template<typename K, typename V, class KeyVector, class ValueVector, class IndexVector>
void MaxSearchTree<K, V, KeyVector, ValueVector, IndexVector>::update(const MaxSearchTree<K, V, KeyVector, ValueVector, IndexVector>::iterator& it, const V& new_value) {
    
    if (debug_mst) {
        std::cerr << "performing update\n";
    }
    
    // TODO: these could be simplified into only the second case, but the first one seems
    // more likely to me, and it is a simpler loop...
    
    if (new_value > value[subtree_max[it.i]]) {
        // we're making a new max for the subtree rooted here
        subtree_max[it.i] = it.i;
        
        // this node might become the new max on the path to the root
        size_t here = it.i;
        while (here != 0) {
            here = parent(here);
            if (new_value > value[subtree_max[here]]) {
                subtree_max[here] = it.i;
            }
            else {
                break;
            }
        }
        value[it.i] = new_value;
    }
    else {
        // this value is less than the max
        value[it.i] = new_value;
        if (subtree_max[it.i] == it.i) {
            // the max was formerly at this node, but it might not be now
            reidentify_subtree_max(it.i);
            
            // we may also need to update the max toward the root
            size_t here = it.i;
            while (here != 0) {
                here = parent(here);
                if (subtree_max[here] != it.i) {
                    // the other subtree has a larger value anyway, we can stop early
                    break;
                }
                reidentify_subtree_max(here);
            }
        }
    }
}

template<typename K, typename V, class KeyVector, class ValueVector, class IndexVector>
typename MaxSearchTree<K, V, KeyVector, ValueVector, IndexVector>::iterator MaxSearchTree<K, V, KeyVector, ValueVector, IndexVector>::range_max(const K& lo, const K& hi) const {
    
    if (debug_mst) {
        std::cerr << "searching for range max\n";
    }
    
    // traverse downward until finding the first value in the interval
    size_t cursor = 0;
    while (cursor < size() && (key[cursor] < lo || key[cursor] >= hi)) {
        if (key[cursor] >= lo) {
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
    
    if (cursor >= size()) {
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
    while (left_cursor < size()) {
        if (key[left_cursor] >= lo) {
            if (value[left_cursor] > value[max_idx]) {
                max_idx = left_cursor;
            }
            size_t r = right(left_cursor);
            if (r < size() && value[subtree_max[r]] > value[max_idx]) {
                max_idx = subtree_max[r];
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
    while (right_cursor < size()) {
        if (key[right_cursor] < hi) {
            if (value[right_cursor] > value[max_idx]) {
                max_idx = right_cursor;
            }
            size_t l = left(right_cursor);
            if (l < size() && value[subtree_max[l]] > value[max_idx]) {
                max_idx = subtree_max[l];
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

template<typename K, typename V, class KeyVector, class ValueVector, class IndexVector>
size_t MaxSearchTree<K, V, KeyVector, ValueVector, IndexVector>::size() const {
    return key.size();
}

template<typename K, typename V, class KeyVector, class ValueVector, class IndexVector>
bool MaxSearchTree<K, V, KeyVector, ValueVector, IndexVector>::empty() const {
    return key.empty();
}

template<typename K, typename V, class KeyVector, class ValueVector, class IndexVector>
MaxSearchTree<K, V, KeyVector, ValueVector, IndexVector>::iterator::iterator(const MaxSearchTree<K, V, KeyVector, ValueVector, IndexVector>& iteratee, size_t i) :
    iteratee(&iteratee), i(i)
{
    // nothing to do
}

template<typename K, typename V, class KeyVector, class ValueVector, class IndexVector>
typename MaxSearchTree<K, V, KeyVector, ValueVector, IndexVector>::iterator& MaxSearchTree<K, V, KeyVector, ValueVector, IndexVector>::iterator::operator++() {
    size_t r = MaxSearchTree<K, V, KeyVector, ValueVector, IndexVector>::right(i);
    if (r < iteratee->size()) {
        // go to leftmost node of right subtree
        i = r;
        while (MaxSearchTree<K, V, KeyVector, ValueVector, IndexVector>::left(i) < iteratee->size()) {
            i = MaxSearchTree<K, V, KeyVector, ValueVector, IndexVector>::left(i);
        }
    }
    else if (i == 0) {
        // we are the root, and there is no right branch, this is the end
        i = iteratee->size();
    }
    else {
        // walk until we cross an edge rightward
        while (true) {
            size_t p = MaxSearchTree<K, V, KeyVector, ValueVector, IndexVector>::parent(i);
            if (i == MaxSearchTree<K, V, KeyVector, ValueVector, IndexVector>::left(p)) {
                // edge is rightward
                i = p;
                break;
            }
            if (p == 0) {
                // we walked only leftward edges to the root, we must have been at the end
                i = iteratee->size();
                break;
            }
            i = p;
        }
    }
    return *this;
}

template<typename K, typename V, class KeyVector, class ValueVector, class IndexVector>
const std::pair<K, V> MaxSearchTree<K, V, KeyVector, ValueVector, IndexVector>::iterator::operator*() const {
    return std::pair<K, V>(iteratee->key[i], iteratee->value[i]);
}

template<typename K, typename V, class KeyVector, class ValueVector, class IndexVector>
bool MaxSearchTree<K, V, KeyVector, ValueVector, IndexVector>::iterator::operator==(const iterator& other) const {
    return other.iteratee == this->iteratee && other.i == this->i;
}

template<typename K, typename V, class KeyVector, class ValueVector, class IndexVector>
bool MaxSearchTree<K, V, KeyVector, ValueVector, IndexVector>::iterator::operator!=(const iterator& other) const {
    return !(*this == other);
}


template<typename K, typename V, class KeyVector, class ValueVector, class IndexVector>
size_t MaxSearchTree<K, V, KeyVector, ValueVector, IndexVector>::memory_size() const {
    return get_vec_memory_size(key) + get_vec_memory_size(value) + get_vec_memory_size(subtree_max);
}
    
}

#endif /* centrolign_max_search_tree_hpp */
