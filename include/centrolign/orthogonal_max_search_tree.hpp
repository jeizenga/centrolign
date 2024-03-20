#ifndef centrolign_orthogonal_max_search_tree_hpp
#define centrolign_orthogonal_max_search_tree_hpp

#include "centrolign/max_search_tree.hpp"

#include <vector>
#include <cstdint>
#include <iostream>

namespace centrolign {




/*
 * Static-topology search tree that supports 2D orthogonal range max queries
 * and value updates
 */
template<typename K1, typename K2, typename V, typename UIntSize = size_t>
class OrthogonalMaxSearchTree {
public:
    
    // vector will be sorted, key pairs must be unique
    OrthogonalMaxSearchTree(std::vector<std::tuple<K1, K2, V>>& values);
    OrthogonalMaxSearchTree() = default;
    ~OrthogonalMaxSearchTree() = default;
    
    class iterator; // forward declaration
    
    // standard iteration interface
    iterator begin() const;
    iterator end() const;
    
    // return iterator that has this key
    iterator find(const K1& key1, const K2& key2) const;
    // change value at the key the iterator points to
    void update(const iterator& it, const V& value);
    // returns iterator to max value in key range [lo1, hi1) x [lo2, hi2)
    iterator range_max(const K1& lo1, const K1& hi1,
                       const K2& lo2, const K2& hi2) const;
    
    class iterator {
    public:
        iterator() = default;
        ~iterator() = default;
        
        iterator& operator++();
        const std::tuple<K1, K2, V>& operator*() const;
        const std::tuple<K1, K2, V>* operator->() const;
        bool operator==(const iterator& other) const;
        bool operator!=(const iterator& other) const;
        
    private:
        friend class OrthogonalMaxSearchTree<K1, K2, V, UIntSize>;
        
        // internal constructor
        iterator(const OrthogonalMaxSearchTree<K1, K2, V, UIntSize>& iteratee, UIntSize i);
        
        const OrthogonalMaxSearchTree<K1, K2, V, UIntSize>* iteratee = nullptr;
        UIntSize i = 0;
    };
    
    size_t memory_size() const;
    
private:
    
    static const bool debug_mst = false;
    
    struct OuterNode {
        OuterNode() = default;
        ~OuterNode() = default;
        std::tuple<K1, K2, V> key_value;
        MaxSearchTree<K2, std::pair<V, UIntSize>, UIntSize> cross_tree;
    };
    
    static inline UIntSize left(UIntSize x);
    static inline UIntSize right(UIntSize x);
    static inline UIntSize parent(UIntSize x);
    
    std::vector<OuterNode> nodes;
    
    friend class iterator;
    
};






/*
 * Template implementations
 */


// this is useful for debugging
//inline std::ostream& operator<<(std::ostream& out, const std::tuple<int64_t, size_t, size_t, size_t>& k) {
//    return out << "(" << std::get<0>(k) << ","  << std::get<1>(k) << ","  << std::get<2>(k) << ","  << std::get<3>(k) << ")";
//}

template<typename K1, typename K2, typename V, typename UIntSize>
OrthogonalMaxSearchTree<K1, K2, V, UIntSize>::OrthogonalMaxSearchTree(std::vector<std::tuple<K1, K2, V>>& values) {
    
    // TODO: repetitive with MaxSearchTree
    
    if (size_t(std::numeric_limits<UIntSize>::max()) < values.size()) {
        throw std::runtime_error("OrthogonalMaxSearchTree integer size is too small to represent its indexes.");
    }
    
    // handle this as a special case so we can
    if (values.empty()) {
        return;
    }
    
    // comparing only on key, not value
    auto cmp = [](const std::tuple<K1, K2, V>& a, const std::tuple<K1, K2, V>& b) {
        return (std::get<0>(a) < std::get<0>(b) ||
                (std::get<0>(a) == std::get<0>(b) && std::get<1>(a) < std::get<1>(b)));
    };
    if (!std::is_sorted(values.begin(), values.end(), cmp)) {
        std::stable_sort(values.begin(), values.end(), cmp);
    }
    
    // figure out the height of the tree
    size_t height = 0;
    while ((1 << height) - 1 < values.size()) { // max capacity of a tree with this height
        ++height;
    }
    
    if (debug_mst) {
        std::cerr << "building tree of height " << height << " for " << values.size() << " key-value pairs\n";
    }
    
    // init the implicit tree
    nodes.resize(values.size());
    
    // the next item in the vector we will assign to a node
    UIntSize vec_idx = 0;
    // records of (node_idx, queued left)
    std::vector<std::pair<UIntSize, bool>> stack;
    stack.emplace_back(0, false);
    
    // the index of the node for the key-value pair in this ordinal position
    std::vector<UIntSize> indexes(values.size());
        
    // in-order traversal of the tree through a stack to assign keys
    while (!stack.empty()) {
        auto& top = stack.back();
        if (!top.second) {
            // we haven't traversed the left branch yet
            top.second = true;
            UIntSize l = left(top.first);
            if (l < nodes.size()) {
                stack.emplace_back(l, false);
            }
        }
        else {
            // the left is traversed, this is the node's position
            UIntSize node_idx = top.first;
            indexes[vec_idx] = node_idx;
            nodes[node_idx].key_value = values[vec_idx++];
            stack.pop_back();
            // queue up the right branch as well
            UIntSize r = right(node_idx);
            if (r < nodes.size()) {
                stack.emplace_back(r, false);
            }
        }
    }
    
    // FIXME: this will break for duplicate K1 values
    
    // now a depth-first traversal divvying up the values
    std::vector<std::tuple<UIntSize, std::vector<std::tuple<K1, K2, V>>, std::vector<UIntSize>>> outer_stack;
    outer_stack.emplace_back(0, values, std::move(indexes));
    while (!outer_stack.empty()) {
        
        UIntSize n = std::get<0>(outer_stack.back());
        auto subtree_values = std::move(std::get<1>(outer_stack.back()));
        auto subtree_indexes = std::move(std::get<2>(outer_stack.back()));
        outer_stack.pop_back();
        
        const auto& pivot = nodes[n].key_value;
        
        // the parts of the data that will go to the left and right children
        std::vector<std::tuple<K1, K2, V>> left_vals, right_vals;
        std::vector<UIntSize> left_indexes, right_indexes;
        
        // shim into MaxSearchTree
        std::vector<std::pair<K2, std::pair<V, UIntSize>>> max_tree_vals;
        
        max_tree_vals.reserve(subtree_values.size());
        left_vals.reserve(subtree_values.size() / 2);
        right_vals.reserve(subtree_values.size() / 2);
        left_indexes.reserve(subtree_values.size() / 2);
        right_indexes.reserve(subtree_values.size() / 2);
        
        for (size_t i = 0; i < subtree_values.size(); ++i) {
            const auto& value = subtree_values[i];
            max_tree_vals.emplace_back(std::get<1>(value), std::make_pair(std::get<2>(value), subtree_indexes[i]));
            if (value < pivot) {
                left_vals.emplace_back(value);
                left_indexes.emplace_back(subtree_indexes[i]);
            }
            else if (value > pivot) {
                right_vals.emplace_back(value);
                right_indexes.emplace_back(subtree_indexes[i]);
            }
        }
        // construct the axis 2 search tree
        nodes[n].cross_tree = std::move(MaxSearchTree<K2, std::pair<V, UIntSize>, UIntSize>(max_tree_vals));
        
        // queue up the children
        UIntSize l = left(n);
        if (l < nodes.size()) {
            outer_stack.emplace_back(l, std::move(left_vals), std::move(left_indexes));
            UIntSize r = right(n);
            if (r < nodes.size()) {
                outer_stack.emplace_back(r, std::move(right_vals), std::move(right_indexes));
            }
        }
    }
    
    if (debug_mst) {
        for (size_t i = 0; i < nodes.size(); ++i) {
            std::cerr << "node " << i << ": ";
            std::cerr << "? ";
            //std::cerr << std::get<0>(nodes[i].key_value);
            std::cerr << ", " << std::get<1>(nodes[i].key_value) << ": " << std::get<2>(nodes[i].key_value);
            std::cerr << '\n';
        }
    }
}

template<typename K1, typename K2, typename V, typename UIntSize>
inline UIntSize OrthogonalMaxSearchTree<K1, K2, V, UIntSize>::left(UIntSize x) {
    return 2 * x + 1;
}

template<typename K1, typename K2, typename V, typename UIntSize>
inline UIntSize OrthogonalMaxSearchTree<K1, K2, V, UIntSize>::right(UIntSize x) {
    return 2 * x + 2;
}

template<typename K1, typename K2, typename V, typename UIntSize>
inline UIntSize OrthogonalMaxSearchTree<K1, K2, V, UIntSize>::parent(UIntSize x) {
    return (x - 1) / 2;
}

template<typename K1, typename K2, typename V, typename UIntSize>
typename OrthogonalMaxSearchTree<K1, K2, V, UIntSize>::iterator OrthogonalMaxSearchTree<K1, K2, V, UIntSize>::begin() const {
    UIntSize i = 0;
    UIntSize l = left(i);
    while (l < nodes.size()) {
        i = l;
        l = left(i);
    }
    return iterator(*this, i);
}

template<typename K1, typename K2, typename V, typename UIntSize>
typename OrthogonalMaxSearchTree<K1, K2, V, UIntSize>::iterator OrthogonalMaxSearchTree<K1, K2, V, UIntSize>::end() const {
    return iterator(*this, nodes.size());
}

template<typename K1, typename K2, typename V, typename UIntSize>
typename OrthogonalMaxSearchTree<K1, K2, V, UIntSize>::iterator OrthogonalMaxSearchTree<K1, K2, V, UIntSize>::find(const K1& key1, const K2& key2) const {
    UIntSize cursor = 0;
    while (cursor < nodes.size()) {
        if (std::get<0>(nodes[cursor].key_value) == key1 && std::get<1>(nodes[cursor].key_value) == key2) {
            return iterator(*this, cursor);
        }
        else if (std::make_pair(std::get<0>(nodes[cursor].key_value),
                                std::get<1>(nodes[cursor].key_value)) > std::make_pair(key1, key2)) {
            cursor = left(cursor);
        }
        else {
            cursor = right(cursor);
        }
    }
    return end();
}


template<typename K1, typename K2, typename V, typename UIntSize>
void OrthogonalMaxSearchTree<K1, K2, V, UIntSize>::update(const OrthogonalMaxSearchTree<K1, K2, V, UIntSize>::iterator& it, const V& value) {
    
    if (debug_mst) {
        std::cerr << "updating key ";
        //std::cerr << std::get<0>(*it) << ", " << std::get<1>(*it);
        std::cerr << ", val " << std::get<2>(*it) << " to new value " << value << '\n';
    }
    
    // update all the cross trees that include this outer node
    
    std::get<2>(nodes[it.i].key_value) = value;
    for (UIntSize cursor = it.i; cursor < nodes.size(); cursor = parent(cursor)) {
        // note: counting on cursor underflowing at 0 to end iteration
        if (debug_mst) {
            std::cerr << "updating search tree at outer node " << cursor << " with key ";
            //std::cerr << std::get<0>(nodes[cursor].key_value) << ", " << std::get<1>(nodes[cursor].key_value);
            std::cerr << '\n';
        }
        
        // there may be multiple key pairs that have the same key2 value, so find the right one
        // FIXME: this can break the O(log n) run time for this subroutine
        auto cross_range = nodes[cursor].cross_tree.equal_range(std::get<1>(*it));
        while (cross_range.first->second.second != it.i) {
            ++cross_range.first;
        }
        nodes[cursor].cross_tree.update(cross_range.first, std::make_pair(value, it.i));
    }
}

template<typename K1, typename K2, typename V, typename UIntSize>
typename OrthogonalMaxSearchTree<K1, K2, V, UIntSize>::iterator
OrthogonalMaxSearchTree<K1, K2, V, UIntSize>::range_max(const K1& lo1, const K1& hi1,
                                                        const K2& lo2, const K2& hi2) const {
    
    if (debug_mst) {
        std::cerr << "start initial search for value in range ";
        //std::cerr << lo1 << ", " << hi1;
        std::cerr << "\n";
        if (!nodes.empty()) {
            std::cerr << "initial value:\n";
            //std::cerr << std::get<0>(nodes.front().key_value) << '\n';
        }
    }
    
    // traverse downward until finding the first value in the interval
    UIntSize cursor = 0;
    while (cursor < nodes.size() &&
           (std::get<0>(nodes[cursor].key_value) < lo1 ||
            std::get<0>(nodes[cursor].key_value) >= hi1)) {
        
        if (debug_mst) {
            std::cerr << "cursor at " << cursor << " with key1 ";
            //std::cerr << std::get<0>(nodes[cursor].key_value);
            std::cerr << '\n';
        }
        if (std::get<0>(nodes[cursor].key_value) >= hi1) {
            cursor = left(cursor);
        }
        else {
            cursor = right(cursor);
        }
    }
    
    if (cursor >= nodes.size()) {
        // we never found a value that was in the interval
        return end();
    }
    
    if (debug_mst) {
        std::cerr << "found contained key1";
        //std::cerr << " " << std::get<0>(nodes[cursor].key_value);
        std::cerr << " at index " << cursor << " with key2 " << std::get<1>(nodes[cursor].key_value) << '\n';
    }
    
    // the max will be either at an index (in outer tree) or iterator (in an inner tree)
    bool max_at_idx = false;
    bool max_at_iter = false;
    UIntSize max_idx = -1;
    typename MaxSearchTree<K2, std::pair<V, UIntSize>, UIntSize>::iterator max_iter;
    if (std::get<1>(nodes[cursor].key_value) >= lo2 && std::get<1>(nodes[cursor].key_value) < hi2) {
        max_at_idx = true;
        max_idx = cursor;
        if (debug_mst) {
            std::cerr << "initial index key2 " << std::get<1>(nodes[cursor].key_value) << " is in range, opt is " << std::get<2>(nodes[cursor].key_value) << '\n';
        }
    }
    
    // check if a new value beats the current best
    auto is_opt = [&](const V& val) -> bool {
        if (max_at_idx) {
            return val > std::get<2>(nodes[max_idx].key_value);
        }
        else if (max_at_iter) {
            return val > max_iter->second.first;
        }
        else {
            return true;
        }
    };
    
    UIntSize right_cursor = right(cursor);
    UIntSize left_cursor = left(cursor);
    if (debug_mst) {
        std::cerr << "start left search at " << left_cursor << '\n';
    }
    
    // split off the leftward traversal where rightward off-path edges are included entirely
    while (left_cursor < nodes.size()) {
        if (std::get<0>(nodes[left_cursor].key_value) >= lo1) {
            if (debug_mst) {
                std::cerr << "outer node " << left_cursor << " is in range with value " << std::get<2>(nodes[left_cursor].key_value) << '\n';
            }
            if (std::get<1>(nodes[left_cursor].key_value) >= lo2 &&
                std::get<1>(nodes[left_cursor].key_value) < hi2 &&
                is_opt(std::get<2>(nodes[left_cursor].key_value))) {
                
                max_idx = left_cursor;
                max_at_idx = true;
                max_at_iter = false;
                if (debug_mst) {
                    std::cerr << "find new opt at outer node " << max_idx << " with value " << std::get<2>(nodes[max_idx].key_value) << '\n';
                }
            }
            UIntSize r = right(left_cursor);
            if (r < nodes.size()) {
                // we're guaranteed that the key1 range is satisfied, so the second element
                // in this range is non-binding
                auto iter = nodes[r].cross_tree.range_max(lo2, hi2);
                if (debug_mst) {
                    std::cerr << "rightward tree has range max";
                    if (iter == nodes[r].cross_tree.end()) {
                        std::cerr << " none\n";
                    }
                    else {
                        //std::cerr << iter->first.first << ", " << iter->first.second;
                        std::cerr <<": " << iter->second.first << " (" << iter->second.second << ")" << '\n';
                    }
                    //std::cerr << "in ranges [" << lo2 << ", " << lo1 << ") x [" << hi2 << ", " << hi1 << ")\n";
                }
                if (iter != nodes[r].cross_tree.end() && is_opt(iter->second.first)) {
                    
                    max_iter = iter;
                    max_at_idx = false;
                    max_at_iter = true;
                    if (debug_mst) {
                        std::cerr << "find new opt at inner node with value ";
                        //std::cerr << max_iter->second;
                        std::cerr << '\n';
                    }
                }
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
        std::cerr << "start right search at " << right_cursor << '\n';
    }
    // and next the rightward traversal where leftward off-path edges are included entirely
    while (right_cursor < nodes.size()) {
        if (std::get<0>(nodes[right_cursor].key_value) < hi1) {
            if (debug_mst) {
                std::cerr << "outer node " << right_cursor << " is in range with value " << std::get<2>(nodes[right_cursor].key_value) << '\n';
            }
            if (std::get<1>(nodes[right_cursor].key_value) >= lo2 &&
                std::get<1>(nodes[right_cursor].key_value) < hi2 &&
                is_opt(std::get<2>(nodes[right_cursor].key_value))) {
                
                max_idx = right_cursor;
                max_at_idx = true;
                max_at_iter = false;
                if (debug_mst) {
                    std::cerr << "find new opt at outer node " << max_idx << " with value " << std::get<2>(nodes[max_idx].key_value) << '\n';
                }
            }
            UIntSize l = left(right_cursor);
            if (l < nodes.size()) {
                // we're guaranteed that the key1 range is satisfied, so the second element
                // in this range is non-binding
                auto iter = nodes[l].cross_tree.range_max(lo2, hi2);
                if (debug_mst) {
                    std::cerr << "leftward tree has range max";
                    if (iter == nodes[l].cross_tree.end()) {
                        std::cerr << " none\n";
                    }
                    else {
                        //std::cerr << iter->first.first << ", " << iter->first.second;
                        std::cerr <<": " << iter->second.first << " (" << iter->second.second << ")" << '\n';
                    }
                    //std::cerr << "in ranges [" << lo2 << ", " << lo1 << ") x [" << hi2 << ", " << hi1 << ")\n";
                }
                if (iter != nodes[l].cross_tree.end() && is_opt(iter->second.first)) {
                    
                    max_iter = iter;
                    max_at_idx = false;
                    max_at_iter = true;
                    if (debug_mst) {
                        std::cerr << "find new opt at inner node with value ";
                        //std::cerr << max_iter->second;
                        std::cerr << '\n';
                    }
                }
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
        std::cerr << "max is from index? " << max_at_idx << " " << max_idx << ", or iter? " << max_at_iter << '\n';
    }
    
    if (max_at_idx) {
        // we're already at an outer tree node
        return iterator(*this, max_idx);
    }
    else if (max_at_iter) {
        // find the outer node that corresponds to the inner node we detected the max at
        return iterator(*this, max_iter->second.second);
    }
    else {
        return end();
    }
}
template<typename K1, typename K2, typename V, typename UIntSize>
OrthogonalMaxSearchTree<K1, K2, V, UIntSize>::iterator::iterator(const OrthogonalMaxSearchTree<K1, K2, V, UIntSize>& iteratee, UIntSize i) :
    iteratee(&iteratee), i(i)
{
    // nothing to do
}


template<typename K1, typename K2, typename V, typename UIntSize>
typename OrthogonalMaxSearchTree<K1, K2, V, UIntSize>::iterator& OrthogonalMaxSearchTree<K1, K2, V, UIntSize>::iterator::operator++() {
    UIntSize r = OrthogonalMaxSearchTree<K1, K2, V, UIntSize>::right(i);
    if (r < iteratee->nodes.size()) {
        // go to leftmost node of right subtree
        i = r;
        while (OrthogonalMaxSearchTree<K1, K2, V, UIntSize>::left(i) < iteratee->nodes.size()) {
            i = OrthogonalMaxSearchTree<K1, K2, V, UIntSize>::left(i);
        }
    }
    else if (i == 0) {
        // we are the root, and there is no right branch, this is the end
        i = iteratee->nodes.size();
    }
    else {
        // walk until we cross an edge rightward
        while (true) {
            UIntSize p = OrthogonalMaxSearchTree<K1, K2, V, UIntSize>::parent(i);
            if (i == OrthogonalMaxSearchTree<K1, K2, V, UIntSize>::left(p)) {
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

template<typename K1, typename K2, typename V, typename UIntSize>
const std::tuple<K1, K2, V>& OrthogonalMaxSearchTree<K1, K2, V, UIntSize>::iterator::operator*() const {
    return iteratee->nodes[i].key_value;
}

template<typename K1, typename K2, typename V, typename UIntSize>
const std::tuple<K1, K2, V>* OrthogonalMaxSearchTree<K1, K2, V, UIntSize>::iterator::operator->() const {
    return &(iteratee->nodes[i].key_value);
}

template<typename K1, typename K2, typename V, typename UIntSize>
bool OrthogonalMaxSearchTree<K1, K2, V, UIntSize>::iterator::operator==(const iterator& other) const {
    return other.iteratee == this->iteratee && other.i == this->i;
}

template<typename K1, typename K2, typename V, typename UIntSize>
bool OrthogonalMaxSearchTree<K1, K2, V, UIntSize>::iterator::operator!=(const iterator& other) const {
    return !(*this == other);
}

template<typename K1, typename K2, typename V, typename UIntSize>
size_t OrthogonalMaxSearchTree<K1, K2, V, UIntSize>::memory_size() const {
    size_t mem_size = (nodes.capacity() - nodes.size()) * sizeof(OrthogonalMaxSearchTree<K1, K2, V, UIntSize>::OuterNode);
    for (const auto& node : nodes) {
        mem_size += sizeof(node.key_value) + node.cross_tree.memory_size();
    }
    return mem_size;
}

}

#endif /* centrolign_orthogonal_max_search_tree_hpp */
