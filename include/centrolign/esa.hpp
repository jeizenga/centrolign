#ifndef centrolign_path_esa_hpp
#define centrolign_path_esa_hpp

#include <vector>
#include <cstdint>
#include <algorithm>

#include "centrolign/modify_graph.hpp"

namespace centrolign {


/*
 * A conceptual node in the suffix tree
 */
struct SANode {
    SANode(size_t begin, size_t end) : begin(begin), end(end) { }
    SANode() = default;
    ~SANode() = default;
    inline bool is_leaf() const;
    inline bool operator<(const SANode& other) const;
    inline bool operator==(const SANode& other) const;
    inline bool operator!=(const SANode& other) const;
    
    size_t begin = -1;
    size_t end = -1;
};

/*
 * Shared components for GESA and PathESA
 */
class ESA {
public:
    
    
protected:
    
    static const bool debug_esa;
    
    void construct_child_array();
    void construct_suffix_links();
    
    inline SANode root() const;
    inline size_t first_l_index(const SANode& node) const;
    
    inline bool child_array_is_up(size_t i) const;
    inline bool child_array_is_down(size_t i) const;
    inline bool child_array_is_l_index(size_t i) const;
    
    inline std::pair<size_t, size_t> st_node_annotation_idx(const SANode& node) const;
    
    std::vector<size_t> child_array;
    std::vector<size_t> lcp_array;
    std::vector<SANode> suffix_links;
    
};

/*
 * Template implementations
 */


bool SANode::is_leaf() const {
    return begin == end;
}

bool SANode::operator<(const SANode& other) const {
    return begin < other.begin;
}

bool SANode::operator==(const SANode& other) const {
    return begin == other.begin && end == other.end;
}

bool SANode::operator!=(const SANode& other) const {
    return !(*this == other);
}

inline SANode ESA::root() const {
    return SANode(0, lcp_array.size() - 1);
}

inline size_t ESA::first_l_index(const SANode& node) const {
    if (node == root()) {
        // a special case
        // TODO: which we could maybe get around by storing a special .up value at the end?
        return child_array.front();
    }
    else if (child_array_is_down(node.begin) && child_array[node.begin] <= node.end) {
        // TODO: i'm not really sure why we also need the <= end check, although this
        // agrees with the paper...
        
        // we prefer the down value, because the up value could belong to an ancestor
        // interval
        return child_array[node.begin];
    }
    else {
        return child_array[node.end];
    }
}

inline std::pair<size_t, size_t> ESA::st_node_annotation_idx(const SANode& node) const {
    if (node.is_leaf()) {
        return std::pair<size_t, size_t>(1, node.begin);
    }
    else {
        return std::pair<size_t, size_t>(0, first_l_index(node));
    }
}

inline bool ESA::child_array_is_down(size_t i) const {
    return i < child_array.size() ? (child_array[i] > i && lcp_array[child_array[i]] != lcp_array[i]) : false;
}

inline bool ESA::child_array_is_up(size_t i) const {
    // other pointers are all downward
    return i < child_array.size() ? child_array[i] <= i : false;
}

inline bool ESA::child_array_is_l_index(size_t i) const {
    return i < child_array.size() ? (child_array[i] > i && lcp_array[child_array[i]] == lcp_array[i]) : false;
}

}

#endif /* centrolign_path_esa_hpp */
