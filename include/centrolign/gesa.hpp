#ifndef centrolign_gesa_hpp
#define centrolign_gesa_hpp

#include <vector>
#include <utility>
#include <algorithm>
#include <unordered_set>
#include <iostream>
#include <array>

#include "centrolign/graph.hpp"
#include "centrolign/modify_graph.hpp"
#include "centrolign/determinize.hpp"
#include "centrolign/path_graph.hpp"

namespace centrolign {

/*
 * A conceptual node in the suffix tree over prefixes
 */
struct GESANode {
    GESANode(size_t begin, size_t end) : begin(begin), end(end) { }
    GESANode() = default;
    ~GESANode() = default;
    inline bool is_leaf() const;
    inline bool operator<(const GESANode& other) const;
    inline bool operator==(const GESANode& other) const;
    inline bool operator!=(const GESANode& other) const;
    
    size_t begin = -1;
    size_t end = -1;
};

/*
 * An enhanced suffix array built over a prefix-sorted automaton
 */
class GESA {
public:
    
    // construct GESA for a single graph
    template<class BGraph>
    GESA(const BGraph& graph);
    // construct GESA for multiple graphs
    template<class BGraph>
    GESA(const std::vector<const BGraph*>& graphs);
    GESA() = default;
    ~GESA() = default;
    
    // the number of components used to construct
    size_t component_size() const;
    char sentinel(size_t component, bool beginning) const;
    
    // the length of the shared prefix of an internal node, or for a leaf the length
    // of its minimal unique prefix
    inline size_t depth(const GESANode& node) const;
    
protected:
    
    static bool debug_gesa;
    
    void print(std::ostream& out) const;
    
    template<class BGraph>
    GESA(const BGraph* const* const graphs, size_t num_graphs);
        
    void construct_child_array();
    inline bool child_array_is_up(size_t i) const;
    inline bool child_array_is_down(size_t i) const;
    inline bool child_array_is_l_index(size_t i) const;
    void construct_suffix_links();
    void compute_subtree_counts();
    
    // note: only valid at internal nodes
    inline size_t first_l_index(const GESANode& node) const;
    // valid at any node
    inline std::pair<size_t, size_t> st_node_annotation_idx(const GESANode& node) const;
     
    std::vector<GESANode> children(const GESANode& parent) const;
    
    GESANode root() const;
    
    GESANode link(const GESANode& node) const;
    
    std::vector<std::pair<char, char>> component_sentinels;
    std::vector<size_t> component_sizes;
    
    std::vector<uint64_t> ranked_node_ids;
    std::vector<size_t> lcp_array;
    std::vector<size_t> child_array;
    std::array<std::vector<GESANode>, 2> suffix_links;
    std::vector<std::array<std::vector<uint64_t>, 2>> component_subtree_counts;
    std::vector<std::vector<uint64_t>> edges;
    
    
};

/*
 * Template and inline implementations
 */

bool GESANode::is_leaf() const {
    return begin == end;
}

bool GESANode::operator<(const GESANode& other) const {
    return begin < other.begin;
}

bool GESANode::operator==(const GESANode& other) const {
    return begin == other.begin && end == other.end;
}

bool GESANode::operator!=(const GESANode& other) const {
    return !(*this == other);
}

template<class BGraph>
GESA::GESA(const std::vector<const BGraph*>& graphs) : GESA(graphs.data(), graphs.size()) {
    // nothing to do besides dispatch
}

template<class BGraph>
GESA::GESA(const BGraph& graph) : GESA(&(&graph), 1) {
    // nothing to do besides dispatch
}

template<class BGraph>
GESA::GESA(const BGraph* const* const graphs, size_t num_graphs) :
    component_sentinels(num_graphs),
    component_sizes(num_graphs)
{
    
    // find the maximum character among all the graphs
    // TODO: in some cases we'll already know this, in which case this feels
    // a bit wasteful
    char max_char = 0;
    for (size_t i = 0; i < num_graphs; ++i) {
        const BGraph& graph = *graphs[i];
        for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
            max_char = std::max(max_char, graph.base(node_id));
        }
    }
            
    // make a single graph with all of the components, with sentinels attached
    BaseGraph joined;
    std::vector<uint64_t> path_id_range(num_graphs + 1, 0);
    for (size_t i = 0; i < num_graphs; ++i) {
        // add a new component
        const BGraph& graph = *graphs[i];
        uint64_t first_path_id = joined.path_size();
        uint64_t first_node_id = joined.node_size();
        append_component(joined, graph);
        
        // record the end of the path range
        path_id_range[i + 1] = joined.path_size();
        component_sizes[i] = joined.node_size() - first_node_id;
    }
    
    if (debug_gesa) {
        std::cerr << "joined graph for GESA indexing:\n";
        print_graph(joined, std::cerr);
    }
    
    // initialize a path graph
    PathGraph path_graph(joined);
    
    // do prefix doubling until prefix sorted
    while (!path_graph.is_prefix_sorted()) {
        path_graph = PathGraph(path_graph);
    }
    // TODO: if i want to reduce to a prefix range sorted graph it would happen here
        
    // reassign node IDs to be ordered by prefix rank
    path_graph.order_by_rank();
    
    // use the original graph to add the edges in
    path_graph.construct_edges(joined);
    
    // get the node array (analog to suffix array)
    ranked_node_ids.resize(path_graph.node_size());
    for (size_t i = 0; i < path_graph.node_size(); ++i) {
        ranked_node_ids[i] = path_graph.from(i);
    }
    // TODO: do i need to have the M and F bit vectors for edges?
    // TODO: i think i don't need the BWT, but it could be used here
    // - for now I'm just maintaining the edges explicitly
    
    // steal the LCP array that is built alongside the path graph
    lcp_array = std::move(path_graph.lcp_array);
    // add an initial 0 to match the definitions in Abouelhoda, et al. (2004)
    // TODO: possible to not do this, but idk if it's worth reformulating the algorithms
    lcp_array.insert(lcp_array.begin(), 0);
    
    // steal the edges as well
    edges.resize(path_graph.node_size());
    for (size_t i = 0; i < edges.size(); ++i) {
        edges[i] = std::move(path_graph.edges[i].next);
    }
    
    // build structure for traversing downward and assigning annotations
    construct_child_array();
    
    // get the suffix links
    construct_suffix_links();
    
    // get subtree counts for each
    compute_subtree_counts();
    
    if (debug_gesa) {
        print(std::cerr);
    }
}

inline bool GESA::child_array_is_down(size_t i) const {
    return i < child_array.size() ? (child_array[i] > i && lcp_array[child_array[i]] != lcp_array[i]) : false;
}

inline bool GESA::child_array_is_up(size_t i) const {
    // other pointers are all downward
    return i < child_array.size() ? child_array[i] <= i : false;
}

inline bool GESA::child_array_is_l_index(size_t i) const {
    return i < child_array.size() ? (child_array[i] > i && lcp_array[child_array[i]] == lcp_array[i]) : false;
}

inline size_t GESA::depth(const GESANode& node) const {
    if (node.is_leaf()) {
        // longer of the two adjacent matches
        size_t l = lcp_array[node.begin];
        if (node.begin + 1 < lcp_array.size()) {
            l = std::max(l, lcp_array[node.begin + 1]);
        }
        // add 1 to make the prefix unique
        return l + 1;
    }
    else {
        return lcp_array[first_l_index(node)];
    }
}

inline size_t GESA::first_l_index(const GESANode& node) const {
    if (node == root()) {
        // a special case
        // TODO: which we could maybe get around by storing a special .up value at the end?
        return child_array.front();
    }
    else if (child_array_is_down(node.begin)) {
        // we prefer the down value, because the up value could belong to an ancestor
        // interval
        return child_array[node.begin];
    }
    else {
        return child_array[node.end];
    }
}

inline std::pair<size_t, size_t> GESA::st_node_annotation_idx(const GESANode& node) const {
    if (node.is_leaf()) {
        return std::pair<size_t, size_t>(0, node.begin);
    }
    else {
        return std::pair<size_t, size_t>(1, first_l_index(node));
    }
}

}

#endif /* centrolign_gesa_hpp */
