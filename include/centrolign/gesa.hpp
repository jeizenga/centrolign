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
#include "centrolign/path_graph.hpp"
#include "centrolign/logging.hpp"

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
    inline size_t component_size() const;
    
    // the length of the shared prefix of an internal node, or for a leaf the length
    // of its minimal unique prefix
    inline size_t depth(const GESANode& node) const;
    
    // return the locations of minimal sequences that occur on multiple components,
    // but at most max_count many times on any component, paired with the length of
    // the minimal match
    std::vector<std::pair<GESANode, size_t>> minimal_rare_matches(size_t max_count) const;
    
    // returns a vector of the number of counts among the unique prefixes in each component
    const std::vector<uint64_t>& component_counts(const GESANode& node) const;
    
    // walk the label of a node out in each of the graphs and return the resulting match
    // each match consists of the component index and a list of node IDs in the original graph
    std::vector<std::pair<size_t, std::vector<uint64_t>>> walk_matches(const GESANode& node,
                                                                       size_t length) const;
    
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
    void label_edges(size_t doubling_steps, const BaseGraph& joined);
    
    // note: only valid at internal nodes
    inline size_t first_l_index(const GESANode& node) const;
    // valid at any node
    inline std::pair<size_t, size_t> st_node_annotation_idx(const GESANode& node) const;
     
    std::vector<GESANode> children(const GESANode& parent) const;
    
    inline char label(const GESANode& node) const;
    inline GESANode child(const GESANode& parent, char label) const;
    
    inline GESANode root() const;
    
    inline GESANode link(const GESANode& node) const;
    
    std::vector<size_t> component_ranges;
    
    std::vector<uint64_t> ranked_node_ids;
    std::vector<size_t> lcp_array;
    std::vector<size_t> child_array;
    std::vector<GESANode> suffix_links;
    // TODO: the leaf annotations are pretty pointless here...
    std::array<std::vector<std::vector<uint64_t>>, 2> component_subtree_counts;
    // TODO: this could be replaced by the M and F bit vectors with rank/select support for Psi
    // TODO: do i really need any more than 1 edge for my purposes?
    std::vector<std::vector<uint64_t>> edges;
    // TODO: this could also be replaced by binary search on the component range vector,
    // which might even be faster for small numbers of components due to memory access
    std::vector<uint16_t> node_to_comp;
    // the label of the downward edge to the node
    std::array<std::vector<unsigned char>, 2> edge_label;
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
    component_ranges(num_graphs + 1, 0)
{
            
    // make a single graph with all of the components
    BaseGraph joined;
    for (size_t i = 0; i < num_graphs; ++i) {
        // add a new component
        const BGraph& graph = *graphs[i];
        append_component(joined, graph);
        
        // record the end of the node range
        component_ranges[i + 1] = joined.node_size();
    }
    
    node_to_comp.resize(component_ranges.back(), 0);
    uint64_t node_id = 0;
    for (size_t i = 0, j = 0; i + 1 < component_ranges.size(); ++i) {
        for (size_t n = component_ranges[i + 1]; j < n; ++j) {
            node_to_comp[node_id] = i;
            ++node_id;
        }
    }
    
    if (debug_gesa) {
        std::cerr << "joined graph for GESA indexing:\n";
        print_graph(joined, std::cerr);
    }
    
    logging::log(logging::Debug, "Initializing path graph");
    
    // initialize a path graph
    PathGraph path_graph(joined);
    
    // do prefix doubling until prefix sorted
    while (!path_graph.is_prefix_sorted()) {
        
        logging::log(logging::Debug, "Performing doubling step " + std::to_string(path_graph.doubling_step + 1));
        
        path_graph = PathGraph(path_graph);
    }
    // TODO: if i want to reduce to a prefix range sorted graph it would happen here
        
    logging::log(logging::Debug, "Ranks stabilized, constructing edges");
    
    // reassign node IDs to be ordered by prefix rank
    path_graph.order_by_rank();
    
    // use the original graph to add the edges in
    path_graph.construct_edges(joined);
    
    // get the node array (analog to suffix array)
    ranked_node_ids.resize(path_graph.node_size());
    for (size_t i = 0; i < path_graph.node_size(); ++i) {
        ranked_node_ids[i] = path_graph.from(i);
    }
    
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
    
    logging::log(logging::Debug, "Constructing child array");
    
    // build structure for traversing downward and assigning annotations
    construct_child_array();
    
    logging::log(logging::Debug, "Constructing suffix links");
    
    // get the suffix links
    construct_suffix_links();
    
    logging::log(logging::Debug, "Labeling edges");
    
    label_edges(path_graph.doubling_step, joined);
    
    logging::log(logging::Debug, "Computing subtree counts");
    
    // get subtree counts for each component
    compute_subtree_counts();
    
    logging::log(logging::Debug, "Finished constructing GESA");
    
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

inline std::pair<size_t, size_t> GESA::st_node_annotation_idx(const GESANode& node) const {
    if (node.is_leaf()) {
        return std::pair<size_t, size_t>(1, node.begin);
    }
    else {
        return std::pair<size_t, size_t>(0, first_l_index(node));
    }
}

size_t GESA::component_size() const {
    return component_ranges.size() - 1;
}

GESANode GESA::root() const {
    return GESANode(0, lcp_array.size() - 1);
}

GESANode GESA::link(const GESANode& node) const {
    size_t i, j;
    std::tie(i, j) = st_node_annotation_idx(node);
    return suffix_links[j];
}

inline char GESA::label(const GESANode& node) const {
    size_t i, j;
    std::tie(i, j) = st_node_annotation_idx(node);
    return edge_label[i][j];
}

inline GESANode GESA::child(const GESANode& parent, char _label) const {
    for (auto& c : children(parent)) {
        if (label(c) == _label) {
            return c;
        }
    }
    // there is no corresponding child
    return GESANode(-1, -1);
}
}

#endif /* centrolign_gesa_hpp */
