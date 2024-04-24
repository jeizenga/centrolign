#ifndef centrolign_gesa_hpp
#define centrolign_gesa_hpp

#include <vector>
#include <utility>
#include <algorithm>
#include <unordered_set>
#include <iostream>
#include <array>
#include <tuple>
#include <limits>

#include "centrolign/graph.hpp"
#include "centrolign/modify_graph.hpp"
#include "centrolign/path_graph.hpp"
#include "centrolign/logging.hpp"
#include "centrolign/esa.hpp"

namespace centrolign {

/*
 * An enhanced suffix array built over a prefix-sorted automaton
 */
class GESA : public ESA {
public:
    
    // construct GESA for a single graph
    template<class BGraph>
    GESA(const BGraph& graph, size_t size_limit = std::numeric_limits<size_t>::max());
    // construct GESA for a single graph with back-translated node IDs
    template<class BGraph>
    GESA(const BGraph& graph, const std::vector<uint64_t>& back_translation,
         size_t size_limit = std::numeric_limits<size_t>::max());
    // construct GESA for multiple graphs
    template<class BGraph>
    GESA(const std::vector<const BGraph*>& graphs, size_t size_limit = std::numeric_limits<size_t>::max());
    // construct GESA for multiple graphs with back-translated node IDs
    template<class BGraph>
    GESA(const std::vector<const BGraph*>& graphs,
         const std::vector<const std::vector<uint64_t>*>& back_translations,
         size_t size_limit = std::numeric_limits<size_t>::max());
    GESA() = default;
    ~GESA() = default;
    
    // return the locations of minimal sequences that occur on multiple components,
    // but at most max_count many times on any component, along with the length of
    // the minimal match, and the counts on each component
    // optionally choose either CSS or RUQ index internally TODO: ugly interface
    std::vector<std::tuple<SANode, size_t, std::vector<uint64_t>>> minimal_rare_matches(size_t max_count,
                                                                                        bool use_css = true) const;
    
    // walk the label of a node out in each of the graphs and return the resulting match
    // each match consists of the component index and a list of node IDs in the original graph
    std::vector<std::pair<size_t, std::vector<uint64_t>>> walk_matches(const SANode& node,
                                                                       size_t length) const;
    
    size_t memory_size() const;
    
protected:
    
    static bool debug_gesa;
    
    void print(std::ostream& out) const;
    
    template<class BGraph>
    GESA(const BGraph* const* const graphs, size_t num_graphs,
         const std::vector<uint64_t>* const* const back_translations, size_t size_limit);
        
    inline void construct_suffix_links();
    
    void label_edges(size_t doubling_steps, const BaseGraph& joined, const PathGraph& path_graph);
         
    inline char label(const SANode& node) const;
    inline SANode child(const SANode& parent, char label) const;
    
    // TODO: this could be replaced by the M and F bit vectors with rank/select support for Psi
    // TODO: do i really need any more than 1 edge for my purposes?
    std::vector<std::vector<uint64_t>> edges;
    // the label of the downward edge to the node
    std::array<std::vector<unsigned char>, 2> edge_label;
};

/*
 * Custom exception for when path graph exceeds a size limit
 */
class GESASizeException : public std::exception {
public:
    GESASizeException() noexcept = default;
    ~GESASizeException() noexcept = default;
    const char* what() const noexcept;
    // the number of times that a node of the input grapha appeared as the source of
    // a PathGraph node in the incomplete PathGraph
    const std::vector<std::vector<uint64_t>>& from_counts() const;
    // the number of times that a node of the input grapha appeared as the source of
    // a PathGraph node in the previous, complete PathGraph
    const std::vector<std::vector<uint64_t>>& previous_from_counts() const;
    // which doubling step during which the size limit was breached
    size_t doubling_step() const;
private:
    friend class GESA;
    GESASizeException(const PathGraphSizeException& ex,
                      const std::vector<uint16_t>& node_to_comp,
                      const std::vector<size_t>& component_ranges) noexcept;
    
    std::vector<std::vector<uint64_t>> curr_counts;
    std::vector<std::vector<uint64_t>> prev_counts;
    size_t step;
};


/*
 * Template and inline implementations
 */

template<class BGraph>
GESA::GESA(const std::vector<const BGraph*>& graphs, size_t size_limit) :
    GESA(graphs.data(), graphs.size(), std::vector<const std::vector<uint64_t>*>(graphs.size(), nullptr).data(), size_limit) {
    // nothing to do besides dispatch
}

template<class BGraph>
GESA::GESA(const BGraph& graph, size_t size_limit) :
    GESA(&(&graph), 1, std::vector<const std::vector<uint64_t>*>(1, nullptr).data(), size_limit) {
    // nothing to do besides dispatch
}

template<class BGraph>
GESA::GESA(const BGraph& graph, const std::vector<uint64_t>& back_translation, size_t size_limit) :
    GESA(&(&graph), 1, std::vector<const std::vector<uint64_t>*>(1, &back_translation).data(), size_limit) {
    // nothing to do besides dispatch
}

template<class BGraph>
GESA::GESA(const std::vector<const BGraph*>& graphs,
           const std::vector<const std::vector<uint64_t>*>& back_translations, size_t size_limit) :
    GESA(graphs.data(), graphs.size(), back_translations.data(), size_limit) {
    // nothing to do besides dispatch
}

template<class BGraph>
GESA::GESA(const BGraph* const* const graphs, size_t num_graphs,
           const std::vector<uint64_t>* const* const back_translations, size_t size_limit)
{
    
    // the node IDs of the i-th input graph will be in the range [i]:[i+1]
    std::vector<size_t> component_ranges(num_graphs + 1, 0);
            
    // make a single graph with all of the components
    BaseGraph joined;
    for (size_t i = 0; i < num_graphs; ++i) {
        // add a new component
        const BGraph& graph = *graphs[i];
        append_component(joined, graph);
        
        // record the end of the node range
        component_ranges[i + 1] = joined.node_size();
    }
    
    std::vector<uint16_t> node_to_comp(joined.node_size());
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
        
        try {
            path_graph = PathGraph(path_graph, size_limit);
        }
        catch (PathGraphSizeException& ex) {
            // convert path graph nodes to original nodes
            throw GESASizeException(ex, node_to_comp, component_ranges);
        }
        logging::log(logging::Debug, "Path graph reached size " + std::to_string(path_graph.node_size()));
    }
    // TODO: if i want to reduce to a prefix range sorted graph it would happen here
        
    logging::log(logging::Debug, "Ranks stabilized, constructing edges");
    
    // reassign node IDs to be ordered by prefix rank, remove redundancies, and construct edges
    path_graph.finish(joined);
    
    // map nodes in the path graph to nodes in the input graphs
    component_ranked_ids.resize(num_graphs);
    nearest_comp_rank.resize(num_graphs);
    leaf_to_comp.resize(path_graph.node_size());
    for (size_t i = 0; i < nearest_comp_rank.size(); ++i) {
        nearest_comp_rank[i].resize(path_graph.node_size());
    }
    
    logging::log(logging::Debug, "Mapping ranks to original nodes");
    
    for (uint64_t path_node_id = 0; path_node_id < path_graph.node_size(); ++path_node_id) {
        for (size_t c = 0; c < nearest_comp_rank.size(); ++c) {
            nearest_comp_rank[c][path_node_id] = component_ranked_ids[c].size();
        }
        uint64_t node_id = path_graph.from(path_node_id);
        size_t comp = node_to_comp[node_id];
        leaf_to_comp[path_node_id] = comp;
        uint64_t original_node_id = node_id - component_ranges[comp];
        if (back_translations[comp]) {
            // apply the back translation to get matches in the original graph
            original_node_id = back_translations[comp]->at(original_node_id);
        }
        component_ranked_ids[comp].push_back(original_node_id);
    }
    // add the final past-the-last entries
    for (size_t c = 0; c < nearest_comp_rank.size(); ++c) {
        nearest_comp_rank[c].push_back(component_ranked_ids[c].size());
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
    
    label_edges(path_graph.doubling_step, joined, path_graph);
    
    logging::log(logging::Debug, "Finished constructing GESA");
    
    if (debug_gesa) {
        print(std::cerr);
    }
}

inline void GESA::construct_suffix_links() {
    auto advance = [&](size_t node) -> size_t {
        const auto& node_edges = edges[node];
        if (node_edges.empty()) {
            return -1;
        }
        else {
            return node_edges.front();
        }
    };
    ESA::construct_suffix_links(advance);
}

inline char GESA::label(const SANode& node) const {
    size_t i, j;
    std::tie(i, j) = st_node_annotation_idx(node);
    return edge_label[i][j];
}

inline SANode GESA::child(const SANode& parent, char _label) const {
    for (auto& c : children(parent)) {
        if (label(c) == _label) {
            return c;
        }
    }
    // there is no corresponding child
    return SANode(-1, -1);
}

}

#endif /* centrolign_gesa_hpp */
