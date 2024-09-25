#ifndef centrolign_cactus_hpp
#define centrolign_cactus_hpp

#include <vector>
#include <cstdint>
#include <utility>
#include <tuple>

#include "centrolign/chain_cycle_graph.hpp"
#include "centrolign/modify_graph.hpp"
#include "centrolign/compacted_graph.hpp"
#include "centrolign/adjacency_graph.hpp"
#include "centrolign/utility.hpp"

namespace centrolign {

/**
 * The cactus graph for a sequence graph
 */
template<class BGraph>
class CactusGraph {
public:
    
    // construct with a top-level chain
    CactusGraph(const BGraph& graph, const SentinelTableau& tableau);
    CactusGraph() = default;
    ~CactusGraph() = default;
    
    // standard digraph interface
    inline size_t node_size() const;
    inline const std::vector<uint64_t>& next(uint64_t node_id) const;
    inline const std::vector<uint64_t>& previous(uint64_t node_id) const;
    inline size_t next_size(uint64_t node_id) const;
    inline size_t previous_size(uint64_t node_id) const;
    
    // the cactus node corresponding to the artificial link between the sink and source nodes
    inline uint64_t get_origin() const;
    // get the sequence of nodes in the original underlying graph (requires original graph to still be held in memory)
    inline std::vector<uint64_t> next_edge_label(uint64_t node_id, size_t edge_idx) const;
    inline std::vector<uint64_t> previous_edge_label(uint64_t node_id, size_t edge_idx) const;
    // get the length of nodes in the original underlying graph
    inline size_t next_edge_label_size(uint64_t node_id, size_t edge_idx) const;
    inline size_t previous_edge_label_size(uint64_t node_id, size_t edge_idx) const;
    // get the edge index of this edge in the opposing adjacency list
    inline size_t next_reverse_edge_index(uint64_t node_id, size_t edge_idx) const;
    inline size_t previous_reverse_edge_index(uint64_t node_id, size_t edge_idx) const;
    
private:
    
    struct Node {
        Node() = default;
        ~Node() = default;
        std::vector<uint64_t> previous;
        std::vector<uint64_t> next;
        // records of (node ID in adj graph, index of adj graph edge, paired edge in previous list)
        std::vector<std::tuple<uint64_t, size_t, size_t>> next_origin;
        // paired edge in next list
        std::vector<size_t> previous_origin;
     };
    
    inline uint64_t edge_to_compacted_id(uint64_t node_id, bool next, size_t edge_idx) const;
    
    inline std::vector<uint64_t> walk_compacted_node(uint64_t compacted_id) const;
    
    const BGraph* underlying_graph = nullptr;
    
    CompactedGraph compacted_graph;
    
    AdjacencyGraph adjacency_graph;
    
    uint64_t origin = -1;
    
    std::vector<Node> nodes;
};

/**
 * Cactus tree of a cactus graph, in which each simple cycle isepresented by a special node.
 * The tree is rooted at a chain that includes the source and sink nodes, and the edgs are
 * all oriented away from the root
 */
class CactusTree {
public:
    
    template<class BGraph>
    CactusTree(const CactusGraph<BGraph>& cactus_graph);
    CactusTree() = default;
    ~CactusTree() = default;
    
    // the root chain of the tree
    uint64_t get_root() const;
    // is this node an adjacency component or a chain?
    inline bool is_chain_node(uint64_t node_id) const;
    // get the chain associated with this node as a sequence of edges of the form (node ID, is outgoing, edge index)
    // edges are oriented consistently so that taking only the node IDs from these edge records also defines the cycle
    inline const std::vector<std::tuple<uint64_t, bool, size_t>>& chain(uint64_t node_id) const;
    // if this node is not a chain, get the cactus graph node ID for this node
    inline uint64_t label(uint64_t node_id) const;
    
    inline size_t node_size() const;
    inline uint64_t get_parent(uint64_t node_id) const;
    inline const std::vector<uint64_t>& get_children(uint64_t node_id) const;
    inline const std::vector<uint64_t>& next(uint64_t node_id) const;
    inline std::vector<uint64_t> previous(uint64_t node_id) const;
    inline size_t next_size(uint64_t node_id) const;
    inline size_t previous_size(uint64_t node_id) const;
    
private:
    
    struct Node {
        
        Node() = default;
        ~Node() = default;
        
        std::vector<std::tuple<uint64_t, bool, size_t>> cycle_edges;
        std::vector<uint64_t> children;
        uint64_t parent = -1;
    };
    
    uint64_t root = -1;
    
    std::vector<Node> nodes;
};



/**
 * Template and inline implementations
 */

template<class BGraph>
CactusGraph<BGraph>::CactusGraph(const BGraph& graph, const SentinelTableau& tableau) : compacted_graph(graph) {
    static const bool debug = false;
    // make a graph of non-branching paths (in the initializer)
    
    // remember the underlying graph
    underlying_graph = &graph;
        
    // find the source and sink nodes again
    uint64_t compacted_src_id = -1, compacted_snk_id = -1;
    for (uint64_t node_id = 0; node_id < compacted_graph.node_size(); ++node_id) {
        if (compacted_graph.front(node_id) == tableau.src_id) {
            assert(compacted_src_id == -1);
            assert(compacted_graph.previous_size(node_id) == 0);
            compacted_src_id = node_id;
        }
        if (compacted_graph.back(node_id) == tableau.snk_id) {
            assert(compacted_snk_id == -1);
            assert(compacted_graph.next_size(node_id) == 0);
            compacted_snk_id = node_id;
        }
    }
    
    // join the source and sink into a cycle
    ChainCycleGraph<CompactedGraph> chain_cycle_graph(compacted_graph, compacted_src_id, compacted_snk_id);
    
    // get the graph over adjacency components
    adjacency_graph = std::move(AdjacencyGraph(chain_cycle_graph));
    
    std::vector<size_t> node_to_comp(adjacency_graph.node_size());
    {
        auto three_edge_comps = three_edge_connected_components(adjacency_graph);
        
        for (size_t i = 0; i < three_edge_comps.size(); ++i) {
            for (auto node_id : three_edge_comps[i]) {
                node_to_comp[node_id] = i;
            }
        }
        nodes.resize(three_edge_comps.size());
    }
    
    for (uint64_t node_id = 0; node_id < adjacency_graph.node_size(); ++node_id) {
        uint64_t comp1 = node_to_comp[node_id];
        const auto& next_edges = adjacency_graph.next_edges(node_id);
        for (size_t i = 0; i < next_edges.size(); ++i) {
            uint64_t comp2 = node_to_comp[next_edges[i].target];
            nodes[comp1].next.emplace_back(comp2);
            nodes[comp1].next_origin.emplace_back(node_id, i, nodes[comp2].previous.size());
            nodes[comp2].previous.emplace_back(comp1);
            nodes[comp2].previous_origin.push_back(nodes[comp1].next.size() - 1);
            if (next_edges[i].label == compacted_src_id) {
                // this is the backdoor adjacency
                assert(origin == -1);
                origin = comp1;
            }
        }
    }
    
    if (origin == -1) {
        assert(node_size() == 0);
    }
    else {
        assert(previous_size(origin) == 1);
        assert(next_size(origin) == 1);
    }
    
    if (debug) {
        std::cerr << "compacted graph\n";
        print_topology(compacted_graph, std::cerr);
        std::cerr << "cyclized graph\n";
        print_topology(chain_cycle_graph, std::cerr);
        std::cerr << "adjacency graph\n";
        print_topology(adjacency_graph, std::cerr);
        std::cerr << "cactus graph\n";
        print_topology(*this, std::cerr);
        std::cerr << "origin: " << origin << '\n';
    }
}

template<class BGraph>
inline uint64_t CactusGraph<BGraph>::get_origin() const {
    return origin;
}


template<class BGraph>
inline uint64_t CactusGraph<BGraph>::edge_to_compacted_id(uint64_t node_id, bool next, size_t edge_idx) const {
    if (!next) {
        // jump to the adjacent node, which actually stores this data
        const auto& node = nodes[node_id];
        node_id = node.previous[edge_idx];
        edge_idx = node.previous_origin[edge_idx];
    }
    const auto& adj_edge = nodes[node_id].next_origin[edge_idx];
    return adjacency_graph.next_edges(std::get<0>(adj_edge))[std::get<1>(adj_edge)].label;
}

template<class BGraph>
inline std::vector<uint64_t> CactusGraph<BGraph>::previous_edge_label(uint64_t node_id, size_t edge_idx) const {
    return walk_compacted_node(edge_to_compacted_id(node_id, false, edge_idx));
}

template<class BGraph>
inline size_t CactusGraph<BGraph>::previous_edge_label_size(uint64_t node_id, size_t edge_idx) const {
    return compacted_graph.label_size(edge_to_compacted_id(node_id, false, edge_idx));
}

template<class BGraph>
inline std::vector<uint64_t> CactusGraph<BGraph>::next_edge_label(uint64_t node_id, size_t edge_idx) const {
    return walk_compacted_node(edge_to_compacted_id(node_id, true, edge_idx));
}

template<class BGraph>
inline size_t CactusGraph<BGraph>::next_edge_label_size(uint64_t node_id, size_t edge_idx) const {
    return compacted_graph.label_size(edge_to_compacted_id(node_id, true, edge_idx));
}

template<class BGraph>
inline size_t CactusGraph<BGraph>::next_reverse_edge_index(uint64_t node_id, size_t edge_idx) const {
    return std::get<2>(nodes[node_id].next_origin[edge_idx]);
}

template<class BGraph>
inline size_t CactusGraph<BGraph>::previous_reverse_edge_index(uint64_t node_id, size_t edge_idx) const {
    return nodes[node_id].previous_origin[edge_idx];
}

template<class BGraph>
inline std::vector<uint64_t> CactusGraph<BGraph>::walk_compacted_node(uint64_t compacted_id) const {
    std::vector<uint64_t> path;
    path.reserve(compacted_graph.label_size(compacted_id));
    path.emplace_back(compacted_graph.front(compacted_id));
    while (path.back() != compacted_graph.back(compacted_id)) {
        path.push_back(underlying_graph->next(path.back()).front());
    }
    return path;
}

template<class BGraph>
inline size_t CactusGraph<BGraph>::node_size() const {
    return nodes.size();
}

template<class BGraph>
inline const std::vector<uint64_t>& CactusGraph<BGraph>::next(uint64_t node_id) const {
    return nodes[node_id].next;
}

template<class BGraph>
inline const std::vector<uint64_t>& CactusGraph<BGraph>::previous(uint64_t node_id) const {
    return nodes[node_id].previous;
}

template<class BGraph>
inline size_t CactusGraph<BGraph>::next_size(uint64_t node_id) const {
    return nodes[node_id].next.size();
}

template<class BGraph>
inline size_t CactusGraph<BGraph>::previous_size(uint64_t node_id) const {
    return nodes[node_id].previous.size();
}


template<class BGraph>
CactusTree::CactusTree(const CactusGraph<BGraph>& cactus_graph) {

    static const bool debug = false;
    
    std::vector<std::vector<std::tuple<uint64_t, bool, size_t>>> cycles;
    
    {
        // check whether a node's edge are on the stack
        std::vector<bool> stacked(cactus_graph.node_size(), false);
        // check whether an edge has been added to the stack
        std::vector<std::vector<bool>> edge_traversed(cactus_graph.node_size());
        for (uint64_t node_id = 0; node_id < cactus_graph.node_size(); ++node_id) {
            edge_traversed[node_id].resize(cactus_graph.next_size(node_id), false);
        }

        // DFS stack records of ((target ID, is next, index in list), next edge idx)
        // each record corresponds the edge list of the node selected in the previous stack frame
        std::vector<std::pair<std::vector<std::tuple<uint64_t, bool, size_t>>, size_t>> stack;
        // initialize the stack with an artificial edge the "backdoor" adjacency
        stack.emplace_back();
        stack.front().first.emplace_back(cactus_graph.get_origin(), false, -1);
        stack.front().second = 0;
        while (!stack.empty()) {
            
            auto& top = stack.back();
            if (top.second == top.first.size()) {
                // we have traversed all of this node's edges
                stack.pop_back();
                if (debug) {
                    std::cerr << "completed edge list, popping stack frame\n";
                }
                continue;
            }

            auto next_edge = top.first[top.second++];
            uint64_t next_id = std::get<0>(next_edge);
            
            if (debug) {
                std::cerr << "stack state:\n";
                for (size_t i = 0; i < stack.size(); ++i) {
                    std::cerr << "frame " << i << ":\n";
                    auto& frame = stack[i];
                    for (size_t j = 0; j < frame.first.size(); ++j) {
                        std::cerr << '\t';
                        if (j + 1 == frame.second) {
                            std::cerr << '*';
                        }
                        std::cerr << std::get<0>(frame.first[j]) << ' ' << std::get<1>(frame.first[j]) << ' ' << (int64_t) std::get<2>(frame.first[j]) << '\n';
                    }
                }
            }
            
            if (stack.size() != 1) {
                // we did not get here from the virtual edge to the origin
                
                // we can get the head of this edge from the previous stack frame
                uint64_t prev_id = std::get<0>(stack[stack.size() - 2].first[stack[stack.size() - 2].second - 1]);
                
                // get the identifier of the edge to here its "next" edge list
                uint64_t edge_src_id;
                size_t edge_idx;
                if (std::get<1>(next_edge)) {
                    edge_src_id = prev_id;
                    edge_idx = std::get<2>(next_edge);
                }
                else {
                    edge_src_id = next_id;
                    edge_idx = cactus_graph.previous_reverse_edge_index(prev_id, std::get<2>(next_edge));
                }
                
                if (edge_traversed[edge_src_id][edge_idx]) {
                    // we already traversed this edge in the other direction
                    if (debug) {
                        std::cerr << "edge has already been traversed in the other direction\n";
                    }
                    continue;
                }
                
                edge_traversed[edge_src_id][edge_idx] = true;
            }
            if (!stacked[next_id]) {
                // we haven't come to this node before
                if (debug) {
                    std::cerr << "encountering node " << next_id << " for the first time\n";
                }

                // create a new stack entry
                stack.emplace_back(std::vector<std::tuple<uint64_t, bool, size_t>>(), 0);
                auto& new_stack_record = stack.back();
                
                // get the edges
                for (bool next : {false, true}) {
                    const auto& edges = next ? cactus_graph.next(next_id) : cactus_graph.previous(next_id);
                    for (size_t i = 0; i < edges.size(); ++i) {
                        new_stack_record.first.emplace_back(edges[i], next, i);
                    }
                }
                stacked[next_id] = true;
            }
            else {
                // we've hit a cycle, which must be a simple cycle because of the properties of cactus graphs and DFS

                if (debug) {
                    std::cerr << "found cycle to " << next_id << " over edge " << std::get<0>(next_edge) << ' ' << std::get<1>(next_edge) << ' ' << (int64_t) std::get<2>(next_edge) << '\n';
                }
                
                cycles.emplace_back();
                auto& cycle = cycles.back();
                // walk the stack records back to get the rest of the edges
                size_t i = stack.size() - 1;
                while (true) {
                    const auto& curr_edge = stack[i].first[stack[i].second - 1];
                    const auto& prev_edge = stack[i - 1].first[stack[i - 1].second - 1];
                    cycle.emplace_back(std::get<0>(prev_edge), std::get<1>(curr_edge), std::get<2>(curr_edge));
                    if (std::get<0>(prev_edge) == next_id) {
                        break;
                    }
                    --i;
                }
                // the cycle was constructed in reverse
                std::reverse(cycle.begin(), cycle.end());
                if (debug) {
                    std::cerr << "derived cycle:\n";
                    for (size_t i = 0; i < cycle.size(); ++i) {
                        auto& edge = cycle[i];
                        if (i) {
                            std::cerr << ',';
                        }
                        std::cerr << ' ' << std::get<0>(edge) << ' ' << std::get<1>(edge) << ' ' << (int64_t) std::get<2>(edge);
                    }
                    std::cerr << '\n';
                }
            }
        }
    }

    // we should at least have the backdoor cycle
    assert(cycles.size() >= 1);
    
    if (debug) {
        std::cerr << "final cycles:\n";
        for (const auto& cycle : cycles) {
            for (size_t i = 0; i < cycle.size(); ++i) {
                auto& edge = cycle[i];
                if (i) {
                    std::cerr << ',';
                }
                std::cerr << ' ' << std::get<0>(edge) << ' ' << std::get<1>(edge) << ' ' << (int64_t) std::get<2>(edge);
            }
            std::cerr << '\n';
        }
    }

    // make a lookup for edge to cycle (if any), also find the root
    std::vector<std::vector<uint64_t>> assigned_cycle(cactus_graph.node_size());
    for (uint64_t node_id = 0; node_id < cactus_graph.node_size(); ++node_id) {
        assigned_cycle[node_id].resize(cactus_graph.next_size(node_id), -1);
    }
    size_t root_cycle = -1;
    for (size_t i = 0; i < cycles.size(); ++i) {
        const auto& cycle = cycles[i];
        if (std::get<0>(cycle.front()) == cactus_graph.get_origin()) {
            root_cycle = i;
        }
        for (const auto& cycle_edge : cycle) {
            // get the identifier of the edge as a "next"
            uint64_t node_id;
            size_t edge_idx;
            if (std::get<1>(cycle_edge)) {
                node_id = std::get<0>(cycle_edge);
                edge_idx = std::get<2>(cycle_edge);
            }
            else {
                node_id = cactus_graph.previous(std::get<0>(cycle_edge))[std::get<2>(cycle_edge)];
                edge_idx = cactus_graph.previous_reverse_edge_index(std::get<0>(cycle_edge), std::get<2>(cycle_edge));
            }
            assigned_cycle[node_id][edge_idx] = i;
        }
    }
    
    assert(root_cycle != -1 || cactus_graph.node_size() == 0);

    // initialize the cactus tree nodes
    nodes.resize(cactus_graph.node_size() + cycles.size());
    for (size_t i = 0, j = cactus_graph.node_size(); i < cycles.size(); ++i, ++j) {
        nodes[j].cycle_edges = std::move(cycles[i]);
    }

    {
        std::vector<bool> stacked(nodes.size(), false);
        // now do a DFS from the root cycle node to orient the edges downward
        std::vector<uint64_t> stack;

        root = cactus_graph.node_size() + root_cycle;
        stack.push_back(root);
        stacked[root] = true;
        
        if (debug) {
            std::cerr << "cactus tree has " << nodes.size() << " nodes, beginning edge construction traversal at root " << root << '\n';
        }

        while (!stack.empty()) {

            if (debug) {
                std::cerr << "stack state\n";
                for (auto n : stack) {
                    std::cerr << '\t' << n << '\n';
                }
            }
            
            uint64_t node_id = stack.back();
            stack.pop_back();
            
            if (debug) {
                std::cerr << "pop node " << node_id << " from stack\n";
            }

            if (node_id >= cactus_graph.node_size()) {
                // this is a cycle node
                if (debug) {
                    std::cerr << "cycle node\n";
                }
                for (const auto& edge : nodes[node_id].cycle_edges) {
                    // get the node on the counterclockwise side of the edge
                    uint64_t next_id = std::get<0>(edge);
                    if (debug) {
                        std::cerr << "\tnode " << next_id << " from edge " << std::get<0>(edge) << ' ' << std::get<1>(edge) << ' ' << std::get<2>(edge) << '\n';
                    }
                    if (stacked[next_id]) {
                        continue;
                    }
                    nodes[node_id].children.push_back(next_id);
                    nodes[next_id].parent = node_id;
                    stack.push_back(next_id);
                    stacked[next_id] = true;
                }
            }
            else {
                // this is an adjacency component node
                if (debug) {
                    std::cerr << "adjacency component node\n";
                }
                for (bool next : {true, false}) {
                    const auto& edges = next ? cactus_graph.next(node_id) : cactus_graph.previous(node_id);
                    for (size_t i = 0; i < edges.size(); ++i) {
                        // get the edge identifier as a "next"
                        uint64_t edge_node_id;
                        uint64_t edge_idx;
                        if (next) {
                            edge_node_id = node_id;
                            edge_idx = i;
                        }
                        else {
                            edge_node_id = edges[i];
                            edge_idx = cactus_graph.previous_reverse_edge_index(node_id, i);
                        }
                        auto cycle = assigned_cycle[edge_node_id][edge_idx];
                        
                        uint64_t next_id;
                        if (cycle == -1) {
                            // this edge is across a bridge, edge connects them directly
                            next_id = edges[i];
                        }
                        else {
                            // this edge is in a cycle, edge goes to cycle node
                            next_id = cactus_graph.node_size() + cycle;
                        }
                        
                        if (debug) {
                            std::cerr << '\t' << "edge " << next << ' ' << i << " to " << next_id << '\n';
                        }
                        
                        if (stacked[next_id]) {
                            continue;
                        }
                        nodes[node_id].children.push_back(next_id);
                        nodes[next_id].parent = node_id;
                        stack.push_back(next_id);
                        stacked[next_id] = true;
                    }
                }
            }
        }
    }
    
    if (debug) {
        std::cerr << "cactus tree\n";
        print_topology(*this, std::cerr);
    }
}



uint64_t CactusTree::get_root() const {
    return root;
}

inline size_t CactusTree::node_size() const {
    return nodes.size();
}

inline uint64_t CactusTree::get_parent(uint64_t node_id) const {
    return nodes[node_id].parent;
}

inline const std::vector<uint64_t>& CactusTree::get_children(uint64_t node_id) const {
    return nodes[node_id].children;
}

inline bool CactusTree::is_chain_node(uint64_t node_id) const {
    return !nodes[node_id].cycle_edges.empty();
}

inline const std::vector<std::tuple<uint64_t, bool, size_t>>& CactusTree::chain(uint64_t node_id) const {
    return nodes[node_id].cycle_edges;
}

inline uint64_t CactusTree::label(uint64_t node_id) const {
    return is_chain_node(node_id) ? -1 : node_id;
}

inline const std::vector<uint64_t>& CactusTree::next(uint64_t node_id) const {
    return get_children(node_id);
}

inline std::vector<uint64_t> CactusTree::previous(uint64_t node_id) const {
    return std::vector<uint64_t>(1, get_parent(node_id));
}

inline size_t CactusTree::next_size(uint64_t node_id) const {
    return get_children(node_id).size();
}

inline size_t CactusTree::previous_size(uint64_t node_id) const {
    return node_id == root ? 0 : 1;
}



}

#endif /* centrolign_cactus_hpp */
