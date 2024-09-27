#ifndef centrolign_structure_tree_hpp
#define centrolign_structure_tree_hpp

#include <vector>
#include <utility>
#include <cstdint>
#include <limits>
#include <iostream>
#include <stdexcept>
#include <unordered_map>
#include <algorithm>

#include "centrolign/topological_order.hpp"
#include "centrolign/modify_graph.hpp"
#include "centrolign/utility.hpp"

namespace centrolign {

/*
 * Data structure provides navigation methods for chainable two-sided graph features
 */
class TwoDisconnectedStructureTree {
public:
    
    TwoDisconnectedStructureTree() = default;
    TwoDisconnectedStructureTree(const TwoDisconnectedStructureTree& other) = default;
    TwoDisconnectedStructureTree(TwoDisconnectedStructureTree&& tree) = default;
    ~TwoDisconnectedStructureTree() = default;
    TwoDisconnectedStructureTree& operator=(const TwoDisconnectedStructureTree& other) = default;
    TwoDisconnectedStructureTree& operator=(TwoDisconnectedStructureTree&& other) = default;
    
    // number of chains
    inline size_t chain_size() const;
    // number of 2-disconnected structures
    inline size_t structure_size() const;
    
    // access structures by graph node ID (returns -1 is there is none)
    inline uint64_t structure_beginning_at(uint64_t node_id) const;
    inline uint64_t structure_ending_at(uint64_t node_id) const;
    
    // returns pair of graph node IDs
    inline const std::pair<uint64_t, uint64_t>& structure_boundaries(uint64_t struct_id) const;
    
    // returns vector of chain IDs, in no particular order
    inline const std::vector<uint64_t>& chains_inside(uint64_t struct_id) const;
    
    // returns chain ID
    inline uint64_t chain_containing(uint64_t struct_id) const;
    
    // returns vector of structure IDs, in order from start to end
    inline const std::vector<uint64_t>& structures_inside(uint64_t chain_id) const;
    
    // returns structure ID (returns -1 is there is none)
    inline uint64_t structure_containing(uint64_t chain_id) const;
    
    // returns all nodes as (feature ID, feature is chain) pairs
    std::vector<std::pair<uint64_t, bool>> postorder() const;
    
protected:
    
    // helper function for derived constructors
    template<class Derived, class Graph>
    void initialize(const Graph& graph, const SentinelTableau* tableau);
    
    // CRTP method used to find 2-disconnected graph features in derived classes
    template<class Derived, class Graph>
    std::vector<std::pair<uint64_t, uint64_t>> find_2_disc_structures(const Graph& graph, const SentinelTableau* tableau);
    
    struct Chain {
        Chain() = default;
        ~Chain() = default;
        
        std::vector<uint64_t> structure_ids;
        uint64_t parent = -1;
    };
    
    struct TwoDisconnectedStructure {
        TwoDisconnectedStructure(const std::pair<uint64_t, uint64_t>& structure) : boundaries(structure) {}
        TwoDisconnectedStructure() = default;
        ~TwoDisconnectedStructure() = default;
        
        std::pair<uint64_t, uint64_t> boundaries;
        uint64_t parent = -1;
        std::vector<uint64_t> chain_ids;
    };
    
    std::vector<TwoDisconnectedStructure> structures;
    std::vector<Chain> chains;
    std::vector<uint64_t> structure_endings;
    std::vector<uint64_t> structure_beginnings;
};

/*
 * A graph for a 2-disconnected structure where contained chains have been abstracted
 * to single nodes
 */
class NetGraph {
public:
    
    template<class Graph>
    NetGraph(const Graph& graph, const TwoDisconnectedStructureTree& structures,
             uint64_t struct_id);
    
    // make a net graph for the outer contents not contained in any 2-disconnected structure
    // note: linear in size of graph to construct
    template<class Graph>
    NetGraph(const Graph& graph, const TwoDisconnectedStructureTree& structures);
    
    template<class Graph>
    NetGraph(const Graph& graph, const TwoDisconnectedStructureTree& structures,
             const SentinelTableau& tableau);
    
    NetGraph() = default;
    ~NetGraph() = default;
    
    
    // pair indicating the feature ID, and whether it is a chain (false => it is a node)
    inline std::pair<uint64_t, bool> label(uint64_t node_id) const;
    
    inline size_t node_size() const;
    inline const std::vector<uint64_t>& next(uint64_t node_id) const;
    inline size_t next_size(uint64_t node_id) const;
    inline size_t previous_size(uint64_t node_id) const;
    
private:
    
    template<class Graph>
    NetGraph(const Graph& graph, const TwoDisconnectedStructureTree& structures,
             const SentinelTableau* tableau);
    
    struct Node {
        Node() = default;
        ~Node() = default;

        std::pair<uint64_t, bool> feature_id;
        std::vector<uint64_t> edges;
        size_t in_degree = 0; // for topo sort interface
    };
    
    inline uint64_t add_node(uint64_t feature_id, bool is_chain);
    inline void add_edge(uint64_t node_id1, uint64_t node_id2);
    
    std::vector<Node> nodes;
    
    void print(std::ostream& out) const;
};




/*
 * Template implementations
 */

template<class Derived, class Graph>
void TwoDisconnectedStructureTree::initialize(const Graph& graph, const SentinelTableau* tableau) {
    
    static const bool debug = false;
    
    structure_beginnings.resize(graph.node_size(), -1);
    structure_endings.resize(graph.node_size(), -1);
    
    // identify and collect the 2-disconnected structures
    for (const auto& structure : find_2_disc_structures<Derived, Graph>(graph, tableau)) {
        if (tableau && (tableau->src_id == structure.first || tableau->snk_id == structure.second
                        || tableau->snk_id == structure.first || tableau->src_id == structure.second)) {
            // structures including the sentinels are not interesting to us
            continue;
        }
        
        if (debug) {
            std::cerr << "recording bubble " << structure.first << " " << structure.second << '\n';
        }
        
        structure_beginnings[structure.first] = structures.size();
        structure_endings[structure.second] = structures.size();
        
        structures.emplace_back(structure);
    }
    
    // form them into chains
    for (uint64_t struct_id = 0; struct_id < structures.size(); ++struct_id) {
        if (structures[struct_id].parent != -1) {
            // we already chained this one
            continue;
        }
        
        uint64_t chain_id = chains.size();
        chains.emplace_back();
        auto& chain = chains.back();
        
        // extend to the left
        chain.structure_ids.push_back(struct_id);
        structures[struct_id].parent = chain_id;
        uint64_t struct_id_here = structure_ending_at(structure_boundaries(struct_id).first);
        while (struct_id_here != -1) {
            chain.structure_ids.push_back(struct_id_here);
            structures[struct_id_here].parent = chain_id;
            struct_id_here = structure_ending_at(structure_boundaries(struct_id_here).first);
        }
        // we build this backwards, now reverse it
        std::reverse(chain.structure_ids.begin(), chain.structure_ids.end());
        // extend to the right
        struct_id_here = structure_beginning_at(structure_boundaries(struct_id).second);
        while (struct_id_here != -1) {
            chain.structure_ids.push_back(struct_id_here);
            structures[struct_id_here].parent = chain_id;
            struct_id_here = structure_beginning_at(structure_boundaries(struct_id_here).second);
        }
    }
    
    // identify the chains within structures
    
    // we'll use a global bank for the traversed status so that we don't need to use
    // hash sets in the inner loop
    std::vector<bool> traversed(graph.node_size(), false);
    for (uint64_t struct_id = 0; struct_id < structures.size(); ++struct_id) {
        
        auto& structure = structures[struct_id];
                
        // DFS to find contained chains
        std::vector<uint64_t> stack(1, structure.boundaries.first);
        while (!stack.empty()) {
            auto node_id = stack.back();
            stack.pop_back();
            for (auto next_id : graph.next(node_id)) {
                if (next_id == structure.boundaries.second || traversed[next_id]) {
                    // don't go into the end of the structure
                    continue;
                }
                traversed[next_id] = true;
                auto next_struct_id = structure_beginning_at(next_id);
                if (next_struct_id != -1) {
                    // record this as a contained chain
                    auto chain_id = chain_containing(next_struct_id);
                    chains[chain_id].parent = struct_id;
                    structure.chain_ids.push_back(chain_id);

                    // jump to the end of the chain
                    auto final_struct_id = structures_inside(chain_id).back();
                    auto final_node_id = structure_boundaries(final_struct_id).second;
                    traversed[final_node_id] = true;
                    stack.push_back(final_node_id);
                }
                else {
                    // just a regular node
                    stack.push_back(next_id);
                }
            }
        }
    }
}

inline size_t TwoDisconnectedStructureTree::chain_size() const {
    return chains.size();
}

inline size_t TwoDisconnectedStructureTree::structure_size() const {
    return structures.size();
}

inline uint64_t TwoDisconnectedStructureTree::structure_ending_at(uint64_t node_id) const {
    return structure_endings[node_id];
}

inline uint64_t TwoDisconnectedStructureTree::structure_beginning_at(uint64_t node_id) const {
    return structure_beginnings[node_id];
}

inline const std::pair<uint64_t, uint64_t>& TwoDisconnectedStructureTree::structure_boundaries(uint64_t struct_id) const {
    return structures[struct_id].boundaries;
}

inline const std::vector<uint64_t>& TwoDisconnectedStructureTree::chains_inside(uint64_t struct_id) const {
    return structures[struct_id].chain_ids;
}

inline uint64_t TwoDisconnectedStructureTree::chain_containing(uint64_t struct_id) const {
    return structures[struct_id].parent;
}

inline const std::vector<uint64_t>& TwoDisconnectedStructureTree::structures_inside(uint64_t chain_id) const {
    return chains[chain_id].structure_ids;
}

inline uint64_t TwoDisconnectedStructureTree::structure_containing(uint64_t chain_id) const {
    return chains[chain_id].parent;
}


template<class Derived, class Graph>
std::vector<std::pair<uint64_t, uint64_t>> TwoDisconnectedStructureTree::find_2_disc_structures(const Graph& graph,
                                                                                                const SentinelTableau* tableau) {
    return static_cast<Derived*>(this)->find_2_disc_structures_impl(graph, tableau);
}



template<class Graph>
NetGraph::NetGraph(const Graph& graph, const TwoDisconnectedStructureTree& structures,
                   uint64_t struct_id) {
        
    static const bool debug = false;
    
    uint64_t start, end;
    std::tie(start, end) = structures.structure_boundaries(struct_id);
    
    if (debug) {
        std::cerr << "building netgraph for superbubble " << struct_id << " with boundaries " << start << ", " << end << '\n';
    }
    
    std::unordered_map<uint64_t, uint64_t> forward_translation;
    forward_translation[start] = add_node(start, false);
    
    // DFS to find contained chains
    std::vector<uint64_t> stack(1, start);
    while (!stack.empty()) {
        auto node_id = stack.back();
        stack.pop_back();
        if (debug) {
            std::cerr << "dequeue " << node_id << '\n';
        }
        
        if (node_id == end) {
            // we stop going at the end of the bubble
            continue;
        }
        
        for (auto next_id : graph.next(node_id)) {
            if (debug) {
                std::cerr << "follow edge to " << next_id << '\n';
            }
            auto it = forward_translation.find(next_id);
            if (it != forward_translation.end()) {
                // already traversed and initialized
                add_edge(forward_translation[node_id], it->second);
                if (debug) {
                    std::cerr << "add edge between net nodes " << forward_translation[node_id] << " and " << it->second << '\n';
                }
            }
            else {
                auto next_struct_id = structures.structure_beginning_at(next_id);
                if (next_struct_id != -1 && next_id != end) {
                    // new chain node
                    auto chain_id = structures.chain_containing(next_struct_id);
                    auto chain_net_id = add_node(chain_id, true);
                    
                    // jump to the end of the chain
                    auto final_struct_id = structures.structures_inside(chain_id).back();
                    auto final_node_id = structures.structure_boundaries(final_struct_id).second;
                    
                    // both ends of the chain project to this net graph node
                    forward_translation[next_id] = chain_net_id;
                    forward_translation[final_node_id] = chain_net_id;
                    
                    add_edge(forward_translation[node_id], chain_net_id);
                    
                    stack.push_back(final_node_id);
                    if (debug) {
                        std::cerr << "add new net node " << chain_net_id << " from chain " << chain_id << '\n';
                    }
                }
                else {
                    // new non-chain node
                    auto net_id = add_node(next_id, false);
                    forward_translation[next_id] = net_id;
                    
                    add_edge(forward_translation[node_id], net_id);
                    
                    stack.push_back(next_id);
                    if (debug) {
                        std::cerr << "add new net node " << net_id << '\n';
                    }
                }
            }
        }
    }
    
    if (debug) {
        std::cerr << "final topology for bubble " << struct_id << " net graph" << "\n";
        print(std::cerr);
    }
}

template<class Graph>
NetGraph::NetGraph(const Graph& graph, const TwoDisconnectedStructureTree& structures) : NetGraph(graph, structures, nullptr) {
    
}

template<class Graph>
NetGraph::NetGraph(const Graph& graph, const TwoDisconnectedStructureTree& structures,
                   const SentinelTableau& tableau) : NetGraph(graph, structures, &tableau) {
    
}

template<class Graph>
NetGraph::NetGraph(const Graph& graph, const TwoDisconnectedStructureTree& structures,
                   const SentinelTableau* tableau) {
    
    static const bool debug = false;
    
    // mark all the nodes that are contained in a superbubble
    std::vector<bool> contained(graph.node_size());
    for (uint64_t struct_id = 0; struct_id < structures.structure_size(); ++struct_id) {
        NetGraph net_graph(graph, structures, struct_id);
        for (uint64_t net_id = 0; net_id < net_graph.node_size(); ++net_id) {
            uint64_t feature_id;
            bool is_chain;
            std::tie(feature_id, is_chain) = net_graph.label(net_id);
            if (!is_chain) {
                contained[feature_id] = true;
                if (debug) {
                    std::cerr << feature_id << " marked contained from bubble " << struct_id << '\n';
                }
            }
        }
    }
    
    if (debug) {
        std::cerr << "making nodes\n";
    }
    
    // add chain nodes
    std::unordered_map<std::pair<uint64_t, bool>, uint64_t> forward_translation;
    for (uint64_t chain_id = 0; chain_id < structures.chain_size(); ++chain_id) {
        if (structures.structure_containing(chain_id) == -1) {
            // top level chain
            forward_translation[std::make_pair(chain_id, true)] = add_node(chain_id, true);
        }
    }
    for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
        if (!contained[node_id] && (!tableau || (node_id != tableau->src_id && node_id != tableau->snk_id))) {
            forward_translation[std::make_pair(node_id, false)] = add_node(node_id, false);
        }
    }
    
    if (debug) {
        std::cerr << "made nodes:\n";
        for (const auto& t : forward_translation) {
            std::cerr << '\t' << t.first.first << ' ' << t.first.second << ": " << t.second << '\n';
        }
        std::cerr << "making edges\n";
    }
    
    // make the edges
    for (uint64_t net_id = 0; net_id < node_size(); ++net_id) {
        uint64_t feature_id;
        bool is_chain;
        std::tie(feature_id, is_chain) = label(net_id);
        if (is_chain) {
            feature_id = structures.structure_boundaries(structures.structures_inside(feature_id).back()).second;
        }
        for (auto next_id : graph.next(feature_id)) {
            if (tableau && next_id == tableau->snk_id) {
                continue;
            }
            auto struct_id = structures.structure_beginning_at(next_id);
            uint64_t next_net_id;
            if (struct_id == -1) {
                next_net_id = forward_translation.at(std::make_pair(next_id, false));
            }
            else {
                auto chain_id = structures.chain_containing(struct_id);
                next_net_id = forward_translation.at(std::make_pair(chain_id, true));
            }
            add_edge(net_id, next_net_id);
        }
    }
    
    if (debug) {
        std::cerr << "final topology for outside-bubble net graph\n";
        print(std::cerr);
    }
}

inline uint64_t NetGraph::add_node(uint64_t feature_id, bool is_chain) {
    nodes.emplace_back();
    auto& node = nodes.back();
    node.feature_id.first = feature_id;
    node.feature_id.second = is_chain;
    return nodes.size() - 1;
}

inline void NetGraph::add_edge(uint64_t node_id1, uint64_t node_id2) {
    nodes[node_id1].edges.push_back(node_id2);
    nodes[node_id2].in_degree++;
}

inline std::pair<uint64_t, bool> NetGraph::label(uint64_t node_id) const {
    return nodes[node_id].feature_id;
}

inline size_t NetGraph::node_size() const {
    return nodes.size();
}

inline const std::vector<uint64_t>& NetGraph::next(uint64_t node_id) const {
    return nodes[node_id].edges;
}

inline size_t NetGraph::next_size(uint64_t node_id) const {
    return nodes[node_id].edges.size();
}

inline size_t NetGraph::previous_size(uint64_t node_id) const {
    return nodes[node_id].in_degree;
}

}

#endif /* centrolign_structure_tree_hpp */
