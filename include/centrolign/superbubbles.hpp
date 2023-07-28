#ifndef centrolign_superbubbles_hpp
#define centrolign_superbubbles_hpp

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
 * Data structure that finds superbubbles and provides navigation
 * methods for the superbubble/chain tree
 */
class SuperbubbleTree {
public:
    
    // construct from a single-source, single-sink graph
    template<class Graph>
    SuperbubbleTree(const Graph& graph);
    
    // construct from a graph with sentinel sink and source nodes, which
    // will be removed from the superbubble decomposition
    template<class Graph>
    SuperbubbleTree(const Graph& graph, const SentinelTableau& tableau);
    
    SuperbubbleTree() = default;
    SuperbubbleTree(const SuperbubbleTree& other) = default;
    SuperbubbleTree(SuperbubbleTree&& tree) = default;
    ~SuperbubbleTree() = default;
    SuperbubbleTree& operator=(const SuperbubbleTree& other) = default;
    SuperbubbleTree& operator=(SuperbubbleTree&& other) = default;
    
    // number of chains
    inline size_t chain_size() const;
    // number of superbubbles
    inline size_t superbubble_size() const;
    
    // access superbubbles by graph node ID (returns -1 is there is none)
    inline uint64_t superbubble_beginning_at(uint64_t node_id) const;
    inline uint64_t superbubble_ending_at(uint64_t node_id) const;
    
    // returns pair of graph node IDs
    inline const std::pair<uint64_t, uint64_t>& superbubble_boundaries(uint64_t superbubble_id) const;
    
    // returns vector of chain IDs, in no particular order
    inline const std::vector<uint64_t>& chains_inside(uint64_t superbubble_id) const;
    
    // returns chain ID
    inline uint64_t chain_containing(uint64_t superbubble_id) const;
    
    // returns vector of superbubble IDs, in order from start to end
    inline const std::vector<uint64_t>& superbubbles_inside(uint64_t chain_id) const;
    
    // returns superbubble ID (returns -1 is there is none)
    inline uint64_t superbubble_containing(uint64_t chain_id) const;
    
    // returns all nodes as (feature ID, feature is chain) pairs
    std::vector<std::pair<uint64_t, bool>> postorder() const;
    
protected:
    
    template<class Graph>
    SuperbubbleTree(const Graph& graph, const SentinelTableau* tableau);
    
    // find superbubbles with algorithm of Gartner, et al. (2018) for DAGs
    template<class Graph>
    static std::vector<std::pair<uint64_t, uint64_t>> find_superbubbles(const Graph& graph);
    
    struct Chain {
        Chain() = default;
        ~Chain() = default;
        
        std::vector<uint64_t> superbubble_ids;
        uint64_t parent = -1;
    };
    
    struct Superbubble {
        Superbubble(const std::pair<uint64_t, uint64_t>& superbubble) : boundaries(superbubble) {}
        Superbubble() = default;
        ~Superbubble() = default;
        
        std::pair<uint64_t, uint64_t> boundaries;
        uint64_t parent = -1;
        std::vector<uint64_t> chain_ids;
    };
    
    std::vector<Superbubble> superbubbles;
    std::vector<Chain> chains;
    std::vector<uint64_t> superbubble_endings;
    std::vector<uint64_t> superbubble_beginnings;
};

/*
 * A graph for a superbubble where contained chains have been abstracted
 * to single nodes
 */
class NetGraph {
public:
    
    template<class Graph>
    NetGraph(const Graph& graph, const SuperbubbleTree& superbubbles,
             uint64_t superbubble_id);
    
    // make a net graph for the outer contents not contained in any superbubble
    // note: linear in size of graph to construct
    template<class Graph>
    NetGraph(const Graph& graph, const SuperbubbleTree& superbubbles);
    template<class Graph>
    NetGraph(const Graph& graph, const SuperbubbleTree& superbubbles,
             const SentinelTableau& tableau);
    
    NetGraph() = default;
    ~NetGraph() = default;
    
    
    // pair indicating the feature ID, and whether it is a chain (false -> it is a node)
    inline std::pair<uint64_t, bool> label(uint64_t node_id) const;
    
    inline size_t node_size() const;
    inline const std::vector<uint64_t>& next(uint64_t node_id) const;
    inline size_t next_size(uint64_t node_id) const;
    inline size_t previous_size(uint64_t node_id) const;
    
private:
    
    template<class Graph>
    NetGraph(const Graph& graph, const SuperbubbleTree& superbubbles,
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


template<class Graph>
SuperbubbleTree::SuperbubbleTree(const Graph& graph) : SuperbubbleTree(graph, nullptr) {
    // nothing else
}

template<class Graph>
SuperbubbleTree::SuperbubbleTree(const Graph& graph,
                                 const SentinelTableau& tableau) : SuperbubbleTree(graph, &tableau) {
    // nothing else
}


template<class Graph>
SuperbubbleTree::SuperbubbleTree(const Graph& graph, const SentinelTableau* tableau) {
    
    static const bool debug = false;
    
    superbubble_beginnings.resize(graph.node_size(), -1);
    superbubble_endings.resize(graph.node_size(), -1);
    
    // identify and collect the superbubbles
    for (const auto& sup_bub : find_superbubbles(graph)) {
        if (tableau && (tableau->src_id == sup_bub.first || tableau->snk_id == sup_bub.second)) {
            // superbubbles including the sentinels are not interesting to us
            continue;
        }
        
        if (debug) {
            std::cerr << "recording bubble " << sup_bub.first << " " << sup_bub.second << '\n';
        }
        
        superbubble_beginnings[sup_bub.first] = superbubbles.size();
        superbubble_endings[sup_bub.second] = superbubbles.size();
        
        superbubbles.emplace_back(sup_bub);
    }
    
    // form them into chains
    for (uint64_t bub_id = 0; bub_id < superbubbles.size(); ++bub_id) {
        if (superbubbles[bub_id].parent != -1) {
            // we already chained this one
            continue;
        }
        
        uint64_t chain_id = chains.size();
        chains.emplace_back();
        auto& chain = chains.back();
        
        // extend to the left
        chain.superbubble_ids.push_back(bub_id);
        superbubbles[bub_id].parent = chain_id;
        uint64_t bub_id_here = superbubble_ending_at(superbubble_boundaries(bub_id).first);
        while (bub_id_here != -1) {
            chain.superbubble_ids.push_back(bub_id_here);
            superbubbles[bub_id_here].parent = chain_id;
            bub_id_here = superbubble_ending_at(superbubble_boundaries(bub_id_here).first);
        }
        // we build this backwards, now reverse it
        std::reverse(chain.superbubble_ids.begin(), chain.superbubble_ids.end());
        // extend to the right
        bub_id_here = superbubble_beginning_at(superbubble_boundaries(bub_id).second);
        while (bub_id_here != -1) {
            chain.superbubble_ids.push_back(bub_id_here);
            superbubbles[bub_id_here].parent = chain_id;
            bub_id_here = superbubble_beginning_at(superbubble_boundaries(bub_id_here).second);
        }
    }
    
    // identify the chains within superbubbles
    
    // we'll use a global bank for the traversed status so that we don't need to use
    // hash sets in the inner loop
    std::vector<bool> traversed(graph.node_size(), false);
    for (uint64_t bub_id = 0; bub_id < superbubbles.size(); ++bub_id) {
        
        auto& superbubble = superbubbles[bub_id];
                
        // DFS to find contained chains
        std::vector<uint64_t> stack(1, superbubble.boundaries.first);
        while (!stack.empty()) {
            auto node_id = stack.back();
            stack.pop_back();
            for (auto next_id : graph.next(node_id)) {
                if (next_id == superbubble.boundaries.second || traversed[next_id]) {
                    // don't go into the end of the superbubble
                    continue;
                }
                traversed[next_id] = true;
                auto next_bub_id = superbubble_beginning_at(next_id);
                if (next_bub_id != -1) {
                    // record this as a contained chain
                    auto chain_id = chain_containing(next_bub_id);
                    chains[chain_id].parent = bub_id;
                    superbubble.chain_ids.push_back(chain_id);

                    // jump to the end of the chain
                    auto final_bub_id = superbubbles_inside(chain_id).back();
                    auto final_node_id = superbubble_boundaries(final_bub_id).second;
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

template<class Graph>
std::vector<std::pair<uint64_t, uint64_t>> SuperbubbleTree::find_superbubbles(const Graph& graph) {
    
    static const bool debug = false;
    
    std::vector<std::pair<uint64_t, uint64_t>> superbubbles;
    
    // TODO: i'm pretty sure the depth-first Kahn's algorithm also works here,
    // so i'm not using the DFS postorder like in the paper
    auto order = topological_order(graph);
    std::vector<size_t> index(order.size());
    size_t num_sources = 0;
    size_t num_sinks = 0;
    for (size_t i = 0; i < order.size(); ++i) {
        index[order[i]] = i;
        if (graph.previous_size(i) == 0) {
            ++num_sources;
        }
        if (graph.next_size(i) == 0) {
            ++num_sinks;
        }
    }
    
    if (num_sources != 1 || num_sinks != 1) {
        throw std::runtime_error("Can only find superbubbles in single-source, single-sink graphs");
    }
    
    if (debug) {
        std::cerr << "got postorder:\n";
        for (auto n : order) {
            std::cerr << ' ' << n;
        }
        std::cerr << '\n';
    }
    
    // indexes of the candidate ends
    std::vector<int64_t> candidate_stack;
    // the furthest backward reach (filled in as we go at candidate ends)
    std::vector<int64_t> backward_reach(order.size(),
                                        std::numeric_limits<int64_t>::max());
    for (int64_t i = order.size() - 1; i >= 0; --i) {
        
        if (debug) {
            std::cerr << "at index " << i << ", node " << order[i] << ", stack:\n";
            for (auto t : candidate_stack) {
                std::cerr << ' ' << t;
            }
            std::cerr << '\n';
        }
        
        int64_t forward_reach = -1;
        for (auto node_id : graph.next(order[i])) {
            forward_reach = std::max<int64_t>(forward_reach, index[node_id]);
        }
        if (forward_reach == i + 1) {
            candidate_stack.push_back(i + 1);
            if (debug) {
                std::cerr << "identify " << (i + 1) << " (node " << order[i + 1] <<  ") as a candidate end\n";
            }
        }
        
        while (!candidate_stack.empty() && forward_reach > candidate_stack.back()) {
            // this node reaches outside the candidate interval, so the nearest candidate
            // is not a superbubble
            auto invalid_candidate = candidate_stack.back();
            candidate_stack.pop_back();
            if (!candidate_stack.empty()) {
                // give backward reach info to next nearest candidate
                backward_reach[candidate_stack.back()] = std::min(backward_reach[candidate_stack.back()],
                                                                  backward_reach[invalid_candidate]);
            }
            if (debug) {
                std::cerr << "candidate " << invalid_candidate << " (node " << order[invalid_candidate] <<  ") fails forward reach test against " << forward_reach << '\n';
            }
        }
        if (debug){
            if (!candidate_stack.empty()) {
                std::cerr << "comparing position " << i << " to backward reach of " << candidate_stack.back() << ": " << backward_reach[candidate_stack.back()] << '\n';
            }
        }
        
        if (!candidate_stack.empty() && backward_reach[candidate_stack.back()] == i) {
            // we found a superbubble interval
            auto confirmed_candidate = candidate_stack.back();
            superbubbles.emplace_back(order[i], order[confirmed_candidate]);
            candidate_stack.pop_back();
            if (!candidate_stack.empty()) {
                // give backward reach info to candidate parent
                backward_reach[candidate_stack.back()] = std::min(backward_reach[candidate_stack.back()],
                                                                  backward_reach[confirmed_candidate]);
            }
            if (debug) {
                std::cerr << "candidate " << confirmed_candidate << " (node " << order[confirmed_candidate] <<  ") is confirmed as an exit\n";
            }
        }
        for (auto node_id : graph.previous(order[i])) {
            backward_reach[i] = std::min<int64_t>(backward_reach[i], index[node_id]);
        }
        if (!candidate_stack.empty()) {
            backward_reach[candidate_stack.back()] = std::min(backward_reach[candidate_stack.back()],
                                                              backward_reach[i]);
        }
    }
    
    return superbubbles;
}

inline size_t SuperbubbleTree::chain_size() const {
    return chains.size();
}

inline size_t SuperbubbleTree::superbubble_size() const {
    return superbubbles.size();
}

inline uint64_t SuperbubbleTree::superbubble_ending_at(uint64_t node_id) const {
    return superbubble_endings[node_id];
}

inline uint64_t SuperbubbleTree::superbubble_beginning_at(uint64_t node_id) const {
    return superbubble_beginnings[node_id];
}

inline const std::pair<uint64_t, uint64_t>& SuperbubbleTree::superbubble_boundaries(uint64_t superbubble_id) const {
    return superbubbles[superbubble_id].boundaries;
}

inline const std::vector<uint64_t>& SuperbubbleTree::chains_inside(uint64_t superbubble_id) const {
    return superbubbles[superbubble_id].chain_ids;
}

inline uint64_t SuperbubbleTree::chain_containing(uint64_t superbubble_id) const {
    return superbubbles[superbubble_id].parent;
}

inline const std::vector<uint64_t>& SuperbubbleTree::superbubbles_inside(uint64_t chain_id) const {
    return chains[chain_id].superbubble_ids;
}

inline uint64_t SuperbubbleTree::superbubble_containing(uint64_t chain_id) const {
    return chains[chain_id].parent;
}





template<class Graph>
NetGraph::NetGraph(const Graph& graph, const SuperbubbleTree& superbubbles,
                   uint64_t superbubble_id) {
        
    static const bool debug = false;
    
    
    uint64_t start, end;
    std::tie(start, end) = superbubbles.superbubble_boundaries(superbubble_id);
    
    if (debug) {
        std::cerr << "building netgraph for superbubble " << superbubble_id << " with boundaries " << start << ", " << end << '\n';
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
                auto next_bub_id = superbubbles.superbubble_beginning_at(next_id);
                if (next_bub_id != -1 && next_id != end) {
                    // new chain node
                    auto chain_id = superbubbles.chain_containing(next_bub_id);
                    auto chain_net_id = add_node(chain_id, true);
                    
                    // jump to the end of the chain
                    auto final_bub_id = superbubbles.superbubbles_inside(chain_id).back();
                    auto final_node_id = superbubbles.superbubble_boundaries(final_bub_id).second;
                    
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
        std::cerr << "final topology for bubble " << superbubble_id << " net graph" << "\n";
        print(std::cerr);
    }
}

template<class Graph>
NetGraph::NetGraph(const Graph& graph, const SuperbubbleTree& superbubbles) : NetGraph(graph, superbubbles, nullptr) {
    
}

template<class Graph>
NetGraph::NetGraph(const Graph& graph, const SuperbubbleTree& superbubbles,
                   const SentinelTableau& tableau) : NetGraph(graph, superbubbles, &tableau) {
    
}

template<class Graph>
NetGraph::NetGraph(const Graph& graph, const SuperbubbleTree& superbubbles,
                   const SentinelTableau* tableau) {
    
    static const bool debug = false;
    
    // mark all the nodes that are contained in a superbubble
    std::vector<bool> contained(graph.node_size());
    for (uint64_t bub_id = 0; bub_id < superbubbles.superbubble_size(); ++bub_id) {
        NetGraph net_graph(graph, superbubbles, bub_id);
        for (uint64_t net_id = 0; net_id < net_graph.node_size(); ++net_id) {
            uint64_t feature_id;
            bool is_chain;
            std::tie(feature_id, is_chain) = net_graph.label(net_id);
            if (!is_chain) {
                contained[feature_id] = true;
                if (debug) {
                    std::cerr << feature_id << " marked contained from bubble " << bub_id << '\n';
                }
            }
        }
    }
    
    if (debug) {
        std::cerr << "making nodes\n";
    }
    
    // add chain nodes
    std::unordered_map<std::pair<uint64_t, bool>, uint64_t> forward_translation;
    for (uint64_t chain_id = 0; chain_id < superbubbles.chain_size(); ++chain_id) {
        if (superbubbles.superbubble_containing(chain_id) == -1) {
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
            feature_id = superbubbles.superbubble_boundaries(superbubbles.superbubbles_inside(feature_id).back()).second;
        }
        for (auto next_id : graph.next(feature_id)) {
            if (tableau && next_id == tableau->snk_id) {
                continue;
            }
            auto bub_id = superbubbles.superbubble_beginning_at(next_id);
            uint64_t next_net_id;
            if (bub_id == -1) {
                next_net_id = forward_translation.at(std::make_pair(next_id, false));
            }
            else {
                auto chain_id = superbubbles.chain_containing(bub_id);
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

#endif /* centrolign_superbubbles_hpp */
