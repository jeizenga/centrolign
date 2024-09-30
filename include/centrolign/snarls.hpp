#ifndef centrolign_snarls_hpp
#define centrolign_snarls_hpp

#include "centrolign/structure_tree.hpp"
#include "centrolign/is_acyclic.hpp"
#include "centrolign/cactus.hpp"

namespace centrolign {

/*
 * Data structure that finds snarls and provides navigation
 * methods for the snarl/chain tree
 */
class SnarlTree : public TwoDisconnectedStructureTree {
public:
    
    // construct from a graph with sentinel sink and source nodes, which
    // will be removed from the snarl decomposition
    template<class Graph>
    SnarlTree(const Graph& graph, const SentinelTableau& tableau);
    
    SnarlTree() = default;
    SnarlTree(const SnarlTree& other) = default;
    SnarlTree(SnarlTree&& other) = default;
    ~SnarlTree() = default;
    SnarlTree& operator=(const SnarlTree& other) = default;
    SnarlTree& operator=(SnarlTree&& other) = default;
    
    // does the chain contain any cycles
    inline bool chain_is_acyclic(uint64_t chain_id) const;
    // does the snarl contain any cycles
    inline bool snarl_is_acyclic(uint64_t snarl_id) const;
    // does the snarl's net graph contain any cycles
    inline bool net_graph_is_acyclic(uint64_t snarl_id) const;
    
protected:
    
    // parent class calls this to find snarls using the cactus tree
    template<class Graph>
    static std::vector<std::pair<uint64_t, uint64_t>> find_2_disc_structures_impl(const Graph& graph,
                                                                                  const SentinelTableau* tableau);
    
    std::vector<bool> chain_acyclic;
    std::vector<bool> snarl_acyclic;
    std::vector<bool> net_graph_acyclic;
    
    friend class TwoDisconnectedStructureTree;
};





/*
 * Template implementations
 */

template<class Graph>
SnarlTree::SnarlTree(const Graph& graph, const SentinelTableau& tableau) {
    
    // form the tree and chains
    TwoDisconnectedStructureTree::initialize<SnarlTree, Graph>(graph, &tableau);
    
    static const bool debug = false;
    
    
    // identify which snarls, net graphs, and chains are acyclic
    chain_acyclic.resize(chain_size());
    snarl_acyclic.resize(structure_size());
    net_graph_acyclic.resize(structure_size());
    
    for (auto feature : postorder()) {
        uint64_t feature_id;
        bool is_chain;
        std::tie(feature_id, is_chain) = feature;
        
        if (debug) {
            std::cerr << "check acyclicity at feature " << feature_id << ' ' << is_chain << '\n';
        }
        
        if (is_chain) {
            
            bool acyclic = true;
            for (auto snarl_id : structures_inside(feature_id)) {
                if (!snarl_acyclic[snarl_id]) {
                    acyclic = false;
                    break;
                }
            }
            
            chain_acyclic[feature_id] = acyclic;
        }
        else {
            {
                if (debug) {
                    std::cerr << "testing net graph acyclicity\n";
                }
                NetGraph net_graph(graph, *this, feature_id);
                net_graph_acyclic[feature_id] = is_acyclic(net_graph);
            }
            if (net_graph_acyclic[feature_id]) {
                bool acyclic = true;
                for (auto chain_id : chains_inside(feature_id)) {
                    if (!chain_acyclic[chain_id]) {
                        acyclic = false;
                        break;
                    }
                }
                snarl_acyclic[feature_id] = acyclic;
            }
            else {
                snarl_acyclic[feature_id] = false;
            }
        }
    }
    
    if (debug) {
        std::cerr << "final ayclicity:\n";
        std::cerr << "snarl\n";
        for (auto b : snarl_acyclic) {
            std::cerr << '\t' << b << '\n';
        }
        std::cerr << "net graph\n";
        for (auto b : net_graph_acyclic) {
            std::cerr << '\t' << b << '\n';
        }
        std::cerr << "chain\n";
        for (auto b : chain_acyclic) {
            std::cerr << '\t' << b << '\n';
        }
    }
}

template<class Graph>
std::vector<std::pair<uint64_t, uint64_t>> SnarlTree::find_2_disc_structures_impl(const Graph& graph,
                                                                                  const SentinelTableau* tableau) {
    
    assert(tableau != nullptr);

    std::vector<std::pair<uint64_t, uint64_t>> snarls;
    
    CactusGraph<Graph> cactus_graph(graph, *tableau);
    CactusTree cactus_tree(cactus_graph);
    
    // get the walk of the edge in the forward orientation and add its trivial snarls to the output
    auto get_edge_walk = [&](const std::tuple<uint64_t, bool, size_t>& edge) -> std::vector<uint64_t> {
        std::vector<uint64_t> walk;
        if (std::get<1>(edge)) {
            walk = std::move(cactus_graph.next_edge_label(std::get<0>(edge), std::get<2>(edge)));
        }
        else {
            walk = std::move(cactus_graph.previous_edge_label(std::get<0>(edge), std::get<2>(edge)));
        }
        for (size_t i = 1; i < walk.size(); ++i) {
            snarls.emplace_back(walk[i - 1], walk[i]);
        }
        return walk;
    };
    
    // traverse the tree from root down
    std::vector<uint64_t> stack(1, cactus_tree.get_root());
    while (!stack.empty()) {
        auto node_id = stack.back();
        stack.pop_back();
        
        if (cactus_tree.is_chain_node(node_id)) {
            
            // add the chain pairs
            const auto& chain = cactus_tree.chain(node_id);
            auto prev_walk = get_edge_walk(chain.front());
            for (size_t i = 1; i < chain.size(); ++i) {
                auto walk = get_edge_walk(chain[i]);
                if (std::get<1>(chain[i - 1]) == std::get<1>(chain[i])) {
                    // the snarl is consistently oriented
                    if (std::get<1>(chain[i])) {
                        snarls.emplace_back(prev_walk.back(), walk.front());
                    }
                    else {
                        snarls.emplace_back(walk.back(), prev_walk.front());
                    }
                }
                
                prev_walk = std::move(walk);
            }
        }
        
        for (auto next_id : cactus_tree.get_children(node_id)) {
            stack.push_back(next_id);
        }
    }
    
    return snarls;
}


inline bool SnarlTree::chain_is_acyclic(uint64_t chain_id) const {
    return chain_acyclic[chain_id];
}

inline bool SnarlTree::snarl_is_acyclic(uint64_t snarl_id) const {
    return snarl_acyclic[snarl_id];
}

inline bool SnarlTree::net_graph_is_acyclic(uint64_t snarl_id) const {
    return net_graph_acyclic[snarl_id];
}


}

#endif /* centrolign_snarls_hpp */
