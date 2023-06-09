#ifndef centrolign_simplifier_hpp
#define centrolign_simplifier_hpp

#include <deque>
#include <limits>
#include <iostream>
#include <sstream>

#include "centrolign/graph.hpp"
#include "centrolign/modify_graph.hpp"
#include "centrolign/trie.hpp"
#include "centrolign/step_index.hpp"
#include "centrolign/superbubbles.hpp"

namespace centrolign {

/*
 * Struct to summarize results of expanding a graph
 */
struct ExpandedGraph {
    // the expanded graph
    BaseGraph graph;
    // the translation of expanded graph node IDs to original node IDs
    std::vector<uint64_t> back_translation;
    // a revised sentinel tableau
    SentinelTableau tableau;
    
    ExpandedGraph() = default;
    ~ExpandedGraph() = default;
};

/*
 * Object that simplifies complex regions of the graph
 */
class Simplifier {
public:
    
    Simplifier() = default;
    ~Simplifier() = default;
    
    ExpandedGraph simplify(const BaseGraph& graph, const SentinelTableau& tableau) const;
    
    // simplify the regions in front of the target nodes, up to a max distance
    ExpandedGraph targeted_simplify(const BaseGraph& graph, const SentinelTableau& tableau,
                                    const std::vector<uint64_t>& node_ids, size_t distance) const;
    
    std::vector<std::vector<uint64_t>> identify_target_nodes(const std::vector<std::vector<uint64_t>>& node_counts) const;
    
    /*
     * Configurable parameters
     */
    
    // consider windows of chains with at least this min distance across
    size_t min_dist_window = 128;
    // exclude superbubble from pruning if they contain a path of this size
    size_t preserve_bubble_size = 32;
    // largest number of walks allowed before separating paths
    uint64_t max_walks = 24;
    // the minimum proportion of nodes we will include as targets for resimplification
    double min_resimplify_fraction = 0.01;
    // target all nodes above this count during resimplification, even if above the
    // minimum proportion
    uint64_t max_resimplify_count = 1000;
    
private:
    
    static bool debug;
    
    std::vector<std::vector<uint64_t>> mergeable_nodes(const Trie& trie) const;
    
    ExpandedGraph perform_simplification(const BaseGraph& graph,
                                         const SentinelTableau& tableau,
                                         const StepIndex& step_index,
                                         std::vector<std::pair<Trie, uint64_t>>& interval_rev_tries,
                                         const std::vector<size_t>& node_to_trie) const;
    
    void simplify_chain_interval(const BaseGraph& graph, const StepIndex& step_index,
                                 const SuperbubbleTree& superbubbles,
                                 std::vector<std::pair<Trie, uint64_t>>& interval_rev_tries,
                                 std::vector<size_t>& node_to_trie,
                                 uint64_t chain_id, size_t begin, size_t end) const;
    
    // a specialized (non-general) multiprecision int
    // assumes that factors are the same for * and / and that they come
    // in the same relative order, with * happening first for each factor
    struct window_prod_t {
    public:
        // default value is 1
        window_prod_t() : factors(1, 1) { }
        window_prod_t(uint64_t val) : factors(1, val) { }
        window_prod_t(const window_prod_t& prod) : factors(prod.factors) { }
        
        inline window_prod_t operator*(uint64_t val) const;
        inline window_prod_t& operator*=(uint64_t val);
        inline window_prod_t operator/(uint64_t val) const;
        inline window_prod_t& operator/=(uint64_t val);
        
        inline bool operator<(uint64_t val) const;
        inline bool operator<=(uint64_t val) const;
        inline bool operator>(uint64_t val) const;
        inline bool operator>=(uint64_t val) const;
        
        inline window_prod_t& operator=(uint64_t val);
                
        inline std::string str() const;
        
    private:
        std::deque<uint64_t> factors;
    };
};




/*
 * Template and inline implementations
 */

inline Simplifier::window_prod_t Simplifier::window_prod_t::operator*(uint64_t val) const {
    window_prod_t result(*this);
    result *= val;
    return result;
}

inline Simplifier::window_prod_t& Simplifier::window_prod_t::operator*=(uint64_t val) {
    if (factors.back() <= std::numeric_limits<uint64_t>::max() / val) {
        // we won't overflow
        factors.back() *= val;
    }
    else {
        // add a new factor to avoid overflow
        factors.push_back(val);
    }
    return *this;
}

inline Simplifier::window_prod_t Simplifier::window_prod_t::operator/(uint64_t val) const {
    window_prod_t result(*this);
    result /= val;
    return result;
}

inline Simplifier::window_prod_t& Simplifier::window_prod_t::operator/=(uint64_t val) {
    factors.front() /= val;
    if (factors.size() > 1 && factors[1] <= std::numeric_limits<uint64_t>::max() / factors.front()) {
        // we can combine the first two factors without overflow
        factors[1] *= factors.front();
        factors.pop_front();
    }
    return *this;
}

inline bool Simplifier::window_prod_t::operator<(uint64_t val) const {
    return factors.size() == 1 && factors.front() < val;
}

inline bool Simplifier::window_prod_t::operator<=(uint64_t val) const {
    return !(*this > val);
}

inline bool Simplifier::window_prod_t::operator>(uint64_t val) const {
    // there is an invariant that the factors only splill over into 2
    // if they are above the uint64_t max
    return factors.size() > 1 || factors.front() > val;
}

inline bool Simplifier::window_prod_t::operator>=(uint64_t val) const {
    return !(*this < val);
}

inline Simplifier::window_prod_t& Simplifier::window_prod_t::operator=(uint64_t val) {
    factors.clear();
    factors.push_back(val);
    return *this;
}

std::string Simplifier::window_prod_t::str() const {
    std::stringstream sstrm;
    if (factors.size() == 1) {
        sstrm << factors.front();
    }
    else {
        sstrm << '(';
        for (size_t i = 0; i < factors.size(); ++i) {
            if (i) {
                sstrm << '*';
            }
            sstrm << factors[i];
        }
        sstrm << ')';
    }
    return sstrm.str();
}

}

#endif /* centrolign_simplifier_hpp */
