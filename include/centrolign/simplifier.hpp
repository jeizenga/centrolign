#ifndef centrolign_simplifier_hpp
#define centrolign_simplifier_hpp

#include <deque>
#include <limits>

#include "centrolign/graph.hpp"
#include "centrolign/modify_graph.hpp"

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
 * Object that
 */
class Simplifier {
public:
    
    Simplifier() = default;
    ~Simplifier() = default;
    
    ExpandedGraph simplify(const BaseGraph& graph, const SentinelTableau& tableau) const;
    
    /*
     * Configurable parameters
     */
    
    // consider windows of chains with at least this min distance across
    size_t min_dist_window = 128;
    // exclude superbubble from pruning if they contain a path of this size
    size_t preserve_bubble_size = 32;
    // largest number of walks allowed before separating paths
    uint64_t max_walks = 24;
    
private:
    
    // a specialized (non-general) multiprecision int
    // assumes that factors are the same for * and / and that they come
    // in the same relative order, with * happening first for each factor
    struct long_product_t {
    public:
        // default value is 1
        long_product_t() : factors(1, 1) { }
        long_product_t(uint64_t val) : factors(1, val) { }
        long_product_t(const long_product_t& prod) : factors(prod.factors) { }
        
        inline long_product_t operator*(uint64_t val) const;
        inline long_product_t& operator*=(uint64_t val);
        inline long_product_t operator/(uint64_t val) const;
        inline long_product_t& operator/=(uint64_t val);
        
        inline bool operator<(uint64_t val) const;
        inline bool operator<=(uint64_t val) const;
        inline bool operator>(uint64_t val) const;
        inline bool operator>=(uint64_t val) const;
        
        inline long_product_t& operator=(uint64_t val);
        
    private:
        std::deque<uint64_t> factors;
    };
};

/*
 * Template and inline implementations
 */


inline Simplifier::long_product_t Simplifier::long_product_t::operator*(uint64_t val) const {
    long_product_t result(*this);
    result *= val;
    return result;
}

inline Simplifier::long_product_t& Simplifier::long_product_t::operator*=(uint64_t val) {
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

inline Simplifier::long_product_t Simplifier::long_product_t::operator/(uint64_t val) const {
    long_product_t result(*this);
    result /= val;
    return result;
}

inline Simplifier::long_product_t& Simplifier::long_product_t::operator/=(uint64_t val) {
    factors.front() /= val;
    if (factors.size() > 1 && factors[1] <= std::numeric_limits<uint64_t>::max() / factors.front()) {
        // we can combine the first two factors without overflow
        factors[1] *= factors.front();
        factors.pop_front();
    }
    return *this;
}

inline bool Simplifier::long_product_t::operator<(uint64_t val) const {
    return factors.size() == 1 && factors.front() < val;
}

inline bool Simplifier::long_product_t::operator<=(uint64_t val) const {
    return !(*this > val);
}

inline bool Simplifier::long_product_t::operator>(uint64_t val) const {
    // there is an invariant that the factors only splill over into 2
    // if they are above the uint64_t max
    return factors.size() > 1 || factors.front() > val;
}

inline bool Simplifier::long_product_t::operator>=(uint64_t val) const {
    return !(*this < val);
}

inline Simplifier::long_product_t& Simplifier::long_product_t::operator=(uint64_t val) {
    factors.clear();
    factors.push_back(val);
    return *this;
}

}

#endif /* centrolign_simplifier_hpp */
