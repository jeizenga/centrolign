#ifndef centrolign_utility_hpp
#define centrolign_utility_hpp

#include <vector>
#include <cstdint>
#include <iterator>
#include <string>
#include <iostream>
#include <utility>
#include <fstream>

#include "centrolign/graph.hpp"

namespace centrolign {

// return a vector of 0,1,...,size-1
std::vector<size_t> range_vector(size_t size);

// convert order into index or vice versa
template<typename Int>
std::vector<Int> invert(const std::vector<Int>& order);

// returns the highest 1-bit (or 0 in case of 0)
uint32_t hi_bit(uint64_t x);

// saturating int addition
uint64_t sat_add(uint64_t a, uint64_t b);
// saturating int multiplication
uint64_t sat_mult(uint64_t a, uint64_t b);

// convert ACGTN sequences into 01234 sequences
uint8_t encode_base(char nt);
std::string encode_seq(const std::string& seq);
// convert 01234 sequences into ACTGN
char decode_base(uint8_t nt);
std::string decode_seq(const std::string& seq);

// read from stdin or from a file
std::istream* get_input(const std::string& stream_name, std::ifstream& openable);

// return pairs of (name, sequence) from a FASTA
std::vector<std::pair<std::string, std::string>> parse_fasta(std::istream& in);

int64_t parse_int(const std::string& str);
double parse_double(const std::string& str);

std::string path_to_string(const BaseGraph& graph, const std::vector<uint64_t>& path);

template<class BGraph>
void print_graph(const BGraph& graph, std::ostream& out);

// adapter that uses reverse iteration instead of forward
template<class ReverseIterable>
class ReverseForEachAdapter {
public:
    ReverseForEachAdapter(const ReverseIterable& iterable) : iterable(iterable) {}
    ReverseForEachAdapter() = delete;
    ~ReverseForEachAdapter() = default;
    using const_iterator = typename ReverseIterable::const_reverse_iterator;
    const_iterator begin() const { return iterable.rbegin(); }
    const_iterator end() const { return iterable.rend(); }
private:
    const ReverseIterable& iterable;
};



/*
 * Template implementations
 */


template<typename Int>
std::vector<Int> invert(const std::vector<Int>& order) {
    std::vector<Int> inverted(order.size());
    for (Int i = 0; i < order.size(); ++i) {
        inverted[order[i]] = i;
    }
    return inverted;
}


template<class BGraph>
void print_graph(const BGraph& graph, std::ostream& out) {
    for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
        char base = graph.label(node_id);
        if (base <= 4) {
            base = decode_base(base);
        }
        out << node_id << ": " << base << '\n';
        for (auto next_id : graph.next(node_id)) {
            out << "\t-> " << next_id << '\n';
        }
    }
    if (graph.path_size() != 0) {
        out << "paths:\n";
        for (uint64_t path_id = 0; path_id < graph.path_size(); ++path_id) {
            bool first = true;
            for (auto node_id : graph.path(path_id)) {
                if (!first) {
                    out << ", ";
                }
                out << node_id;
                first = false;
            }
            out << '\n';
        }
    }
}

}

// http://stackoverflow.com/questions/4870437/pairint-int-pair-as-key-of-unordered-map-issue#comment5439557_4870467
// https://github.com/Revolutionary-Games/Thrive/blob/fd8ab943dd4ced59a8e7d1e4a7b725468b7c2557/src/util/pair_hash.h
// taken from boost
#ifndef OVERLOAD_PAIR_HASH
#define OVERLOAD_PAIR_HASH
namespace std {
namespace
{

// Code from boost
// Reciprocal of the golden ratio helps spread entropy
//     and handles duplicates.
// See Mike Seymour in magic-numbers-in-boosthash-combine:
//     http://stackoverflow.com/questions/4948780

template <class T>
inline void hash_combine(size_t& seed, T const& v)
{
    seed ^= std::hash<T>()(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
}

// Recursive template code derived from Matthieu M.
template <class Tuple, size_t Index = std::tuple_size<Tuple>::value - 1>
struct HashValueImpl
{
    static void apply(size_t& seed, Tuple const& tuple)
    {
        HashValueImpl<Tuple, Index-1>::apply(seed, tuple);
        hash_combine(seed, std::get<Index>(tuple));
    }
};

template <class Tuple>
struct HashValueImpl<Tuple,0>
{
    static void apply(size_t& seed, Tuple const& tuple)
    {
        hash_combine(seed, std::get<0>(tuple));
    }
};
}

template <typename A, typename B>
struct hash<pair<A,B> > {
    size_t operator()(const pair<A,B>& x) const {
        size_t hash_val = std::hash<A>()(x.first);
        hash_combine(hash_val, x.second);
        return hash_val;
    }
};

// from http://stackoverflow.com/questions/7110301/generic-hash-for-tuples-in-unordered-map-unordered-set
template <typename ... TT>
struct hash<std::tuple<TT...>>
{
    size_t
    operator()(std::tuple<TT...> const& tt) const
    {
        size_t seed = 0;
        HashValueImpl<std::tuple<TT...> >::apply(seed, tt);
        return seed;
    }
    
};

}
#endif  // OVERLOAD_PAIR_HASH

#endif /* centrolign_utility_hpp */
