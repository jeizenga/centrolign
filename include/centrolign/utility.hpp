#ifndef centrolign_utility_hpp
#define centrolign_utility_hpp

#include <vector>
#include <cstdint>
#include <iterator>
#include <string>
#include <iostream>
#include <utility>
#include <fstream>
#include <limits>
#include <cmath>
#include <stdexcept>
#include <functional>

#include "centrolign/graph.hpp"
#include "centrolign/logging.hpp"

namespace centrolign {

// return a vector of 0,1,...,size-1
std::vector<size_t> range_vector(size_t size);

// convert order into index or vice versa
template<typename Int>
std::vector<Int> invert(const std::vector<Int>& order);

// reorder a vector according to a permutation (converts permutation to identity)
template<typename Int, class T>
void reorder(std::vector<T>& vec, std::vector<Int>& to_index);

template<class T>
void filter_by_index(std::vector<T>& vec, const std::function<bool(size_t)>& remove);

// returns the highest 1-bit (or 0 in case of 0)
uint32_t hi_bit(uint64_t x);

// saturating int addition
inline uint64_t sat_add(uint64_t a, uint64_t b);
// saturating int multiplication
inline uint64_t sat_mult(uint64_t a, uint64_t b);
// overflow-resistant addition of log-tranformed values
inline double add_log(double log_x, double log_y);
// Holder's generalized means, implementation only valid for positive values
template<class Iterator>
inline double generalized_mean(Iterator begin, Iterator end, double p);

// median between random access iterators, reorders the range, returns midpoint for even-length ranges
template<class Iterator>
typename std::iterator_traits<Iterator>::value_type median(Iterator begin, Iterator end);

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
bool parse_bool(const std::string& str);

template<typename IntType>
std::string to_hex(const IntType& value);

// return the current resident set size in bytes
int64_t current_memory_usage();
// return the maximum resident set size of the program so far in bytes
int64_t max_memory_usage();
// format the memory usage as a reader-friendly string
std::string format_memory_usage(int64_t mem);
// log the current memory usage
void log_memory_usage(logging::LoggingLevel level);


std::vector<std::string> tokenize(const std::string& str, char delim = '\t');
std::string join(const std::vector<std::string>& values, const std::string& sep);

std::string path_to_string(const BaseGraph& graph, const std::vector<uint64_t>& path);

template<class BGraph>
void print_graph(const BGraph& graph, std::ostream& out);

template<class Graph>
void print_topology(const Graph& graph, std::ostream& out);

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

// add a vector-like interface to a slice of another vector-like object
template<class VectorLike>
class SliceView {
public:
    SliceView(const SliceView<VectorLike>& slice, size_t begin, size_t end) : vec(slice.vec), _begin(slice._begin + begin), _size(end - begin) {
        if (end < begin || end > slice.size()) {
            throw std::runtime_error("Invalid slice " + std::to_string(begin) + ":" + std::to_string(end) + " in object of size " + std::to_string(slice.size()));
        }
    }
    SliceView(const VectorLike& vec, size_t begin, size_t end) : vec(vec), _begin(begin), _size(end - begin) {
        if (end < begin || end > vec.size()) {
            throw std::runtime_error("Invalid slice " + std::to_string(begin) + ":" + std::to_string(end) + " in object of size " + std::to_string(vec.size()));
        }
    }
    SliceView() = delete;
    ~SliceView() = default;
    
    using value_type = typename VectorLike::value_type;
    
    typename VectorLike::const_iterator begin() const {
        return vec.begin() + _begin;
    }
    typename VectorLike::const_iterator end() const {
        return vec.begin() + _begin + _size;
    }
    const typename VectorLike::value_type& operator[](size_t i) const {
        return vec[_begin + i];
    }
    size_t size() const {
        return _size;
    }
    bool empty() const {
        return _size == 0;
    }
private:
    const VectorLike& vec;
    size_t _begin;
    size_t _size;
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


template<typename Int, class T>
void reorder(std::vector<T>& vec, std::vector<Int>& to_index) {
    for (size_t i = 0; i < vec.size(); ++i) {
        while (to_index[i] != i) {
            std::swap(vec[to_index[i]], vec[i]);
            std::swap(to_index[to_index[i]], to_index[i]);
        }
    }
}

template<class T>
void filter_by_index(std::vector<T>& vec, const std::function<bool(size_t)>& remove) {
    size_t removed = 0;
    for (size_t i = 0; i < vec.size(); ++i) {
        if (remove(i)) {
            ++removed;
        }
        else if (removed) {
            vec[i - removed] = std::move(vec[i]);
        }
    }
    if (removed) {
        vec.resize(vec.size() - removed);
    }
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

template<class Graph>
void print_topology(const Graph& graph, std::ostream& out) {
    for (uint64_t n = 0; n < graph.node_size(); ++n) {
        out << n << ":\n";
        for (auto m : graph.next(n)) {
            out << '\t' << m << '\n';
        }
    }
}

inline uint64_t sat_add(uint64_t a, uint64_t b) {
    if (std::numeric_limits<uint64_t>::max() - a >= b) {
        return a + b;
    }
    else {
        return std::numeric_limits<uint64_t>::max();
    }
}

inline uint64_t sat_mult(uint64_t a, uint64_t b) {
    if (b == 0 || a <= std::numeric_limits<uint64_t>::max() / b) {
        return a * b;
    }
    else {
        return std::numeric_limits<uint64_t>::max();
    }
}

inline double add_log(double log_x, double log_y) {
    if (log_x < log_y) {
        return log_y + log1p(exp(log_x - log_y));
    }
    else {
        return log_x + log1p(exp(log_y - log_x));
    }
}

template<class Iterator>
inline double generalized_mean(Iterator begin, Iterator end, double p) {
    double n = 0.0;
    if (p == 0.0) {
        // sum of logs
        double sum_log = 0.0;
        for (auto it = begin; it != end; ++it) {
            sum_log += log(*it);
            n += 1.0;
        }
        return exp(sum_log / n);
    }
    else {
        // log of sum
        double log_sum = std::numeric_limits<double>::lowest() / 2.0; // a small epsilon
        for (auto it = begin; it != end; ++it) {
            log_sum = add_log(log_sum, p * log(*it));
            n += 1.0;
        }
        return exp((log_sum - log(n)) / p);
    }
}


template<class Iterator>
typename std::iterator_traits<Iterator>::value_type median(Iterator begin, Iterator end) {
    auto mid = begin + (end - begin) / 2;
    std::nth_element(begin, mid, end);
    if ((end - begin) % 2 == 0) {
        return (*mid + *std::max_element(begin, mid)) / 2;
    }
    else {
        return *mid;
    }
}


template<typename IntType>
std::string to_hex(const IntType& value) {
    
    static const char hex_values[16] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F'};
    
    volatile int32_t test_val = 1;
    bool little_endian = (*((uint8_t*) &test_val) == 1);
    
    std::string hex_string(sizeof(IntType) * 2, '\0');
    
    const uint8_t* data = (const uint8_t*) &value;
    
    size_t j = little_endian ? sizeof(IntType) - 1 : 0;
    size_t incr = little_endian ? -1 : 1;
    for (int64_t i = 0; i < sizeof(IntType); ++i, j += incr) {
        hex_string[2 * i] = hex_values[(data[j] & 0xF0) >> 4];
        hex_string[2 * i + 1] = hex_values[data[j] & 0xF];
    }
    
    return hex_string;
}

}

// http://stackoverflow.com/questions/4870437/pairint-int-pair-as-key-of-unordered-map-issue#comment5439557_4870467
// https://github.com/Revolutionary-Games/Thrive/blob/fd8ab943dd4ced59a8e7d1e4a7b725468b7c2557/src/util/pair_hash.h
// taken from boost
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

#endif /* centrolign_utility_hpp */
