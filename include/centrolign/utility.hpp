#ifndef centrolign_utility_hpp
#define centrolign_utility_hpp

#include <vector>
#include <cstdint>
#include <iterator>
#include <string>
#include <iostream>

namespace centrolign {

// return a vector of 0,1,...,size-1
std::vector<size_t> range_vector(size_t size);

// returns the highest 1-bit (or 0 in case of 0)
uint32_t hi_bit(uint64_t x);

// convert ACGTN sequences into 01234 sequences
uint8_t encode_base(char nt);
std::string encode_seq(const std::string& seq);
// convert 01234 sequences into ACTGN
char decode_base(uint8_t nt);
std::string decode_seq(const std::string& seq);

// return pairs of (name, sequence) from a FASTA
std::vector<std::pair<std::string, std::string>> parse_fasta(std::istream& in);

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


template<class BGraph>
void print_graph(const BGraph& graph, std::ostream& out) {
    for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
        out << node_id << ": " << graph.base(node_id) << '\n';
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

#endif /* centrolign_utility_hpp */
