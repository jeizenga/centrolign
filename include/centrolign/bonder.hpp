#ifndef centrolign_bonder_hpp
#define centrolign_bonder_hpp

#include <vector>
#include <array>
#include <utility>

#include "centrolign/anchorer.hpp"
#include "centrolign/partitioner.hpp"
#include "centrolign/utility.hpp"

namespace centrolign {

struct bond_t;

struct bond_t {
    
    bond_t() noexcept = default;
    bond_t(const bond_t& other) noexcept = default;
    bond_t(bond_t&& other) noexcept = default;
    ~bond_t() = default;
    bond_t& operator=(const bond_t& other) noexcept = default;
    bond_t& operator=(bond_t&& other) noexcept = default;
    
    std::string path1;
    std::string path2;
    size_t offset1;
    size_t offset2;
    size_t length;
};

/*
 * Class to identify cycles for the alignment from anchors
 */
class Bonder : public Extractor, public PartitionClient {
public:
    Bonder() = default;
    ~Bonder() = default;
    
    template<class BGraph, class XMerge>
    std::vector<std::vector<bond_t>> identify_bonds(const BGraph& graph1, const BGraph& graph2,
                                                    const SentinelTableau& tableau1, const SentinelTableau& tableau2,
                                                    const XMerge& xmerge1, const XMerge& xmerge2,
                                                    const std::vector<anchor_t>& opt_chain,
                                                    const std::vector<anchor_t>& secondary_chain) const;
    
    double min_opt_proportion = 0.9;
    
    double min_length = 10000.0;
    
protected:
    
    // passed in as records of (length, opt segment score, secondary segment score)
    std::vector<std::pair<size_t, size_t>> longest_partition(const std::vector<std::tuple<double, double, double>>& shared_subanchors,
                                                             const std::vector<std::tuple<double, double, double>>& intervening_segments) const;
    
};



/*
 * Template implementations
 */

template<class BGraph, class XMerge>
std::vector<std::vector<bond_t>> Bonder::identify_bonds(const BGraph& graph1, const BGraph& graph2,
                                                        const SentinelTableau& tableau1, const SentinelTableau& tableau2,
                                                        const XMerge& xmerge1, const XMerge& xmerge2,
                                                        const std::vector<anchor_t>& opt_chain,
                                                        const std::vector<anchor_t>& secondary_chain) const {
    
    
    std::vector<std::vector<bond_t>> bonds;
    
    std::vector<std::pair<size_t, size_t>> node_locations(graph1.node_size(),
                                                          std::pair<size_t, size_t>(-1, -1));
    
    for (size_t i = 0; i < opt_chain.size(); ++i) {
        auto& anchor = opt_chain[i];
        for (size_t j = 0; j < anchor.walk1.size(); ++j) {
            node_locations[anchor.walk1[j]] = std::make_pair(i, j);
        }
    }
    
//    // pairs (one from each chain) of anchored node pairs (one from each graph) with their score
//    std::vector<std::array<std::pair<std::vector<std::pair<uint64_t, uint64_t>>, double>, 2>> shared_subanchors;
//    std::vector<std::array<std::pair<std::vector<std::pair<uint64_t, uint64_t>>, double>, 2>> parallel_subanchors;

    
    // records of (opt anchor, idx on opt anchor, sec anchor, idx on sec anchor, length)
    std::vector<std::tuple<size_t, size_t, size_t, size_t, size_t>> shared_subanchors;
    
    size_t prev_k = -1, prev_l = -1;
    for (size_t i = 0; i < secondary_chain.size(); ++i) {
        const auto& anchor = secondary_chain[i];
        for (size_t j = 0; j < anchor.walk1.size(); ++j) {
            size_t k, l;
            std::tie(k, l) = node_locations[anchor.walk1[j]];
            if (k != -1) {
                if (prev_k == k && prev_l == l - 1) {
                    ++std::get<4>(shared_subanchors.back());
                }
                else {
                    shared_subanchors.emplace_back(i, j, k, l, 1);
                }
            }
        }
    }
    
    return bonds;
}


}

#endif /* centrolign_bond_identifier_hpp */
