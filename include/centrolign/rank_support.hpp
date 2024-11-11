#ifndef centrolign_rank_support_hpp
#define centrolign_rank_support_hpp

#include <vector>
#include <cmath>
#include <algorithm>
#include <bitset>

#include "centrolign/packed_vector.hpp"

namespace centrolign {

/*
 * Jacobson's rank support for a bit vector
 */
class RankSupport {

public:
    
    template<class BitVector>
    RankSupport(const BitVector& bit_vector) noexcept;
    RankSupport() noexcept = default;
    RankSupport(const RankSupport& other) noexcept = default;
    RankSupport(RankSupport&& other) noexcept = default;
    ~RankSupport() noexcept = default;
    
    RankSupport& operator=(const RankSupport& other) noexcept = default;
    RankSupport& operator=(RankSupport&& other) noexcept = default;
    
    // the number of 1 bits in the bit vector in positions [0, ..., i - 1]
    inline size_t rank(size_t i) const;
    
private:
    
    // log^2 n
    size_t chunk_len = 0;
    // 1/2 * log n
    size_t subchunk_len = 0;
    size_t subchunks_per_chunk = 0; // this is not strictly necessary, but convenient to not have to recompute
    
    PackedVector chunk_rank;
    PackedVector subchunk_subrank;
    PackedVector subchunk_code;
    std::vector<PackedArray> subchunk_table;
};

/*
 * Template implementations
 */

template<class BitVector>
RankSupport::RankSupport(const BitVector& bit_vector) noexcept {
    
    static const bool debug = false;
    
    if (bit_vector.size() == 0) {
        return;
    }
    
    double logn = log2(bit_vector.size());
    chunk_len = std::max<size_t>(ceil(logn * logn), 1);
    subchunk_len = std::max<size_t>(ceil(0.5 * logn), 1);
    subchunks_per_chunk = (chunk_len + subchunk_len - 1) / subchunk_len;
    
    size_t num_chunks = (bit_vector.size() + chunk_len - 1) / chunk_len;
    size_t num_subchunks = num_chunks * subchunks_per_chunk;
    
    if (debug) {
        std::cerr << "bit vector length " << bit_vector.size() << '\n';
        std::cerr << "chunk len " << chunk_len << '\n';
        std::cerr << "subchunk len " << subchunk_len << '\n';
        std::cerr << "subchunks per chunk " << subchunks_per_chunk << '\n';
        std::cerr << "num chunks " << num_chunks << '\n';
        std::cerr << "num subchunk " << num_subchunks << '\n';
        
    }
    
    chunk_rank = std::move(PackedVector(num_chunks));
    subchunk_subrank = std::move(PackedVector(num_subchunks));
    subchunk_code = std::move(PackedVector(num_subchunks));
    
    // exhaustively enumerate the rank queries for all possible subchunks
    size_t num_codes = (1 << subchunk_len);
    subchunk_table.reserve(num_codes);
    for (size_t i = 0; i < num_codes; ++i) {
        subchunk_table.emplace_back(subchunk_len);
        auto& table_entry = subchunk_table.back();
        uint64_t rank = 0;
        for (size_t j = 0; j < subchunk_len; ++j) {
            table_entry.set(j, rank, subchunk_len);
            if (i & (1 << j)) {
                rank++;
            }
        }
    }
    
    size_t rank = 0;
    for (size_t i = 0, c = 0, s = 0; i < bit_vector.size(); i += chunk_len, ++c) {
        // iterate over chunks
        chunk_rank[c] = rank;
        size_t subrank = 0;
        for (size_t j = i, n = std::min<size_t>(bit_vector.size(), i + chunk_len); j < n; j += subchunk_len, ++s) {
            // iterate over sub-chunks
            subchunk_subrank[s] = subrank;
            uint64_t code = 0;
            for (size_t k = j, l = 0, m = std::min<size_t>(n, j + subchunk_len); k < m; ++k, ++l) {
                // iterate over indexes
                if (bit_vector[k]) {
                    code |= (1 << l);
                    ++subrank;
                }
            }
            subchunk_code[s] = code;
        }
        rank += subrank;
    }
    if (debug) {
        std::cerr << "exhaustive table:\n";
        for (size_t i = 0; i < subchunk_table.size(); ++i) {
            std::cerr << i << '\t' << std::bitset<5>(i);
            for (size_t j = 0; j < subchunk_len; ++j) {
                std::cerr << '\t' << subchunk_table[i].at(j);
            }
            std::cerr << '\n';
        }
        std::cerr << "chunk vectors\n";
        for (size_t i = 0, j = 0; i < chunk_rank.size(); ++i) {
            std::cerr << "--- chunk " << i << " rank " << chunk_rank.at(i) << " ---\n";
            for (size_t k = 0; k < subchunks_per_chunk; ++k, ++j) {
                std::cerr << "idx " << j << '\t' << "subidx " << k << '\t' << "subrank " << subchunk_subrank.at(j) << '\t' << "code " << subchunk_code.at(j) << '\t' << "set " << std::bitset<5>(subchunk_code.at(j)) << '\n';
            }
        }
    }
}

inline size_t RankSupport::rank(size_t i) const {
    size_t chunk = i / chunk_len;
    size_t idx_in_chunk = i % chunk_len;
    size_t subchunk_idx = idx_in_chunk / subchunk_len;
    size_t idx_in_subchunk = idx_in_chunk % subchunk_len;
    size_t subchunk = chunk * subchunks_per_chunk + subchunk_idx;
    //std::cerr << "i " << i << " idx in chunk " << idx_in_chunk << " subchunk idx " << subchunk_idx << " idx in subchunk " << idx_in_subchunk << " subchunk " << subchunk << '\n';
    return chunk_rank.at(chunk) + subchunk_subrank.at(subchunk) + subchunk_table[subchunk_code.at(subchunk)].at(idx_in_subchunk);
}

}

#endif /* centrolign_rank_support_hpp */
