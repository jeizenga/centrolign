#ifndef centrolign_stitcher_hpp
#define centrolign_stitcher_hpp

#include <vector>
#include <cstdint>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

#include "centrolign/chain_merge.hpp"
#include "centrolign/alignment.hpp"
#include "centrolign/anchorer.hpp"
#include "centrolign/subgraph_extraction.hpp"
#include "centrolign/logging.hpp"

namespace centrolign {

/*
 * Object that connects an anchor chain into a base-level alignment
 */
class Stitcher {
public:
    
    Stitcher();
    ~Stitcher() = default;
    
    template<class BGraph1, class BGraph2, class XMerge1, class XMerge2>
    Alignment stitch(const std::vector<anchor_t>& anchor_chain,
                     const BGraph1& graph1, const BGraph2& graph2,
                     const SentinelTableau& tableau1, const SentinelTableau& tableau2,
                     const XMerge1& chain_merge1, const XMerge2& chain_merge2) const;
    
    // TODO: think out the defaults better
    // defaults: m = 4, x = 8, o1 = 6, e1 = 4, o2 = 150, e1 = 1
    AlignmentParameters<2> alignment_params;
    
private:
    
    static const bool debug;
    static const bool instrument;
    
    void subalign(const SubGraphInfo& extraction1, const SubGraphInfo& extraction2,
                  Alignment& stitched) const;
    
    void do_instrument(const SubGraphInfo& extraction1, const SubGraphInfo& extraction2,
                       int64_t score) const;
    
};

/*
 * Template implementations
 */

template<class BGraph1, class BGraph2, class XMerge1, class XMerge2>
Alignment Stitcher::stitch(const std::vector<anchor_t>& anchor_chain,
                           const BGraph1& graph1, const BGraph2& graph2,
                           const SentinelTableau& tableau1, const SentinelTableau& tableau2,
                           const XMerge1& chain_merge1, const XMerge2& chain_merge2) const {
    
    if (anchor_chain.empty()) {
        throw std::runtime_error("Stitcher cannot stitch an empty anchor chain");
    }
    
    size_t next_log_idx = 0;
    std::vector<size_t> logging_indexes;
    if (logging::level >= logging::Debug) {
        for (size_t i = 1; i < 10; ++i) {
            logging_indexes.push_back((anchor_chain.size() * i) / 10);
        }
        auto end = std::unique(logging_indexes.begin(), logging_indexes.end());
        logging_indexes.resize(end - logging_indexes.begin());
        
        logging::log(logging::Debug, "Stitching a chain of " + std::to_string(anchor_chain.size()) + " anchors");
    }
    
    Alignment stitched;
    
    // left end alignment
    {
        auto extraction1 = extract_connecting_graph(graph1, tableau1.src_id,
                                                    anchor_chain.front().walk1.front(),
                                                    chain_merge1);
        auto extraction2 = extract_connecting_graph(graph2, tableau2.src_id,
                                                    anchor_chain.front().walk2.front(),
                                                    chain_merge2);
                
        subalign(extraction1, extraction2, stitched);
    }
    
    for (size_t i = 0; i < anchor_chain.size(); ++i) {
        if (next_log_idx < logging_indexes.size() && i == logging_indexes[next_log_idx]) {
            logging::log(logging::Debug, "Alignment stitching iteration " + std::to_string(i + 1) + " of " + std::to_string(anchor_chain.size()));
            ++next_log_idx;
        }
        
        const auto& anchor = anchor_chain[i];
        if (i != 0) {
            // make intervening alignment
            const auto& prev_anchor = anchor_chain[i - 1];
            
            auto extraction1 = extract_connecting_graph(graph1,
                                                        prev_anchor.walk1.back(),
                                                        anchor.walk1.front(),
                                                        chain_merge1);
            auto extraction2 = extract_connecting_graph(graph2,
                                                        prev_anchor.walk2.back(),
                                                        anchor.walk2.front(),
                                                        chain_merge2);
            
            subalign(extraction1, extraction2, stitched);
        }
        // copy the anchor
        for (size_t j = 0; j < anchor.walk1.size(); ++j) {
            stitched.emplace_back(anchor.walk1[j], anchor.walk2[j]);
        }
        
        if (debug) {
            std::cerr << "anchor " << i << ":\n";
            for (size_t j = 0; j < anchor.walk1.size(); ++j) {
                std::cerr << '\t' << anchor.walk1[j] << '\t' << anchor.walk2[j] << '\n';
            }
        }
    }
    
    // right end alignment
    {
        auto extraction1 = extract_connecting_graph(graph1, anchor_chain.back().walk1.back(),
                                                    tableau1.snk_id, chain_merge1);
        auto extraction2 = extract_connecting_graph(graph2, anchor_chain.back().walk2.back(),
                                                    tableau2.snk_id, chain_merge2);
        
        subalign(extraction1, extraction2, stitched);
    }
    
    return stitched;
}

}

#endif /* centrolign_stitcher_hpp */
