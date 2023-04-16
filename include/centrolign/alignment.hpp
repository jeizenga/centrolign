#ifndef centrolign_alignment_hpp
#define centrolign_alignment_hpp

#include <vector>
#include <cstdint>
#include <array>
#include <limits>
#include <algorithm>

#include "centrolign/topological_order.hpp"
#include "centrolign/utility.hpp"

namespace centrolign {

/*
 * Two nodes aligned to each other, one of which might be a gap
 */
struct AlignedPair {
    AlignedPair() = default;
    ~AlignedPair() = default;
    AlignedPair(uint64_t node_id1, uint64_t node_id2);
    
    static const uint64_t gap;
    
    uint64_t node_id1 = gap;
    uint64_t node_id2 = gap;
    
    bool operator==(const AlignedPair& other) const;
};

typedef std::vector<AlignedPair> Alignment;

// affine gap score parameters
template<int NumPW = 1>
struct AlignmentParameters {
    uint32_t match;
    uint32_t mismatch;
    std::array<uint32_t, NumPW> gap_open;
    std::array<uint32_t, NumPW> gap_extend;
    
    AlignmentParameters() = default;
    ~AlignmentParameters() = default;
};

// PO-POA, can crash if either graph cannot reach one of the sinks from
// one of the sources
template<int NumPW, class Graph>
Alignment po_poa(const Graph& graph1, const Graph& graph2,
                 const std::vector<uint64_t>& sources1,
                 const std::vector<uint64_t>& sources2,
                 const std::vector<uint64_t>& sinks1,
                 const std::vector<uint64_t>& sinks2,
                 const AlignmentParameters<NumPW>& params);


/*
 * Template instantiations
 */


using IntDP = int32_t;

template<int NumPW>
struct cell_t {
    static const IntDP mininf = std::numeric_limits<IntDP>::min();
    cell_t() : M(mininf) {
        for (int i = 0; i < NumPW; ++i) {
            I[i] = mininf;
            D[i] = mininf;
        }
    }
    ~cell_t() = default;
    IntDP M;
    std::array<IntDP, NumPW> I;
    std::array<IntDP, NumPW> D;
};

template<int NumPW, class Graph>
Alignment po_poa(const Graph& graph1, const Graph& graph2,
                 const std::vector<uint64_t>& sources1,
                 const std::vector<uint64_t>& sources2,
                 const std::vector<uint64_t>& sinks1,
                 const std::vector<uint64_t>& sinks2,
                 const AlignmentParameters<NumPW>& params) {
    
    
    auto order1 = topological_order(graph1);
    auto order2 = topological_order(graph2);
    
    std::vector<std::vector<cell_t<NumPW>>> dp(graph1.node_size(),
                                               std::vector<cell_t<NumPW>>(graph2.node_size()));
    
    // initial conditions
    for (auto node_id1 : sources1) {
        for (auto node_id2 : sources2) {
            dp[node_id1][node_id2].M = (graph1.label(node_id1) == graph2.label(node_id2) ? params.match : -params.mismatch);
        }
    }
    
    for (auto node_id1 : order1) {
        for (auto node_id2 : order2) {
            
            auto& cell = dp[node_id1][node_id2];
            
            // insertions
            for (auto prev_id1 : graph1.previous(node_id1)) {
                auto& prev_cell = dp[prev_id1][node_id2];
                
                for (int pw = 0; pw < NumPW; ++pw) {
                    cell.I[pw] = std::max<IntDP>(cell.I[pw],
                                                 std::max<IntDP>(prev_cell.M - params.gap_open[pw] - params.gap_extend[pw],
                                                                 prev_cell.I[pw] - params.gap_extend[pw]));
                }
            }
            
            // deletions
            for (auto prev_id2 : graph2.previous(node_id2)) {
                auto& prev_cell = dp[node_id1][prev_id2];
                
                for (int pw = 0; pw < NumPW; ++pw) {
                    cell.D[pw] = std::max<IntDP>(cell.D[pw],
                                                 std::max<IntDP>(prev_cell.M - params.gap_open[pw] - params.gap_extend[pw],
                                                                 prev_cell.D[pw] - params.gap_extend[pw]));
                }
            }
            
            // matches
            for (auto prev_id1 : graph1.previous(node_id1)) {
                for (auto prev_id2 : graph2.previous(node_id2)) {
                    auto& prev_cell = dp[prev_id1][prev_id2];
                    cell.M = std::max(cell.M, prev_cell.M);
                }
            }
            cell.M += (graph1.label(node_id1) == graph2.label(node_id2) ? params.match : -params.mismatch);
            
            // choose opt
            for (int pw = 0; pw < NumPW; ++pw) {
                cell.M = std::max(cell.M, std::max(cell.I[pw], cell.D[pw]));
            }
        }
    }
    
    // find the global opt
    uint64_t tb_node1 = -1, tb_node2 = -1;
    for (auto node_id1 : sinks1) {
        for (auto node_id2 : sinks2) {
            if (tb_node1 == -1 || dp[node_id1][node_id2].M > dp[tb_node1][tb_node2].M) {
                tb_node1 = node_id1;
                tb_node2 = node_id2;
            }
        }
    }
    int tb_comp = 0; // positive numbers for insertions, negative for deletions
    
    Alignment alignment;
    
    // do traceback to find alignment
    while (tb_node1 != -1 && tb_node2 != -1) {
        auto here1 = tb_node1;
        auto here2 = tb_node2;
        tb_node1 = -1;
        tb_node2 = -1;
        
        const auto& cell = dp[here1][here2];
        if (tb_comp == 0) {
            // check if this is a gap close
            for (int pw = 0; pw < NumPW; ++pw) {
                if (cell.M == cell.I[pw]) {
                    tb_comp = pw;
                    break;
                }
                if (cell.M == cell.D[pw]) {
                    tb_comp = -pw;
                    break;
                }
            }
        }
        
        if (tb_comp == 0) {
            // diagonal
            alignment.emplace_back(here1, here2);
            
            IntDP align_score = (graph1.label(here1) == graph2.label(here2) ? params.match : -params.mismatch);
            for (auto prev1 : graph1.previous(here1)) {
                for (auto prev2 : graph2.previous(here2)) {
                    if (dp[prev1][prev2].M + align_score == cell.M) {
                        tb_node1 = prev1;
                        tb_node2 = prev2;
                        break;
                    }
                }
            }
        }
        else if (tb_comp > 0) {
            // along graph1
            alignment.emplace_back(AlignedPair::gap, here2);
            for (auto prev1 : graph1.previous(here1)) {
                auto& prev_cell = dp[prev1][here2];
                if (cell.I[tb_comp] == prev_cell.M - params.gap_open[tb_comp] - params.gap_extend[tb_comp]) {
                    tb_comp = 0;
                    tb_node1 = prev1;
                    tb_node2 = here2;
                    break;
                }
                if (cell.I[tb_comp] == prev_cell.I[tb_comp] - params.gap_extend[tb_comp]) {
                    tb_node1 = prev1;
                    tb_node2 = here2;
                    break;
                }
            }
        }
        else {
            // along graph2
            alignment.emplace_back(here1, AlignedPair::gap);
            
            for (auto prev2 : graph2.previous(here2)) {
                auto& prev_cell = dp[here1][prev2];
                if (cell.D[-tb_comp] == prev_cell.M - params.gap_open[-tb_comp] - params.gap_extend[-tb_comp]) {
                    tb_comp = 0;
                    tb_node1 = here1;
                    tb_node2 = prev2;
                    break;
                }
                if (cell.D[-tb_comp] == prev_cell.D[-tb_comp] - params.gap_extend[-tb_comp]) {
                    tb_node1 = here1;
                    tb_node2 = prev2;
                    break;
                }
            }
        }
    }
    
    // traceback is constructed in reverse
    std::reverse(alignment.begin(), alignment.end());
    
    return alignment;
}

}

#endif /* centrolign_alignment_hpp */
