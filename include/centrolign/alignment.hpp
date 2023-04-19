#ifndef centrolign_alignment_hpp
#define centrolign_alignment_hpp

#include <vector>
#include <cstdint>
#include <array>
#include <limits>
#include <algorithm>
#include <unordered_set>

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

/*
 * An alignment is a list of aligned pairs and gaps
 */
typedef std::vector<AlignedPair> Alignment;

// translate a subgraph's alignment to the parent's node IDs
void translate(Alignment& alignment,
               const std::vector<uint64_t>& back_translation1,
               const std::vector<uint64_t>& back_translation2);

/*
 * Piecewise-affine gap score parameters
 */
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

// Global Needleman-Wunsch alignment
template<int NumPW>
Alignment align_nw(const std::string& seq1, const std::string& seq2,
                   const AlignmentParameters<NumPW>& params);

/*
 * Template instantiations
 */


using IntDP = int32_t;

template<int NumPW>
struct cell_t {
    static const IntDP mininf = std::numeric_limits<IntDP>::lowest() / 2; // to avoid underflow
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
    
    static const bool debug_popoa = false;
    
    if (debug_popoa) {
        std::cerr << "begining PO-POA\n";
        std::cerr << "sources on 1:\n";
        for (auto s : sources1) {
            std::cerr << ' ' << s;
        }
        std::cerr << '\n';
        std::cerr << "sources on 2:\n";
        for (auto s : sources2) {
            std::cerr << ' ' << s;
        }
        std::cerr << '\n';
        std::cerr << "sinks on 1:\n";
        for (auto s : sinks1) {
            std::cerr << ' ' << s;
        }
        std::cerr << '\n';
        std::cerr << "sinks on 2:\n";
        for (auto s : sinks2) {
            std::cerr << ' ' << s;
        }
        std::cerr << '\n';
    }
    
    // the DP matrix, with an extra final row/column that acts as boundary conditions
    std::vector<std::vector<cell_t<NumPW>>> dp(graph1.node_size() + 1,
                                               std::vector<cell_t<NumPW>>(graph2.node_size() + 1));
    
    
    auto order1 = topological_order(graph1);
    auto order2 = topological_order(graph2);
    
    // initialize in the boundary row/column
    for (auto node_id1 : sources1) {
        // init matches
        for (auto node_id2 : sources2) {
            dp[node_id1][node_id2].M = (graph1.label(node_id1) == graph2.label(node_id2) ? params.match : -params.mismatch);
        }
        // init initial insertions
        for (int pw = 0; pw < NumPW; ++pw) {
            dp[node_id1].back().I[pw] = -params.gap_open[pw] - params.gap_extend[pw];
        }
    }
    for (auto node_id2 : sources2) {
        // init initial deletions
        for (int pw = 0; pw < NumPW; ++pw) {
            dp.back()[node_id2].D[pw] = -params.gap_open[pw] - params.gap_extend[pw];
        }
    }
    
    // DP along initial insertions
    for (auto node_id1 : order1) {
        auto& cell = dp[node_id1].back();
        // find the opt
        for (int pw = 0; pw < NumPW; ++pw) {
            cell.M = std::max(cell.M, cell.I[pw]);
        }
        // extend initial insertion
        for (auto next_id1 : graph1.next(node_id1)) {
            auto& next_cell = dp[next_id1].back();
            for (int pw = 0; pw < NumPW; ++pw) {
                next_cell.I[pw] = std::max<IntDP>(next_cell.I[pw], cell.I[pw] - params.gap_extend[pw]);
            }
        }
        // open initial deletion
        for (auto next_id2 : sources2) {
            auto& next_cell = dp[node_id1][next_id2];
            for (int pw = 0; pw < NumPW; ++pw) {
                next_cell.D[pw] = std::max<IntDP>(next_cell.D[pw],
                                                  cell.M - params.gap_open[pw] - params.gap_extend[pw]);
            }
        }
        // take initial match/mismatch
        for (auto next_id1 : graph1.next(node_id1)) {
            for (auto next_id2 : sources2) {
                auto& next_cell = dp[next_id1][next_id2];
                IntDP align_score = (graph1.label(next_id1) == graph2.label(next_id2) ? params.match : -params.mismatch);
                next_cell.M = std::max(next_cell.M, cell.M + align_score);
            }
        }
    }
    
    // DP along initial deletions
    for (auto node_id2 : order2) {
        // find the opt
        auto& cell = dp.back()[node_id2];
        for (int pw = 0; pw < NumPW; ++pw) {
            cell.M = std::max(cell.M, cell.D[pw]);
        }
        // extend initial deletion
        for (auto next_id2 : graph2.next(node_id2)) {
            auto& next_cell = dp.back()[next_id2];
            for (int pw = 0; pw < NumPW; ++pw) {
                next_cell.D[pw] = std::max<IntDP>(next_cell.D[pw], cell.D[pw] - params.gap_extend[pw]);
            }
        }
        // open initial insertion
        for (auto next_id1 : sources1) {
            auto& next_cell = dp[next_id1][node_id2];
            for (int pw = 0; pw < NumPW; ++pw) {
                next_cell.I[pw] = std::max<IntDP>(next_cell.I[pw],
                                                  cell.M - params.gap_open[pw] - params.gap_extend[pw]);
            }
        }
        // take initial match/mismatch
        for (auto next_id2 : graph2.next(node_id2)) {
            for (auto next_id1 : sources1) {
                auto& next_cell = dp[next_id1][next_id2];
                IntDP align_score = (graph1.label(next_id1) == graph2.label(next_id2) ? params.match : -params.mismatch);
                next_cell.M = std::max(next_cell.M, cell.M + align_score);
            }
        }
    }
    
    
    // DP in the matrix interior
    for (auto node_id1 : order1) {
        for (auto node_id2 : order2) {
            
            auto& cell = dp[node_id1][node_id2];
            // choose opt
            for (int pw = 0; pw < NumPW; ++pw) {
                cell.M = std::max(cell.M, std::max(cell.I[pw], cell.D[pw]));
            }
            
            // extend insertions
            for (auto next_id1 : graph1.next(node_id1)) {
                auto& next_cell = dp[next_id1][node_id2];
                
                for (int pw = 0; pw < NumPW; ++pw) {
                    next_cell.I[pw] = std::max<IntDP>(next_cell.I[pw],
                                                      std::max<IntDP>(cell.M - params.gap_open[pw] - params.gap_extend[pw],
                                                                      cell.I[pw] - params.gap_extend[pw]));
                }
            }
            
            // extend deletions
            for (auto next_id2 : graph2.next(node_id2)) {
                auto& next_cell = dp[node_id1][next_id2];
                
                for (int pw = 0; pw < NumPW; ++pw) {
                    next_cell.D[pw] = std::max<IntDP>(next_cell.D[pw],
                                                      std::max<IntDP>(cell.M - params.gap_open[pw] - params.gap_extend[pw],
                                                                      cell.D[pw] - params.gap_extend[pw]));
                }
            }
            
            // extend matches/mismatches
            for (auto next_id1 : graph1.next(node_id1)) {
                for (auto next_id2 : graph2.next(node_id2)) {
                    auto& next_cell = dp[next_id1][next_id2];
                    IntDP align_score = (graph1.label(next_id1) == graph2.label(next_id2) ? params.match : -params.mismatch);
                    next_cell.M = std::max(next_cell.M, cell.M + align_score);
                }
            }
        }
    }

    if (debug_popoa) {
        std::cerr << "\t";
        for (uint64_t n2 = 0; n2 < graph2.node_size() + 1; ++n2) {
            std::cerr << '\t' << n2;
        }
        std::cerr << "\n\t";
        for (uint64_t n2 = 0; n2 < graph2.node_size() + 1; ++n2) {
            std::cerr << '\t';
            if (n2 < graph2.node_size()) {
                std::cerr << graph2.label(n2);
            }
            else {
                std::cerr << '*';
            }
        }
        std::cerr << '\n';
        for (uint64_t n1 = 0; n1 < graph1.node_size() + 1; ++n1) {
            std::cerr << n1 << '\t';
            if (n1 < graph1.node_size()) {
                std::cerr << graph1.label(n1);
            }
            else {
                std::cerr << '*';
            }
            for (uint64_t n2 = 0; n2 < graph2.node_size() + 1; ++n2) {
                std::cerr << '\t';
                IntDP m = dp[n1][n2].M;
                if (m < std::numeric_limits<IntDP>::lowest() / 4) {
                    std::cerr << '.';
                }
                else {
                    std::cerr << m;
                }
            }
            std::cerr << '\n';
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
    if (debug_popoa) {
        std::cerr << "opt sink is " << tb_node1 << ", " << tb_node2 << " with score " << dp[tb_node1][tb_node2].M << '\n';
    }
    
    Alignment alignment;
    
    std::unordered_set<uint64_t> sources1_set, sources2_set;
    sources1_set.insert(sources1.begin(), sources1.end());
    sources2_set.insert(sources2.begin(), sources2.end());
    
    int tb_comp = 0; // positive numbers for insertions, negative for deletions
    
    // do traceback to find alignment
    while (tb_node1 != -1 && tb_node2 != -1) {
        
        if (debug_popoa) {
            std::cerr << "tb from " << tb_node1 << ", " << tb_node2 << ", comp " << tb_comp << '\n';
        }
        
        auto here1 = tb_node1;
        auto here2 = tb_node2;
        tb_node1 = -1;
        tb_node2 = -1;
        
        const auto& cell = dp[here1][here2];
        if (tb_comp == 0) {
            // check if this is a gap close
            for (int pw = 0; pw < NumPW; ++pw) {
                if (cell.M == cell.I[pw]) {
                    tb_comp = pw + 1;
                    break;
                }
                if (cell.M == cell.D[pw]) {
                    tb_comp = -pw - 1;
                    break;
                }
            }
        }
        
        // figure out which edges we should look at
        std::vector<uint64_t> previous1;
        std::vector<uint64_t> previous2;
        // add graph edges to previous node list, unless we're in the boundary
        if (here1 < graph1.node_size()) {
            previous1 = graph1.previous(here1);
        }
        if (here2 < graph2.node_size()) {
            previous2 = graph2.previous(here2);
        }
        // add edge to the outer boundary
        if (sources1_set.count(here1)) {
            previous1.push_back(graph1.node_size());
        }
        if (sources2_set.count(here2)) {
            previous2.push_back(graph2.node_size());
        }
        
        if (tb_comp == 0) {
            // diagonal
            alignment.emplace_back(here1, here2);
            
            IntDP align_score = (graph1.label(here1) == graph2.label(here2) ? params.match : -params.mismatch);
            for (auto prev1 : previous1) {
                for (auto prev2 : previous2) {
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
            alignment.emplace_back(here1, AlignedPair::gap);
            
            for (auto prev1 : previous1) {
                auto& prev_cell = dp[prev1][here2];
                if (cell.I[tb_comp - 1] == prev_cell.M - params.gap_open[tb_comp - 1] - params.gap_extend[tb_comp - 1]) {
                    tb_comp = 0;
                    tb_node1 = prev1;
                    tb_node2 = here2;
                    break;
                }
                if (cell.I[tb_comp - 1] == prev_cell.I[tb_comp - 1] - params.gap_extend[tb_comp - 1]) {
                    tb_node1 = prev1;
                    tb_node2 = here2;
                    break;
                }
            }
        }
        else {
            // along graph2
            alignment.emplace_back(AlignedPair::gap, here2);
            
            for (auto prev2 : previous2) {
                auto& prev_cell = dp[here1][prev2];
                if (cell.D[-tb_comp - 1] == prev_cell.M - params.gap_open[-tb_comp - 1] - params.gap_extend[-tb_comp - 1]) {
                    tb_comp = 0;
                    tb_node1 = here1;
                    tb_node2 = prev2;
                    break;
                }
                if (cell.D[-tb_comp - 1] == prev_cell.D[-tb_comp - 1] - params.gap_extend[-tb_comp - 1]) {
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

template<int NumPW>
Alignment align_nw(const std::string& seq1, const std::string& seq2,
                   const AlignmentParameters<NumPW>& params) {
    
    std::vector<std::vector<cell_t<NumPW>>> dp(seq1.size() + 1, std::vector<cell_t<NumPW>>(seq2.size() + 1));
    
    // boundary conditions
    dp[0][0].M = 0;
    for (size_t i = 1; i < dp.size(); ++i) {
        for (int pw = 0; pw < NumPW; ++pw) {
            auto& cell = dp[i][0];
            cell.I[pw] = -params.gap_open[pw] - i * params.gap_extend[pw];
            cell.M = std::max(cell.M, cell.I[pw]);
        }
    }
    for (size_t j = 1; j < dp[0].size(); ++j) {
        for (int pw = 0; pw < NumPW; ++pw) {
            auto& cell = dp[0][j];
            cell.D[pw] = -params.gap_open[pw] - j * params.gap_extend[pw];
            cell.M = std::max(cell.M, cell.D[pw]);
        }
    }
    
    for (size_t i = 0; i < seq1.size(); ++i) {
        for (size_t j = 0; j < seq2.size(); ++j) {
            auto& cell = dp[i + 1][j + 1];
            auto& left = dp[i + 1][j];
            auto& up = dp[i][j + 1];
            auto& diag = dp[i][j];
            cell.M = diag.M + (seq1[i] == seq2[j] ? params.match : -params.mismatch);
            for (int pw = 0; pw < NumPW; ++pw) {
                cell.I[pw] = std::max<IntDP>(up.M - params.gap_open[pw] - params.gap_extend[pw],
                                             up.I[pw] - params.gap_open[pw]);
                cell.D[pw] = std::max<IntDP>(left.M - params.gap_open[pw] - params.gap_extend[pw],
                                             left.D[pw] - params.gap_open[pw]);
                cell.M = std::max(cell.M, std::max(cell.I[pw], cell.D[pw]));
            }
        }
    }
    
    Alignment alignment;
    
    size_t i = seq1.size(), j = seq2.size();
    int tb_comp = 0;
    while (i != 0 || j != 0) {
        
        const auto& cell = dp[i][j];
        if (tb_comp == 0) {
            // check if this is a gap close
            for (int pw = 0; pw < NumPW; ++pw) {
                if (cell.M == cell.I[pw]) {
                    tb_comp = pw + 1;
                    break;
                }
                if (cell.M == cell.D[pw]) {
                    tb_comp = -pw - 1;
                    break;
                }
            }
        }
        
        if (tb_comp == 0) {
            // diagonal
            alignment.emplace_back(i - 1, j - 1);
            IntDP align_score = (seq1[i - 1] == seq2[j - 1] ? params.match : -params.mismatch);
            assert(dp[i - 1][j - 1].M + align_score == dp[i][j].M);
            --i;
            --j;
        }
        else if (tb_comp > 0) {
            // along seq1
            alignment.emplace_back(i - 1, AlignedPair::gap);
            if (dp[i][j].I[tb_comp - 1] == dp[i - 1][j].M - params.gap_open[tb_comp - 1] - params.gap_extend[tb_comp - 1]) {
                tb_comp = 0;
            }
            else {
                assert(dp[i][j].I[tb_comp - 1] == dp[i - 1][j].I[tb_comp - 1] - params.gap_extend[tb_comp - 1]);
            }
            --i;
        }
        else {
            // along seq2
            alignment.emplace_back(AlignedPair::gap, j - 1);
            
            if (dp[i][j].D[-tb_comp - 1] == dp[i][j - 1].M - params.gap_open[-tb_comp - 1] - params.gap_extend[-tb_comp - 1]) {
                tb_comp = 0;
            }
            else {
                assert(dp[i][j].D[-tb_comp - 1] == dp[i][j - 1].D[-tb_comp - 1] - params.gap_extend[-tb_comp - 1]);
            }
            --j;
        }
    }
    
    std::reverse(alignment.begin(), alignment.end());
    return alignment;
}

}

#endif /* centrolign_alignment_hpp */
