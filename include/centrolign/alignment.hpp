#ifndef centrolign_alignment_hpp
#define centrolign_alignment_hpp

#include <vector>
#include <cstdint>
#include <array>
#include <limits>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <sstream>
#include <deque>
#include <queue>
#include <stack>
#include <functional>

#include "centrolign/topological_order.hpp"
#include "centrolign/utility.hpp"
#include "centrolign/graph.hpp"
#include "centrolign/minmax_distance.hpp"
#include "centrolign/target_reachability.hpp"
#include "centrolign/superbubble_distance_oracle.hpp"

namespace centrolign {

/*
 * Two nodes aligned to each other, one of which might be a gap
 */
struct AlignedPair {
    AlignedPair() = default;
    ~AlignedPair() = default;
    AlignedPair(uint64_t node_id1, uint64_t node_id2);
    
    static const uint64_t gap; // sentinel for a gap alignment
    
    uint64_t node_id1 = gap;
    uint64_t node_id2 = gap;
    
    bool operator==(const AlignedPair& other) const;
};

/*
 * An alignment is a list of aligned pairs and gaps
 */
typedef std::vector<AlignedPair> Alignment;

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

// reduce to parameters containing only the first TruncPW gap penalties
template<int NumPW, int TruncPW>
AlignmentParameters<TruncPW> truncate_parameters(const AlignmentParameters<NumPW>& params);

// compute the score of an alignment
template<int NumPW, class Graph>
int64_t score_alignment(const Graph& graph1, const Graph& graph2,
                        const Alignment& alignment, const AlignmentParameters<NumPW>& params);

// PO-POA, can crash if either graph cannot reach one of the sinks from
// one of the sources
template<int NumPW, class Graph>
Alignment po_poa(const Graph& graph1, const Graph& graph2,
                 const std::vector<uint64_t>& sources1,
                 const std::vector<uint64_t>& sources2,
                 const std::vector<uint64_t>& sinks1,
                 const std::vector<uint64_t>& sinks2,
                 const AlignmentParameters<NumPW>& params,
                 int64_t* score_out = nullptr);

// Global Needleman-Wunsch alignment
template<int NumPW>
Alignment align_nw(const std::string& seq1, const std::string& seq2,
                   const AlignmentParameters<NumPW>& params);


// forward declarations for two classes that can be used as the backing map
template<int NumPW>
class ArrayBackedMap;
class HashBackedMap;

// graph-to-graph WFA
template<int NumPW, class Graph, class BackingMap = HashBackedMap>
Alignment wfa_po_poa(const Graph& graph1, const Graph& graph2,
                     const std::vector<uint64_t>& sources1,
                     const std::vector<uint64_t>& sources2,
                     const std::vector<uint64_t>& sinks1,
                     const std::vector<uint64_t>& sinks2,
                     const AlignmentParameters<NumPW>& params,
                     int64_t* score_out = nullptr);

// WFA with pruning
template<int NumPW, class Graph, class BackingMap = HashBackedMap>
Alignment pwfa_po_poa(const Graph& graph1, const Graph& graph2,
                      const std::vector<uint64_t>& sources1,
                      const std::vector<uint64_t>& sources2,
                      const std::vector<uint64_t>& sinks1,
                      const std::vector<uint64_t>& sinks2,
                      const AlignmentParameters<NumPW>& params,
                      int64_t prune_limit,
                      int64_t* score_out = nullptr);

template<int NumPW, class Graph, class BackingMap = HashBackedMap>
Alignment deletion_wfa_po_poa(const Graph& short_graph, const Graph& long_graph,
                              const std::vector<uint64_t>& sources_short,
                              const std::vector<uint64_t>& sources_long,
                              const std::vector<uint64_t>& sinks_short,
                              const std::vector<uint64_t>& sinks_long,
                              const AlignmentParameters<NumPW>& params,
                              int64_t* score_out = nullptr);


// translate a subgraph's alignment to the parent's node IDs
void translate(Alignment& alignment,
               const std::vector<uint64_t>& back_translation1,
               const std::vector<uint64_t>& back_translation2);

// cigar with MID ops, gaps in graph1 are called insertions, gaps in
// graph2 are called deletions (i.e. graph1 is "ref")
std::string cigar(const Alignment& alignment);
// cigar with =/X ops instead of M
std::string explicit_cigar(const Alignment& alignment,
                           const std::string& seq1, const std::string& seq2);
template<class BGraph1, class BGraph2>
std::string explicit_cigar(const Alignment& alignment,
                           const BGraph1& graph1, const BGraph2& graph2);

// The pairwise alignment induced on two input sequence from a graph, where
// the node IDs in the AlignedPair's are taken to be sequence indexes
Alignment induced_pairwise_alignment(const BaseGraph& graph, uint64_t path_id1, uint64_t path_id2);






/*
 * Template instantiations
 */



template<int NumPW, int TruncPW>
AlignmentParameters<TruncPW> truncate_parameters(const AlignmentParameters<NumPW>& params) {
    static_assert(TruncPW <= NumPW, "Cannot truncate to larger number of components");
    AlignmentParameters<TruncPW> trunc;
    trunc.match = params.match;
    trunc.mismatch = params.mismatch;
    for (size_t i = 0; i < TruncPW; ++i) {
        trunc.gap_open[i] = params.gap_open[i];
        trunc.gap_extend[i] = params.gap_extend[i];
    }
    return trunc;
}

template<int NumPW, class Graph>
int64_t score_alignment(const Graph& graph1, const Graph& graph2,
                        const Alignment& alignment, const AlignmentParameters<NumPW>& params) {
    
    auto score_gap = [&](size_t len) -> int64_t {
        if (len == 0) {
            return 0;
        }
        int64_t s = params.gap_open[0] + len * params.gap_extend[0];
        for (size_t i = 1; i < params.gap_open.size(); ++i) {
            s = std::min<int64_t>(s, params.gap_open[i] + len * params.gap_extend[i]);
        }
        return -s;
    };
    
    size_t gap_len = 0;
    
    int64_t score = 0;
    for (const auto& aln_pair : alignment) {
        if (aln_pair.node_id1 != AlignedPair::gap && aln_pair.node_id2 != AlignedPair::gap) {
            score += score_gap(gap_len);
            if (graph1.label(aln_pair.node_id1) == graph2.label(aln_pair.node_id2))
            {
                score += params.match;
            }
            else {
                score -= (int64_t) params.mismatch;
            }
            gap_len = 0;
        }
        else {
            ++gap_len;
        }
    }
    
    score += score_gap(gap_len);
    
    return score;
}

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
                 const AlignmentParameters<NumPW>& params,
                 int64_t* score_out) {
    
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
    if (graph1.node_size() != 0 && graph2.node_size() != 0) {
        // among designated sinks
        for (auto node_id1 : sinks1) {
            for (auto node_id2 : sinks2) {
                if (tb_node1 == -1 || dp[node_id1][node_id2].M > dp[tb_node1][tb_node2].M) {
                    tb_node1 = node_id1;
                    tb_node2 = node_id2;
                }
            }
        }
    }
    else if (graph1.node_size() != 0) {
        // in the lead insertion row
        for (auto node_id1 : sinks1) {
            if (tb_node1 == -1 || dp[node_id1][0].M > dp[tb_node1][0].M) {
                tb_node1 = node_id1;
                tb_node2 = 0;
            }
        }
    }
    else if (graph2.node_size() != 0) {
        // in the lead deletion column
        for (auto node_id2 : sinks2) {
            if (tb_node2 == -1 || dp[0][node_id2].M > dp[0][tb_node2].M) {
                tb_node1 = 0;
                tb_node2 = node_id2;
            }
        }
    }
    
    if (debug_popoa) {
        std::cerr << "opt sink is " << tb_node1 << ", " << tb_node2;
        if (tb_node1 != -1 && tb_node2 != -1) {
            std::cerr << " with score " << dp[tb_node1][tb_node2].M;
        }
        std::cerr << '\n';
    }
    
    if (score_out) {
        if (tb_node1 != -1) {
            *score_out = dp[tb_node1][tb_node2].M;
        }
        else {
            *score_out = 0;
        }
    }
    
    std::unordered_set<uint64_t> sources1_set, sources2_set;
    sources1_set.insert(sources1.begin(), sources1.end());
    sources2_set.insert(sources2.begin(), sources2.end());
    
    Alignment alignment;
    
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
                    if (debug_popoa) {
                        std::cerr << "found gap close for insertion\n";
                    }
                    break;
                }
                if (cell.M == cell.D[pw]) {
                    tb_comp = -pw - 1;
                    if (debug_popoa) {
                        std::cerr << "found gap close for deletion\n";
                    }
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
    
    if (debug_popoa) {
        std::cerr << "completed aligment:\n";
        for (const auto& ap : alignment) {
            std::cerr << (int64_t) ap.node_id1 << '\t' << (int64_t) ap.node_id2 << '\n';
        }
    }
    
    return alignment;
}

// convert to WFA style parameters and also reduce with GCD (returns reduction factor as well)
template<int NumPW>
std::pair<AlignmentParameters<NumPW>, uint32_t> to_wfa_params(const AlignmentParameters<NumPW>& params) {
    
    // euclid's algorithm
    std::function<uint32_t(uint32_t, uint32_t)> gcd = [&](uint32_t a, uint32_t b)  {
        if (a < b) {
            std::swap(a, b);
        }
        auto r = a % b;
        if (r == 0) {
            return b;
        }
        else {
            return gcd(b, r);
        }
    };
    
    // convert to equivalent WFA style parameters
    std::pair<AlignmentParameters<NumPW>, uint32_t> to_return;
    auto& wfa_params = to_return.first;
    auto& factor = to_return.second;
    wfa_params.mismatch = 2 * (params.match + params.mismatch);
    factor = wfa_params.mismatch;
    for (int i = 0; i < NumPW; ++i) {
        wfa_params.gap_open[i] = 2 * params.gap_open[i];
        wfa_params.gap_extend[i] = 2 * params.gap_extend[i] + params.match;
        
        factor = gcd(factor, wfa_params.gap_open[i]);
        factor = gcd(factor, wfa_params.gap_extend[i]);
    }
    
    if (factor != 1) {
        // we can divide through by the GCD to reduce the gap between scores
        wfa_params.mismatch /= factor;
        for (int i = 0; i < NumPW; ++i) {
            wfa_params.gap_open[i] /= factor;
            wfa_params.gap_extend[i] /= factor;
        }
    }
    
    return to_return;
}

// wrapper for a matrix
template<int NumPW>
class ArrayBackedMap {
public:
    ArrayBackedMap(size_t size1, size_t size2) {
        mat_size = row_size * (size2 + 1);
        table.resize(mat_size * (size1 + 1),
                     std::tuple<uint64_t, uint64_t, int>(-1, -1, 0));
    }
    ~ArrayBackedMap() = default;
    
    inline bool count(const std::tuple<uint64_t, uint64_t, int>& key) const {
        return std::get<0>(table[index(key)]) != -1;
    }
    
    inline std::tuple<uint64_t, uint64_t, int>& operator[](const std::tuple<uint64_t, uint64_t, int>& key) {
        return table[index(key)];
    }
    
private:
    
    inline size_t index(const std::tuple<uint64_t, uint64_t, int>& key) const {
        return std::get<0>(key) * mat_size + std::get<1>(key) * row_size + std::get<2>(key) + NumPW;
    }
    
    static const size_t row_size = 2 * NumPW + 1;
    size_t mat_size;
    
    std::vector<std::tuple<uint64_t, uint64_t, int>> table;
    
};

// wrapper for hash table
class HashBackedMap {
public:
    HashBackedMap(size_t size1, size_t size2) { }
    ~HashBackedMap() = default;
    
    inline bool count(const std::tuple<uint64_t, uint64_t, int>& key) const {
        return table.count(key);;
    }
    
    inline std::tuple<uint64_t, uint64_t, int>& operator[](const std::tuple<uint64_t, uint64_t, int>& key) {
        return table[key];
    }
    
private:
    
    std::unordered_map<std::tuple<uint64_t, uint64_t, int>, std::tuple<uint64_t, uint64_t, int>> table;

};

// FIXME: does this always find the optimum score? the lengths of the paths are non-constant
// it's possible that it stops early rather than picking up more matches

// returns (-1, -1) unless this iterations completes, then the final nodes
template<class BackingMap, int NumPW, class Graph, class PruneFunc, class UpdateFunc, class NextFunc1, class NextFunc2, class StopFunc>
inline std::pair<uint64_t, uint64_t>
wfa_iteration(std::deque<std::queue<std::tuple<uint64_t, uint64_t, int, uint64_t, uint64_t, int>>>& queue,
              int64_t& queue_min_score, BackingMap& backpointer,
              const Graph& graph1, const Graph& graph2, const AlignmentParameters<NumPW>& wfa_params,
              const PruneFunc& prune_function, const UpdateFunc& update_function,
              const NextFunc1& next_function1, const NextFunc2& next_function2,
              const StopFunc& stop_function) {
    
    static const bool debug = false;
    
    static const std::pair<uint64_t, uint64_t> null(-1, -1);
    
    auto enqueue = [&](uint64_t from_id1, uint64_t from_id2, int from_comp,
                       uint64_t to_id1, uint64_t to_id2, int to_comp, uint64_t penalty) {
        if (debug) {
            std::cerr << '\t' << "enqueue " << to_id1 << ", " << to_id2 << ", comp " << to_comp << " at score " << (queue_min_score + penalty) << '\n';
        }
        while (queue.size() <= penalty) {
            queue.emplace_back();
        }
        queue[penalty].emplace(from_id1, from_id2, from_comp, to_id1, to_id2, to_comp);
    };
    
    // advance to the next empty score bucket
    while (queue.front().empty()) {
        if (debug) {
            std::cerr << "cleared queue at score " << queue_min_score << " advancing in queue with max score " << (queue_min_score + queue.size()) << '\n';
        }
        queue.pop_front();
        ++queue_min_score;
    }
    
    uint64_t from_id1, from_id2, here_id1, here_id2;
    int from_comp, here_comp;
    std::tie(from_id1, from_id2, from_comp, here_id1, here_id2, here_comp) = queue.front().front();
    queue.front().pop();

    auto key = std::make_tuple(here_id1, here_id2, here_comp);
    if (prune_function(key, queue_min_score) || backpointer.count(key)) {
        // we already reached this position at a lower or equal score
        if (debug) {
            std::cerr << "skip " << here_id1 << ", " << here_id2 << ", comp " << here_comp << " because prune? " << prune_function(key, queue_min_score) << ", backpointer? " << backpointer.count(key) << '\n';
        }
        return null;
    }
    
    if (debug) {
        std::cerr << "WFA score " << queue_min_score << " at " << here_id1 << ", " << here_id2 << ", comp " << here_comp << " with backpointer to " << from_id1 << ", " << from_id2 << ", comp " << from_comp << '\n';
    }
    
    update_function(key, queue_min_score);
    
    backpointer[key] = std::make_tuple(from_id1, from_id2, from_comp);
    
    if (stop_function(here_id1, here_id2, here_comp)) {
        if (debug) {
            std::cerr << "hit stop condition at " << here_id1 << ", " << here_id2 << '\n';
        }
        return std::make_pair(here_id1, here_id2);
    }
    
    if (here_comp == 0) {
        // match/mismatch
        for (auto next_id1 : next_function1(here_id1)) {
            for (auto next_id2 : next_function2(here_id2)) {
                uint64_t penalty = graph1.label(next_id1) == graph2.label(next_id2) ? 0 : wfa_params.mismatch;
                enqueue(here_id1, here_id2, here_comp, next_id1, next_id2, 0, penalty);
            }
            // insert open
            for (int i = 0; i < NumPW; ++i) {
                enqueue(here_id1, here_id2, here_comp, next_id1, here_id2, i + 1,
                        wfa_params.gap_open[i] + wfa_params.gap_extend[i]);
            }
        }
        for (auto next_id2 : next_function2(here_id2)) {
            // deletion open
            for (int i = 0; i < NumPW; ++i) {
                enqueue(here_id1, here_id2, here_comp, here_id1, next_id2, -i - 1,
                        wfa_params.gap_open[i] + wfa_params.gap_extend[i]);
            }
        }
    }
    else {
        // gap close
        enqueue(here_id1, here_id2, here_comp, here_id1, here_id2, 0, 0);
        
        if (here_comp > 0) {
            // extend insertion
            for (auto next_id1 : next_function1(here_id1)) {
                enqueue(here_id1, here_id2, here_comp, next_id1, here_id2, here_comp,
                        wfa_params.gap_extend[here_comp - 1]);
            }
        }
        else {
            // extend insertion
            for (auto next_id2 : next_function2(here_id2)) {
                enqueue(here_id1, here_id2, here_comp, here_id1, next_id2, here_comp,
                        wfa_params.gap_extend[-here_comp - 1]);
            }
        }
    }
    
    return null;
}

template<int NumPW>
int64_t convert_wfa_score(const Alignment& alignment, int64_t wfa_score,
                          const AlignmentParameters<NumPW>& params, uint32_t factor) {
    int64_t total_len = 0;
    for (const auto& aln_pair : alignment) {
        if (aln_pair.node_id1 != AlignedPair::gap) {
            ++total_len;
        }
        if (aln_pair.node_id2 != AlignedPair::gap) {
            ++total_len;
        }
    }
    return (int64_t(params.match) * total_len - wfa_score * factor) / 2;
}

template<class BackingMap, class Graph>
Alignment wfa_traceback(BackingMap& backpointer, uint64_t tb_node1, uint64_t tb_node2,
                        const Graph& graph1, const Graph& graph2) {
    
    Alignment alignment;
    
    // traceback using backpointers
    int tb_comp = 0;
    while (tb_node1 != graph1.node_size() || tb_node2 != graph2.node_size()) {
        
        uint64_t next_id1, next_id2;
        int next_comp;
        std::tie(next_id1, next_id2, next_comp) = backpointer[std::make_tuple(tb_node1, tb_node2, tb_comp)];
        
        if (next_id1 != tb_node1 && next_id2 != tb_node2) {
            alignment.emplace_back(tb_node1, tb_node2);
        }
        else if (next_id1 != tb_node1) {
            alignment.emplace_back(tb_node1, AlignedPair::gap);
        }
        else if (next_id2 != tb_node2) {
            alignment.emplace_back(AlignedPair::gap, tb_node2);
        }
        
        tb_node1 = next_id1;
        tb_node2 = next_id2;
        tb_comp = next_comp;
    }
    
    // put the alignment forward
    std::reverse(alignment.begin(), alignment.end());
    
    return alignment;
}

template<class BackingMap, int NumPW, class Graph, class PruneFunc, class UpdateFunc>
Alignment pwfa_po_poa_internal(const Graph& graph1, const Graph& graph2,
                               const std::vector<uint64_t>& sources1,
                               const std::vector<uint64_t>& sources2,
                               const std::vector<uint64_t>& sinks1,
                               const std::vector<uint64_t>& sinks2,
                               const AlignmentParameters<NumPW>& params,
                               const PruneFunc& prune_function,
                               const UpdateFunc& update_function,
                               int64_t* score_out) {
    
    static const bool debug = false;
    
    if (debug) {
        std::cerr << "## new WFA problem in graphs of size " << graph1.node_size() << " and " << graph2.node_size() << "\n";
    }
    
    // convert to equivalent WFA style parameters
    AlignmentParameters<NumPW> wfa_params;
    uint32_t factor;
    std::tie(wfa_params, factor) = to_wfa_params(params);
        
    // records of (node1, node2, component), positive numbers for insertions, negative for deletions
    BackingMap backpointer(graph1.node_size(), graph2.node_size());
    
    // init the queue (the first "from" location is just a placeholder)
    int64_t queue_min_score = 0;
    std::deque<std::queue<std::tuple<uint64_t, uint64_t, int, uint64_t, uint64_t, int>>> queue;
    queue.emplace_back();
    queue.back().emplace(0, 0, 0, graph1.node_size(), graph2.node_size(), 0);
    
    // lambdas to get either source nodes or neighbors
    // TODO: for some reason these aren't recognized as the same type?
    auto get_next1 = [&](uint64_t node_id1) -> const std::vector<uint64_t>& {
        return node_id1 == graph1.node_size() ? sources1 : graph1.next(node_id1);
    };
    auto get_next2 = [&](uint64_t node_id2) -> const std::vector<uint64_t>& {
        return node_id2 == graph2.node_size() ? sources2 : graph2.next(node_id2);
    };
    
    // convert to sets for O(1) membership testing
    std::unordered_set<uint64_t> sink_set1(sinks1.begin(), sinks1.end());
    std::unordered_set<uint64_t> sink_set2(sinks2.begin(), sinks2.end());
    
    auto stop = [&](uint64_t node_id1, uint64_t node_id2, int comp) {
        return ((sink_set1.empty() || sink_set1.count(node_id1)) &&
                (sink_set2.empty() || sink_set2.count(node_id2)) && comp == 0);
    };
    
    std::pair<uint64_t, uint64_t> end(-1, -1);
    while (end == std::pair<uint64_t, uint64_t>(-1, -1)) {
        end = wfa_iteration(queue, queue_min_score, backpointer, graph1, graph2, wfa_params,
                            prune_function, update_function, get_next1, get_next2, stop);
    }
    
    if (debug) {
        std::cerr << "beginning traceback\n";
    }
        
    Alignment alignment = wfa_traceback(backpointer, end.first, end.second, graph1, graph2);
    
    if (score_out) {
        // convert back to conventional score
        *score_out = convert_wfa_score(alignment, queue_min_score, params, factor);
    }
    
    return alignment;
}


template<int NumPW, class Graph, class BackingMap>
Alignment deletion_wfa_po_poa(const Graph& short_graph, const Graph& long_graph,
                              const std::vector<uint64_t>& sources_short,
                              const std::vector<uint64_t>& sources_long,
                              const std::vector<uint64_t>& sinks_short,
                              const std::vector<uint64_t>& sinks_long,
                              const AlignmentParameters<NumPW>& params,
                              int64_t* score_out) {
    
    // convert to equivalent WFA style parameters
    AlignmentParameters<NumPW> wfa_params;
    uint32_t factor;
    std::tie(wfa_params, factor) = to_wfa_params(params);
    
    // figure out the scope of these scoring parameters
    int64_t scope = wfa_params.mismatch;
    for (int i = 0; i < NumPW; ++i) {
        scope = std::max<int64_t>(scope, wfa_params.gap_open[i] + wfa_params.gap_extend[i]);
    }
    
    // records of (node1, node2, component), positive numbers for insertions, negative for deletions
    BackingMap backpointer_fwd(short_graph.node_size(), long_graph.node_size());
    BackingMap backpointer_rev(short_graph.node_size(), long_graph.node_size());

    // we'll use this for efficient distance queries
    SuperbubbleDistanceOracle dist_oracle(long_graph);
    
    // init the queue (the first "from" location is just a placeholder)
    int64_t queue_min_score_fwd = 0, queue_min_score_rev = 0;
    std::deque<std::queue<std::tuple<uint64_t, uint64_t, int, uint64_t, uint64_t, int>>> queue_fwd, queue_rev;
    queue_fwd.emplace_back();
    queue_fwd.back().emplace(0, 0, 0, short_graph.node_size(), long_graph.node_size(), 0);
    queue_rev.emplace_back();
    queue_rev.back().emplace(0, 0, 0, short_graph.node_size(), long_graph.node_size(), 0);
    
    // transition functions for each direction
    auto get_next_short = [&](uint64_t node_id) -> const std::vector<uint64_t>& {
        return node_id == short_graph.node_size() ? sources_short : short_graph.next(node_id);
    };
    auto get_next_long = [&](uint64_t node_id) -> const std::vector<uint64_t>& {
        return node_id == long_graph.node_size() ? sources_long : long_graph.next(node_id);
    };
    auto get_prev_short = [&](uint64_t node_id) -> const std::vector<uint64_t>& {
        return node_id == short_graph.node_size() ? sinks_short : short_graph.previous(node_id);
    };
    auto get_prev_long = [&](uint64_t node_id) -> const std::vector<uint64_t>& {
        return node_id == long_graph.node_size() ? sinks_long : long_graph.previous(node_id);
    };
    
    auto no_prune = [](const std::tuple<uint64_t, uint64_t, int>& pos, int64_t s) { return false; };
    
    // FIXME: how to handle the case of the short graph being empty
    
    // we also need to remember the score at each position short node -> (long node, score)s
    std::unordered_map<uint64_t, std::vector<std::pair<uint64_t, int64_t>>> fwd_score, rev_score;
    
    int64_t stop_score = std::numeric_limits<int64_t>::max();
    
    // functions to check for wavefronts meeting and record info
    // TODO: repetitive, reimplement as wrappers for one function?
    auto update_fwd = [&](const std::tuple<uint64_t, uint64_t, int>& pos, int64_t s) {
        if (std::get<2>(pos) == 0) {
            fwd_score[std::get<0>(pos)].emplace_back(std::get<1>(pos), s);
        }
        if (stop_score == std::numeric_limits<int64_t>::max() &&
            std::get<0>(pos) != short_graph.node_size()) {
            // the wavefronts have met on the short graph
            auto it = rev_score.find(std::get<0>(pos));
            if (it != rev_score.end()) {
                for (const auto& rev_pos : it->second) {
                    if (dist_oracle.min_distance(std::get<1>(pos), rev_pos.first) != -1) {
                        // the wavefronts are also reachable on the long graph
                        stop_score = s + scope;
                    }
                }
            }
        }
    };
    auto update_rev = [&](const std::tuple<uint64_t, uint64_t, int>& pos, int64_t s) {
        if (std::get<2>(pos) == 0) {
            rev_score[std::get<0>(pos)].emplace_back(std::get<1>(pos), s);
        }
        if (stop_score == std::numeric_limits<int64_t>::max() &&
            std::get<0>(pos) != short_graph.node_size()) {
            // the wavefronts have met on the short graph
            auto it = fwd_score.find(std::get<0>(pos));
            if (it != fwd_score.end()) {
                for (const auto& fwd_pos : it->second) {
                    if (dist_oracle.min_distance(fwd_pos.first, std::get<1>(pos)) != -1) {
                        // the wavefronts are also reachable on the long graph
                        stop_score = s + scope;
                    }
                }
            }
        }
    };

    auto stop = [&](uint64_t node_id1, uint64_t node_id2, int comp) {
        return queue_min_score_fwd >= stop_score && queue_min_score_rev >= stop_score;
    };
    
    // do the core WFA iterations
    std::pair<uint64_t, uint64_t> end_fwd(-1, -1), end_rev(-1, -1);
    while (end_fwd == std::pair<uint64_t, uint64_t>(-1, -1) &&
           end_rev == std::pair<uint64_t, uint64_t>(-1, -1)) {
        if (queue_min_score_fwd <= queue_min_score_rev) {
            end_fwd = wfa_iteration(queue_fwd, queue_min_score_fwd, backpointer_fwd,
                                    short_graph, long_graph, wfa_params,
                                    no_prune, update_fwd, get_next_short, get_next_long, stop);
        }
        else {
            end_rev = wfa_iteration(queue_rev, queue_min_score_rev, backpointer_rev,
                                    short_graph, long_graph, wfa_params,
                                    no_prune, update_rev, get_prev_short, get_prev_long, stop);
        }
    }
    
    int64_t opt_score = std::numeric_limits<int64_t>::max();
    uint64_t opt_short_id = -1, opt_long_fwd_id = -1, opt_long_rev_id = -1;
    
    // find the pair of locations with the best score
    for (const auto& fwd_rec : fwd_score) {
        auto it = rev_score.find(fwd_rec.first);
        if (it != rev_score.end()) {
            // TODO: can i avoid the N^2 all pairs behavior?
            for (const auto& fwd_pos : fwd_rec.second) {
                for (const auto& rev_pos : it->second) {
                    size_t dist = dist_oracle.min_distance(fwd_pos.first, rev_pos.first);
                    if (dist != -1) {
                        // the positions are reachable in the long graph
                        int64_t score = wfa_params.gap_open[0] + wfa_params.gap_extend[0] * dist;
                        for (int i = 1; i < NumPW; ++i) {
                            score = std::min<int64_t>(score, wfa_params.gap_open[i] + wfa_params.gap_extend[i] * dist);
                        }
                        score += fwd_pos.second + rev_pos.second;
                        if (score < opt_score) {
                            // this pair is the new opt
                            opt_score = score;
                            opt_short_id = fwd_rec.first;
                            opt_long_fwd_id = fwd_pos.first;
                            opt_long_rev_id = rev_pos.first;
                        }
                    }
                }
            }
        }
    }
    
    Alignment fwd_traceback = wfa_traceback(backpointer_fwd, opt_short_id, opt_long_fwd_id,
                                            short_graph, long_graph);
    Alignment rev_traceback = wfa_traceback(backpointer_rev, opt_short_id, opt_long_rev_id,
                                            short_graph, long_graph);
    auto path_between = shortest_path(long_graph, opt_long_fwd_id, opt_long_rev_id);
    for (size_t i = 1; i + 1 < path_between.size(); ++i) {
        fwd_traceback.emplace_back(AlignedPair::gap, path_between[i]);
    }
    for (auto& aln_pair : ReverseForEachAdapter<Alignment>(rev_traceback)) {
        fwd_traceback.push_back(aln_pair);
    }
    
    if (score_out) {
        *score_out = convert_wfa_score(fwd_traceback, opt_score, params, factor);
    }
    
    return fwd_traceback;
}

template<int NumPW, class Graph, class BackingMap>
Alignment wfa_po_poa(const Graph& graph1, const Graph& graph2,
                     const std::vector<uint64_t>& sources1,
                     const std::vector<uint64_t>& sources2,
                     const std::vector<uint64_t>& sinks1,
                     const std::vector<uint64_t>& sinks2,
                     const AlignmentParameters<NumPW>& params,
                     int64_t* score_out) {
    auto no_prune = [](const std::tuple<uint64_t, uint64_t, int>& pos, int64_t s) { return false; };
    auto no_update = [](const std::tuple<uint64_t, uint64_t, int>& pos, int64_t s) { return; };
    return pwfa_po_poa_internal<BackingMap>(graph1, graph2, sources1, sources2, sinks1, sinks2,
                                            params, no_prune, no_update, score_out);
}


template<int NumPW, class Graph, class BackingMap>
Alignment pwfa_po_poa(const Graph& graph1, const Graph& graph2,
                      const std::vector<uint64_t>& sources1,
                      const std::vector<uint64_t>& sources2,
                      const std::vector<uint64_t>& sinks1,
                      const std::vector<uint64_t>& sinks2,
                      const AlignmentParameters<NumPW>& params,
                      int64_t prune_limit,
                      int64_t* score_out) {
    
    auto dists1 = minmax_distance(graph1, &sources1);
    auto dists2 = minmax_distance(graph2, &sources2);
    auto reachable1 = target_reachability(graph1, sinks1);
    auto reachable2 = target_reachability(graph2, sinks2);
    
    int64_t furthest_distance = std::numeric_limits<int64_t>::min() + prune_limit;
    
    // decide if we're lagging too far behind
    auto prune = [&](const std::tuple<uint64_t, uint64_t, int>& pos, int64_t s) {
        if ((std::get<0>(pos) < graph1.node_size() && !reachable1[std::get<0>(pos)]) ||
            (std::get<1>(pos) < graph2.node_size() && !reachable2[std::get<1>(pos)])) {
            return true;
        }
        auto d1 = std::get<0>(pos) != graph1.node_size() ? dists1[std::get<0>(pos)].second : -1;
        auto d2 = std::get<1>(pos) != graph2.node_size() ? dists2[std::get<1>(pos)].second : -1;
        return d1 + d2 < furthest_distance - prune_limit;
    };

    auto update = [&](const std::tuple<uint64_t, uint64_t, int>& pos, int64_t s) {
        if ((std::get<0>(pos) == graph1.node_size() || reachable1[std::get<0>(pos)]) &&
            (std::get<1>(pos) == graph2.node_size() || reachable2[std::get<1>(pos)])) {
            auto d1 = std::get<0>(pos) != graph1.node_size() ? dists1[std::get<0>(pos)].first : -1;
            auto d2 = std::get<1>(pos) != graph2.node_size() ? dists2[std::get<1>(pos)].first : -1;
            furthest_distance = std::max<int64_t>(furthest_distance, d1 + d2);
        }
    };
    
    return pwfa_po_poa_internal<BackingMap>(graph1, graph2, sources1, sources2, sinks1, sinks2,
                                            params, prune, update, score_out);
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


// cigar with =/X ops instead of M
template<class BGraph1, class BGraph2>
std::string explicit_cigar(const Alignment& alignment,
                           const BGraph1& graph1, const BGraph2& graph2) {
    
    std::stringstream strm;
    
    int curr_len = 0;
    char curr_op = '\0';
    for (const auto& aln_pair : alignment) {
        char op;
        if (aln_pair.node_id1 == AlignedPair::gap) {
            op = 'I';
        }
        else if (aln_pair.node_id2 == AlignedPair::gap) {
            op = 'D';
        }
        else if (graph1.label(aln_pair.node_id1) == graph2.label(aln_pair.node_id2)) {
            op = '=';
        }
        else {
            op = 'X';
        }
        
        if (op == curr_op) {
            ++curr_len;
        }
        else {
            if (curr_len != 0) {
                strm << curr_len << curr_op;
            }
            curr_len = 1;
            curr_op = op;
        }
    }
    
    if (curr_len != 0) {
        strm << curr_len << curr_op;
    }
    
    return strm.str();
}

}

#endif /* centrolign_alignment_hpp */
