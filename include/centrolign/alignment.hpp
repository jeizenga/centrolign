#ifndef centrolign_alignment_hpp
#define centrolign_alignment_hpp

#include <vector>
#include <cstdint>
#include <array>
#include <limits>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <sstream>
#include <deque>
#include <queue>
#include <stack>
#include <functional>
#include <memory>

#include "centrolign/topological_order.hpp"
#include "centrolign/utility.hpp"
#include "centrolign/graph.hpp"
#include "centrolign/minmax_distance.hpp"
#include "centrolign/target_reachability.hpp"
#include "centrolign/superbubble_distance_oracle.hpp"
#include "centrolign/shortest_path.hpp"
#include "centrolign/orthogonal_max_search_tree.hpp"
#include "centrolign/range_min_query.hpp"

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

// Myers' O(ND) edit distance alignment, O(D) space variant
template<class StringLike>
Alignment align_ond(const StringLike& seq1, const StringLike& seq2);

// Hunt-Szymanski longest common subsequence alignment
template<class StringLike>
Alignment align_hs(const StringLike& seq1, const StringLike& seq2);


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

// WFA with the assumption that there is a long deletion of one graph vs the other
template<int NumPW, class Graph, class BackingMap = HashBackedMap>
Alignment deletion_wfa_po_poa(const Graph& short_graph, const Graph& long_graph,
                              const std::vector<uint64_t>& sources_short,
                              const std::vector<uint64_t>& sources_long,
                              const std::vector<uint64_t>& sinks_short,
                              const std::vector<uint64_t>& sinks_long,
                              const AlignmentParameters<NumPW>& params,
                              int64_t* score_out = nullptr);

// find the lowest cost deletion of entire graph, which is treated as graph1
// in the alignment
template<int NumPW, class Graph>
Alignment pure_deletion_alignment(const Graph& graph,
                                  const std::vector<uint64_t>& sources,
                                  const std::vector<uint64_t>& sinks,
                                  const AlignmentParameters<NumPW>& params,
                                  int64_t* score_out = nullptr);

// greedily take exact matches from each end and then add a double deletion
// for the parts that are not matched
template<int NumPW, class Graph>
Alignment greedy_partial_alignment(const Graph& graph1, const Graph& graph2,
                                   const std::vector<uint64_t>& sources1,
                                   const std::vector<uint64_t>& sources2,
                                   const std::vector<uint64_t>& sinks1,
                                   const std::vector<uint64_t>& sinks2,
                                   const AlignmentParameters<NumPW>& params,
                                   int64_t* score_out = nullptr);


// translate a subgraph's alignment to the parent's node IDs
void translate(Alignment& alignment,
               const std::vector<uint64_t>& back_translation1,
               const std::vector<uint64_t>& back_translation2);

// swap the graph1/graph2 label in the alignment
void swap_graphs(Alignment& alignment);

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

// A set of colinear pairwise alignments induced from a cyclic graph, where
// the node IDs in the AlignedPair's are taken to be sequence indexes
std::vector<Alignment> induced_cyclic_pairwise_alignment(const BaseGraph& graph, uint64_t path_id1, uint64_t path_id2);

void output_maf(std::ostream& out, const std::vector<Alignment>& blocks, const BaseGraph& graph, uint64_t path_id1, uint64_t path_id2);

template<int NumPW>
int64_t rescore(const Alignment& aln, const BaseGraph& graph1, const BaseGraph& graph2,
                const AlignmentParameters<NumPW>& params, bool wfa_style);


// Extend the intervals of the path to close as many gaps as possible between the intervals
// subject to the constraint that none of the intervals contains a cycle
std::vector<std::pair<size_t, size_t>>
maximum_noncyclic_extension(const std::vector<uint64_t>& path,
                            const std::vector<std::pair<size_t, size_t>>& intervals);

// A long common subsequence of subintervals that do not contain any repeating values
// TODO: I have not found an efficient algorithm that guarantees longest yet
template<class StringLike>
Alignment long_common_subsequence_nonrepeating(const StringLike& str1, const StringLike& str2);

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

template<bool Forward, class StringLike>
void ond_next(const StringLike& seq1, const StringLike& seq2,
              size_t begin1, size_t end1, size_t begin2, size_t end2,
              std::vector<size_t>& next, std::vector<size_t>& prev, size_t iter) {
    
    static const bool debug = false;
    if (debug){
        std::cerr << "entering next routine in " << (Forward ? "forward" : "reverse") << " iteration " << iter << "\n";
    }
    
    int64_t d_begin = (Forward ? (begin1 - begin2) : (end1 - end2)) - iter;
    next.resize(prev.size() + 2, -1);
    
    int64_t incr = Forward ? 1 : -1;
    for (int64_t d_rel = 0; d_rel < prev.size(); ++d_rel) {
        int64_t a = prev[d_rel];
        if (a == -1) {
            continue;
        }
        int64_t d = d_begin + d_rel;
        int64_t i = (a + d) / 2;
        int64_t j = (a - d) / 2;
        bool inside1 = (i + incr >= begin1 && i + incr <= end1);
        bool inside2 = (j + incr >= begin2 && j + incr <= end2);
//        std::cerr << "extend d " << d << ", d rel " << d_rel << ", a " << a << ", i " << i << ", j " << j << ", in1 " << inside1 << ", in2 " << inside2 << '\n';
        if (inside1) {
            // insertion
            size_t d_idx = Forward ? d_rel + 2 : d_rel;
            if (next[d_idx] == -1 || (Forward && a + incr > next[d_idx]) || (!Forward && a + incr < next[d_rel])) {
                next[d_idx] = a + incr;
            }
        }
        if (inside2) {
            // deletion
            size_t d_idx = Forward ? d_rel : d_rel + 2;
            if (next[d_idx] == -1 || (Forward && a + incr > next[d_idx]) || (!Forward && a + incr < next[d_idx])) {
                next[d_idx] = a + incr;
            }
        }
        if (inside1 && inside2) {
            // mismatch
//            std::cerr << "mismatch into d rel " << d_rel << " of anti diag " << a + 2 * incr << " compared to current " << next[d_rel + 1] << '\n';
//            std::cerr << Forward << " " << (a + 2 * incr < next[d_rel + 1]) << " " << (!Forward && a + 2 * incr < next[d_rel + 1]) << '\n';
            size_t d_idx = d_rel + 1;
            if (next[d_idx] == -1 || (Forward && a + 2 * incr > next[d_idx]) || (!Forward && a + 2 * incr < next[d_idx])) {
                next[d_idx] = a + 2 * incr;
            }
        }
    }
    
    if (debug) {
        for (int64_t d = d_begin - 1; d < d_begin - 1 + (int64_t) next.size(); ++d) {
            if (d != d_begin - 1) {
                std::cerr << '\t';
            }
            std::cerr << d;
        }
        std::cerr << '\n';
        for (size_t i = 0; i < next.size(); ++i) {
            if (i != 0) {
                std::cerr << '\t';
            }
            std::cerr << (int64_t) next[i];
        }
        std::cerr << '\n';
    }
}

template<bool Forward, class StringLike>
void ond_extend(const StringLike& seq1, const StringLike& seq2,
                size_t begin1, size_t end1, size_t begin2, size_t end2,
                std::vector<size_t>& extending, size_t iter) {
    static const bool debug = false;
    if (debug){
        std::cerr << "entering extend routine in " << (Forward ? "forward" : "reverse") << " iteration " << iter << "\n";
    }
    int64_t d_begin = (Forward ? (begin1 - begin2) : (end1 - end2)) - iter;
    for (int64_t d_rel = 0; d_rel < extending.size(); ++d_rel) {
        int64_t d = d_begin + d_rel;
        size_t& a = extending[d_rel];
        if (a == -1) {
            continue;
        }
        // get the coordinates of the next position to match
        int64_t i = (d + a) / 2 + (Forward ? 0 : -1);
        int64_t j = (a - d) / 2 + (Forward ? 0 : -1);
        while (i < end1 && j < end2 &&
               i >= (int64_t) begin1 && j >= (int64_t) begin2 &&
               seq1[i] == seq2[j]) {
            i += (Forward ? 1 : -1);
            j += (Forward ? 1 : -1);
            a += (Forward ? 2 : -2);
        }
    }
    
    if (debug) {
        for (int64_t d = d_begin; d < d_begin + (int64_t) extending.size(); ++d) {
            if (d != d_begin) {
                std::cerr << '\t';
            }
            std::cerr << d;
        }
        std::cerr << '\n';
        for (size_t i = 0; i < extending.size(); ++i) {
            if (i != 0) {
                std::cerr << '\t';
            }
            std::cerr << (int64_t) extending[i];
        }
        std::cerr << '\n';
    }
}

template<class StringLike>
bool ond_meet(const StringLike& seq1, const StringLike& seq2,
              size_t begin1, size_t end1, size_t begin2, size_t end2,
              const std::vector<size_t>& fwd, const std::vector<size_t>& rev,
              size_t fwd_iter, size_t rev_iter, int64_t* diag_out) {
    
    static const bool debug = false;
    if (debug){
        std::cerr << "entering meet routine\n";
    }
    
    int64_t d_begin_fwd = begin1 - begin2 - fwd_iter;
    int64_t d_begin_rev = end1 - end2 - rev_iter;
    for (int64_t d = std::max(d_begin_fwd, d_begin_rev),
         d_end = std::min<int64_t>(d_begin_fwd + fwd.size(), d_begin_rev + rev.size()); d < d_end; ++d) {
        if (fwd[d - d_begin_fwd] >= rev[d - d_begin_rev] && fwd[d - d_begin_fwd] != -1) {
            *diag_out = d;
            if (debug){
                std::cerr << "found meet at diagonal " << d << "\n";
            }
            return true;
        }
    }
    if (debug){
        std::cerr << "did not find meet\n";
    }
    return false;
}

template<bool Forward, class StringLike>
Alignment ond_traceback_middle(const StringLike& seq1, const StringLike& seq2,
                               size_t begin1, size_t end1, size_t begin2, size_t end2,
                               const std::vector<size_t>& dp, const std::vector<size_t>& prev, size_t iter,
                               int64_t diag, int64_t& diag_out, size_t& anti_diag_out) {
    static const bool debug = false;
    Alignment trace;
    
    int64_t d_begin = (Forward ? (begin1 - begin2) : (end1 - end2)) - iter;
    int64_t incr = Forward ? 1 : -1;
    
    int64_t d_rel = diag - d_begin;
    
    if (debug){
        std::cerr << "entering middle traceback routine along " << (Forward ? "forward" : "reverse") << " direction in diagonal " << diag  << " at antidiagonal " << dp[d_rel] << ", relative diag " << d_rel << " in iter " << iter << "\n";
    }
    
    size_t a = dp[d_rel];
    size_t i = (a + diag) / 2 - (Forward ? 1 : 0);
    size_t j = (a - diag) / 2 - (Forward ? 1 : 0);
    if (debug){
        std::cerr << "start trace from i " << i << ", j " << j << ", a " << a << '\n';
    }
    while (Forward ? a > begin1 + begin2 : a < end1 + end2) {
        
        if (d_rel >= 2) {
            // check insertion
            if (prev[d_rel - 2] != -1 && a == prev[d_rel - 2] + incr) {
                a -= incr;
                diag -= 1;
                if (Forward) {
                    trace.emplace_back(i, AlignedPair::gap);
                }
                else {
                    trace.emplace_back(AlignedPair::gap, j);
                }
                break;
            }
        }
        if (d_rel >= 1 && d_rel + 1 < dp.size()) {
            // check mismatch
            if (prev[d_rel - 1] != -1 && a == prev[d_rel - 1] + 2 * incr) {
                a -= 2 * incr;
                trace.emplace_back(i, j);
                break;
            }
        }
        if (d_rel + 2 < dp.size()) {
            // check deletion
            if (prev[d_rel] != -1 && a == prev[d_rel] + incr) {
                a -= incr;
                diag += 1;
                if (Forward) {
                    trace.emplace_back(AlignedPair::gap, j);
                }
                else {
                    trace.emplace_back(i, AlignedPair::gap);
                }
                break;
            }
        }
        
        trace.emplace_back(i, j);
        
        a -= 2 * incr;
        i -= incr;
        j -= incr;
        if (debug){
            std::cerr << "move to i " << i << ", j " << j << ", a " << a << '\n';
        }
    }
    
    if (debug) {
        std::cerr << "complete traceback at d " << diag << ", a " << a << '\n';
    }
    
    diag_out = diag;
    anti_diag_out = a;
    
    if (Forward) {
        std::reverse(trace.begin(), trace.end());
    }
    
    if (debug) {
        std::cerr << "traceback:\n";
        for (auto aln_pair : trace) {
            std::cerr << '\t' << (int64_t) aln_pair.node_id1 << ", " << (int64_t) aln_pair.node_id2 << '\n';
        }
    }
    
    return trace;
}

template<class StringLike>
void align_ond_internal(const StringLike& seq1, const StringLike& seq2,
                        size_t begin1, size_t end1, size_t begin2, size_t end2,
                        Alignment& partial_aln) {
    
    static const bool debug = false;
    static const bool print_progress = false;
    if (debug || (print_progress && seq1.size() + seq2.size() > 1000)){
        std::cerr << "recursive O(ND) call in range [" << begin1 << ", " << end1 << ") x [" << begin2 << ", " << end2 << ")\n";
    }
    // handle these cases, so we don't need to worry about them in the traceback
    // TODO: i no longer *need* this, but i still think it might be a good optimization
    if (begin1 == end1) {
        assert(begin2 != end2);
        for (size_t i = begin2; i < end2; ++i) {
            partial_aln.emplace_back(AlignedPair::gap, i);
        }
        return;
    }
    if (begin2 == end2) {
        assert(begin1 != end1);
        for (size_t i = begin1; i < end1; ++i) {
            partial_aln.emplace_back(i, AlignedPair::gap);
        }
        return;
    }
    
    // init the DP
    std::vector<size_t> fwd_row, rev_row, prev_fwd_row, prev_rev_row;
    size_t fwd_iter = 0, rev_iter = 0;
    fwd_row.emplace_back(begin1 + begin2);
    rev_row.emplace_back(end1 + end2);
    ond_extend<true>(seq1, seq2, begin1, end1, begin2, end2, fwd_row, fwd_iter);
    ond_extend<false>(seq1, seq2, begin1, end1, begin2, end2, rev_row, rev_iter);
    
    // next and extend iterations until we get an overlap
    int64_t meet_diag;
    bool forward = true;
    while (!ond_meet(seq1, seq2, begin1, end1, begin2, end2,
                     fwd_row, rev_row, fwd_iter, rev_iter, &meet_diag)) {
        
        // hack-y code to track progress in very large alignments
        if (print_progress && seq1.size() + seq2.size() > 1000 && (fwd_iter + rev_iter) % 20000 == 0) {
            int64_t d_begin_fwd = (begin1 - begin2) - fwd_iter;
            int64_t d_begin_rev = (end1 - end2) - rev_iter;
            size_t fwd_progress = 0;
            for (size_t idx = 0; idx < fwd_row.size(); ++idx) {
                int64_t d = d_begin_fwd + idx;
                int64_t a = fwd_row[idx];
                size_t i = (a + d) / 2;
                size_t j = (a - d) / 2;
                fwd_progress = std::max(std::min(i - begin1, j - begin2), fwd_progress);
            }
            size_t rev_progress = 0;
            for (size_t idx = 0; idx < rev_row.size(); ++idx) {
                int64_t d = d_begin_rev + idx;
                int64_t a = rev_row[idx];
                size_t i = (a + d) / 2;
                size_t j = (a - d) / 2;
                rev_progress = std::max(std::min(end1 - i, end2 - j), rev_progress);
            }
            std::cerr << "iter " << (fwd_iter + rev_iter) << ", fwd progress " << fwd_progress << ", rev progress " << rev_progress << ", approx completion " << double(fwd_progress + rev_progress) / double(std::min(end1 - begin1, end2 - begin2)) << '\n';
            
        }
        
        if (forward) {
            swap(fwd_row, prev_fwd_row);
            fwd_row.clear();
            ond_next<true>(seq1, seq2, begin1, end1, begin2, end2,
                           fwd_row, prev_fwd_row, fwd_iter);
            ++fwd_iter;
            ond_extend<true>(seq1, seq2, begin1, end1, begin2, end2,
                             fwd_row, fwd_iter);
        }
        else {
            swap(rev_row, prev_rev_row);
            rev_row.clear();
            ond_next<false>(seq1, seq2, begin1, end1, begin2, end2,
                            rev_row, prev_rev_row, rev_iter);
            ++rev_iter;
            ond_extend<false>(seq1, seq2, begin1, end1, begin2, end2,
                              rev_row, rev_iter);
        }
        forward = !forward;
    }
    
    Alignment middle_aln;
    size_t trace_begin1, trace_begin2, trace_end1, trace_end2;
    if (!forward) {
        // the last iteration was forward
        int64_t diag_trace_end;
        size_t anti_diag_trace_end;
        middle_aln = std::move(ond_traceback_middle<true>(seq1, seq2, begin1, end1, begin2, end2, fwd_row,
                                                          prev_fwd_row, fwd_iter, meet_diag, diag_trace_end,
                                                          anti_diag_trace_end));
        int64_t diag_begin = begin1 - begin2 - fwd_iter;
        trace_begin1 = (anti_diag_trace_end + diag_trace_end) / 2;
        trace_begin2 = (anti_diag_trace_end - diag_trace_end) / 2;
        trace_end1 = (fwd_row[meet_diag - diag_begin] + meet_diag) / 2;
        trace_end2 = (fwd_row[meet_diag - diag_begin] - meet_diag) / 2;
    }
    else {
        // the last iteration was reverse
        int64_t diag_trace_end;
        size_t anti_diag_trace_end;
        middle_aln = std::move(ond_traceback_middle<false>(seq1, seq2, begin1, end1, begin2, end2, rev_row,
                                                           prev_rev_row, rev_iter, meet_diag, diag_trace_end,
                                                           anti_diag_trace_end));
        int64_t diag_begin = end1 - end2 - rev_iter;
        trace_begin1 = (rev_row[meet_diag - diag_begin] + meet_diag) / 2;
        trace_begin2 = (rev_row[meet_diag - diag_begin] - meet_diag) / 2;
        trace_end1 = (anti_diag_trace_end + diag_trace_end) / 2;
        trace_end2 = (anti_diag_trace_end - diag_trace_end) / 2;
    }
    
    if (debug) {
        std::cerr << "traceback covered interval [" << trace_begin1 << ", " << trace_end1 << ") x [" << trace_begin2 << ", " << trace_end2 << ")\n";
    }
    
    if (trace_begin1 != begin1 || trace_begin2 != begin2) {
        align_ond_internal(seq1, seq2, begin1, trace_begin1, begin2, trace_begin2, partial_aln);
    }
    for (auto& aln_pair : middle_aln) {
        partial_aln.push_back(aln_pair);
    }
    if (trace_end1 < end1 || trace_end2 < end2) {
        align_ond_internal(seq1, seq2, trace_end1, end1, trace_end2, end2, partial_aln);
    }
}

template<class StringLike>
Alignment align_ond(const StringLike& seq1, const StringLike& seq2) {
    
    Alignment alignment;
    
    align_ond_internal(seq1, seq2, 0, seq1.size(), 0, seq2.size(), alignment);
    return alignment;
}

template<class StringLike>
Alignment align_hs(const StringLike& seq1, const StringLike& seq2) {
    
    static const bool debug = false;
    
    // find the active indexes in the matrix
    std::vector<std::vector<size_t>> active_indexes(seq1.size() + 1);
    {
        // bin seq2 occurrences by value (in reverse order)
        std::vector<std::vector<size_t>> occurrences;
        for (int64_t j = seq2.size() - 1; j >= 0; --j) {
            while (occurrences.size() <= seq2[j]) {
                occurrences.emplace_back();
            }
            occurrences[seq2[j]].push_back(j + 1);
        }
        
        // find matches using bins
        for (size_t i = 0; i < seq1.size(); ++i) {
            if (seq1[i] < occurrences.size()) {
                active_indexes[i + 1] = occurrences[seq1[i]];
            }
        }
    }
    
    if (debug) {
        std::cerr << "active indexes:\n";
        for (size_t i = 0; i < active_indexes.size(); ++i) {
            std::cerr << i << ":";
            for (auto j : active_indexes[i]) {
                std::cerr << ' ' << j;
            }
            std::cerr << '\n';
        }
    }
    
    // the previous match in the LCS
    std::unordered_map<std::pair<size_t, size_t>, std::pair<size_t, size_t>> backpointer;
    // for each score, the coordinate that achieves
    std::vector<std::pair<size_t, size_t>> score_heads(1, std::pair<size_t, size_t>(0, 0));
    
    // init the sparse description of the row
    std::vector<size_t> row(1, 0);
    // sparse DP
    for (size_t i = 1; i <= seq1.size(); ++i) {
        for (auto j : active_indexes[i]) {
            auto it = std::lower_bound(row.begin(), row.end(), j);
            if (it == row.end()) {
                // this match makes the highest LCS we've seen so far
                backpointer[std::make_pair(i, j)] = score_heads.back();
                score_heads.emplace_back(i, j);
                row.push_back(j);
            }
            else if (*it != j) {
                // we've found match that makes the same LCS but at an earlier index
                backpointer[std::make_pair(i, j)] = score_heads[it - row.begin() - 1];
                score_heads[it - row.begin()] = std::make_pair(i, j);
                *it = j;
            }
        }
    }
    
    if (debug) {
        std::vector<std::pair<size_t, size_t>> keys;
        for (const auto& r : backpointer) {
            keys.push_back(r.first);
        }
        std::sort(keys.begin(), keys.end());
        std::cerr << "backpointers:\n";
        for (auto k : keys) {
            std::cerr << '\t' << k.first << ',' << k.second << ": " << backpointer[k].first << ',' << backpointer[k].first << '\n';
        }
    }
    
    // follow back pointers to do traceback
    Alignment traceback;
    for (size_t i = seq1.size(); i > score_heads.back().first; --i) {
        traceback.emplace_back(i - 1, AlignedPair::gap);
    }
    for (size_t j = seq2.size(); j > score_heads.back().second; --j) {
        traceback.emplace_back(AlignedPair::gap, j - 1);
    }
    auto here = score_heads.back();
    while (backpointer.count(here)) {
        traceback.emplace_back(here.first - 1, here.second - 1);
        auto next =  backpointer[here];
        for (size_t i = here.first - 1; i > next.first; --i) {
            traceback.emplace_back(i - 1, AlignedPair::gap);
        }
        for (size_t j = here.second - 1; j > next.second; --j) {
            traceback.emplace_back(AlignedPair::gap, j - 1);
        }
        here = next;
    }
    
    // put traceback in forward order
    std::reverse(traceback.begin(), traceback.end());
    
    return traceback;
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

template<bool Maximize, int NumPW, class Graph>
Alignment po_poa_internal(const Graph& graph1, const Graph& graph2,
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
    
    if (!Maximize) {
        // reset initial state to +infinity
        static const IntDP inf = std::numeric_limits<IntDP>::max() / 2;
        for (auto& row : dp) {
            for (auto& entry : row) {
                entry.M = inf;
                for (int i = 0; i < NumPW; ++i) {
                    entry.I[i] = inf;
                    entry.D[i] = inf;
                }
            }
        }
    }
    
    auto order1 = topological_order(graph1);
    auto order2 = topological_order(graph2);
    
    constexpr IntDP sign = (Maximize ? 1 : -1);
    const auto better = [](IntDP val1, IntDP val2) -> IntDP { return Maximize ? std::max(val1, val2) : std::min(val1, val2); };
    const auto is_better = [](IntDP val1, IntDP val2) -> bool { return Maximize ? (val1 > val2) : (val1 < val2); };
    
    // initialize in the boundary row/column
    for (auto node_id1 : sources1) {
        // init matches
        for (auto node_id2 : sources2) {
            dp[node_id1][node_id2].M = (graph1.label(node_id1) == graph2.label(node_id2) ? params.match : -params.mismatch) * sign;
        }
        // init initial insertions
        for (int pw = 0; pw < NumPW; ++pw) {
            dp[node_id1].back().I[pw] = (-params.gap_open[pw] - params.gap_extend[pw]) * sign;
        }
    }
    for (auto node_id2 : sources2) {
        // init initial deletions
        for (int pw = 0; pw < NumPW; ++pw) {
            dp.back()[node_id2].D[pw] = (-params.gap_open[pw] - params.gap_extend[pw]) * sign;
        }
    }
    
    
    // DP along initial insertions
    for (auto node_id1 : order1) {
        auto& cell = dp[node_id1].back();
        // find the opt
        for (int pw = 0; pw < NumPW; ++pw) {
            cell.M = better(cell.M, cell.I[pw]);
        }
        // extend initial insertion
        for (auto next_id1 : graph1.next(node_id1)) {
            auto& next_cell = dp[next_id1].back();
            for (int pw = 0; pw < NumPW; ++pw) {
                next_cell.I[pw] = better(next_cell.I[pw], cell.I[pw] - params.gap_extend[pw] * sign);
            }
        }
        // open initial deletion
        for (auto next_id2 : sources2) {
            auto& next_cell = dp[node_id1][next_id2];
            for (int pw = 0; pw < NumPW; ++pw) {
                next_cell.D[pw] = better(next_cell.D[pw],
                                         cell.M - (params.gap_open[pw] + params.gap_extend[pw]) * sign);
            }
        }
        // take initial match/mismatch
        for (auto next_id1 : graph1.next(node_id1)) {
            for (auto next_id2 : sources2) {
                auto& next_cell = dp[next_id1][next_id2];
                IntDP align_score = (graph1.label(next_id1) == graph2.label(next_id2) ? params.match : -params.mismatch) * sign;
                next_cell.M = better(next_cell.M, cell.M + align_score);
            }
        }
    }
    
    // DP along initial deletions
    for (auto node_id2 : order2) {
        // find the opt
        auto& cell = dp.back()[node_id2];
        for (int pw = 0; pw < NumPW; ++pw) {
            cell.M = better(cell.M, cell.D[pw]);
        }
        // extend initial deletion
        for (auto next_id2 : graph2.next(node_id2)) {
            auto& next_cell = dp.back()[next_id2];
            for (int pw = 0; pw < NumPW; ++pw) {
                next_cell.D[pw] = better(next_cell.D[pw], cell.D[pw] - params.gap_extend[pw] * sign);
            }
        }
        // open initial insertion
        for (auto next_id1 : sources1) {
            auto& next_cell = dp[next_id1][node_id2];
            for (int pw = 0; pw < NumPW; ++pw) {
                next_cell.I[pw] = better(next_cell.I[pw],
                                         cell.M - (params.gap_open[pw] + params.gap_extend[pw]) * sign);
            }
        }
        // take initial match/mismatch
        for (auto next_id2 : graph2.next(node_id2)) {
            for (auto next_id1 : sources1) {
                auto& next_cell = dp[next_id1][next_id2];
                IntDP align_score = (graph1.label(next_id1) == graph2.label(next_id2) ? params.match : -params.mismatch) * sign;
                next_cell.M = better(next_cell.M, cell.M + align_score);
            }
        }
    }
    
    
    // DP in the matrix interior
    for (auto node_id1 : order1) {
        for (auto node_id2 : order2) {
            
            auto& cell = dp[node_id1][node_id2];
            // choose opt
            for (int pw = 0; pw < NumPW; ++pw) {
                cell.M = better(cell.M, better(cell.I[pw], cell.D[pw]));
            }
            
            // extend insertions
            for (auto next_id1 : graph1.next(node_id1)) {
                auto& next_cell = dp[next_id1][node_id2];
                
                for (int pw = 0; pw < NumPW; ++pw) {
                    next_cell.I[pw] = better(next_cell.I[pw],
                                             better(cell.M - (params.gap_open[pw] + params.gap_extend[pw]) * sign,
                                                    cell.I[pw] - params.gap_extend[pw] * sign));
                }
            }
            
            // extend deletions
            for (auto next_id2 : graph2.next(node_id2)) {
                auto& next_cell = dp[node_id1][next_id2];
                
                for (int pw = 0; pw < NumPW; ++pw) {
                    next_cell.D[pw] = better(next_cell.D[pw],
                                             better(cell.M - (params.gap_open[pw] + params.gap_extend[pw]) * sign,
                                                    cell.D[pw] - params.gap_extend[pw] * sign));
                }
            }
            
            // extend matches/mismatches
            for (auto next_id1 : graph1.next(node_id1)) {
                for (auto next_id2 : graph2.next(node_id2)) {
                    auto& next_cell = dp[next_id1][next_id2];
                    IntDP align_score = (graph1.label(next_id1) == graph2.label(next_id2) ? params.match : -params.mismatch) * sign;
                    next_cell.M = better(next_cell.M, cell.M + align_score);
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
                if (tb_node1 == -1 || is_better(dp[node_id1][node_id2].M, dp[tb_node1][tb_node2].M)) {
                    tb_node1 = node_id1;
                    tb_node2 = node_id2;
                }
            }
        }
    }
    else if (graph1.node_size() != 0) {
        // in the lead insertion row
        for (auto node_id1 : sinks1) {
            if (tb_node1 == -1 || is_better(dp[node_id1][0].M, dp[tb_node1][0].M)) {
                tb_node1 = node_id1;
                tb_node2 = 0;
            }
        }
    }
    else if (graph2.node_size() != 0) {
        // in the lead deletion column
        for (auto node_id2 : sinks2) {
            if (tb_node2 == -1 || is_better(dp[0][node_id2].M, dp[0][tb_node2].M)) {
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
            
            IntDP align_score = (graph1.label(here1) == graph2.label(here2) ? params.match : -params.mismatch) * sign;
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
                if (cell.I[tb_comp - 1] == prev_cell.M - (params.gap_open[tb_comp - 1] + params.gap_extend[tb_comp - 1]) * sign) {
                    tb_comp = 0;
                    tb_node1 = prev1;
                    tb_node2 = here2;
                    break;
                }
                if (cell.I[tb_comp - 1] == prev_cell.I[tb_comp - 1] - params.gap_extend[tb_comp - 1] * sign) {
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
                if (cell.D[-tb_comp - 1] == prev_cell.M - (params.gap_open[-tb_comp - 1] + params.gap_extend[-tb_comp - 1]) * sign) {
                    tb_comp = 0;
                    tb_node1 = here1;
                    tb_node2 = prev2;
                    break;
                }
                if (cell.D[-tb_comp - 1] == prev_cell.D[-tb_comp - 1] - params.gap_extend[-tb_comp - 1] * sign) {
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

template<int NumPW, class Graph>
Alignment po_poa(const Graph& graph1, const Graph& graph2,
                 const std::vector<uint64_t>& sources1,
                 const std::vector<uint64_t>& sources2,
                 const std::vector<uint64_t>& sinks1,
                 const std::vector<uint64_t>& sinks2,
                 const AlignmentParameters<NumPW>& params,
                 int64_t* score_out) {
    // dispatch to internal function
    return po_poa_internal<true>(graph1, graph2, sources1, sources2, sinks1, sinks2, params, score_out);
}

template<int NumPW, class Graph>
Alignment min_po_poa(const Graph& graph1, const Graph& graph2,
                     const std::vector<uint64_t>& sources1,
                     const std::vector<uint64_t>& sources2,
                     const std::vector<uint64_t>& sinks1,
                     const std::vector<uint64_t>& sinks2,
                     const AlignmentParameters<NumPW>& params,
                     int64_t* score_out = nullptr) {
    // dispatch to internal function
    return po_poa_internal<false>(graph1, graph2, sources1, sources2, sinks1, sinks2, params, score_out);
}


template<int NumPW, class Graph>
Alignment pure_deletion_alignment(const Graph& graph,
                                  const std::vector<uint64_t>& sources,
                                  const std::vector<uint64_t>& sinks,
                                  const AlignmentParameters<NumPW>& params,
                                  int64_t* score_out) {
    
    std::vector<uint64_t> path;
    if (graph.node_size() != 0) {
        path = std::move(shortest_path(graph, sources, sinks));
        
    }
    
    Alignment aln;
    aln.reserve(path.size());
    for (auto node_id : path) {
        aln.emplace_back(node_id, AlignedPair::gap);
    }
    
    if (score_out) {
        if (path.size() == 0) {
            *score_out = 0;
        }
        else {
            *score_out = std::numeric_limits<int64_t>::max();
            for (int i = 0; i < NumPW; ++i) {
                *score_out = std::min<int64_t>(*score_out, -params.gap_open[i] - params.gap_extend[i]);
            }
        }
    }
    
    return aln;
}

template<int NumPW, class Graph>
Alignment greedy_partial_alignment(const Graph& graph1, const Graph& graph2,
                                   const std::vector<uint64_t>& sources1,
                                   const std::vector<uint64_t>& sources2,
                                   const std::vector<uint64_t>& sinks1,
                                   const std::vector<uint64_t>& sinks2,
                                   const AlignmentParameters<NumPW>& params,
                                   int64_t* score_out) {
    
    static const bool debug = false;
    if (debug) {
        std::cerr << "beginning greedy partial alignment with graphs of size " << graph1.node_size() << " and " << graph2.node_size() << "\n";
        std::cerr << "sources 1:\n";
        for (auto src : sources1) {
            std::cerr << '\t' << src << '\n';
        }
        std::cerr << "sources 2:\n";
        for (auto src : sources2) {
            std::cerr << '\t' << src << '\n';
        }
        std::cerr << "sinks 1:\n";
        for (auto snk : sinks1) {
            std::cerr << '\t' << snk << '\n';
        }
        std::cerr << "sinks 2:\n";
        for (auto snk : sinks2) {
            std::cerr << '\t' << snk << '\n';
        }
    }
    
    Alignment aln_fwd, aln_rev;
    
    for (bool forward : {true, false}) {
        
        size_t max_path_len = 0;
        std::pair<uint64_t, uint64_t> path_end(-1, -1);
        
        std::unordered_map<std::pair<uint64_t, uint64_t>, std::pair<uint64_t, uint64_t>> backpointer;
        
        // records of (node1, node2, match length)
        std::vector<std::tuple<uint64_t, uint64_t, size_t>> stack;
        for (auto node_id1 : (forward ? sources1 : sinks1)) {
            for (auto node_id2 : (forward ? sources2 : sinks2)) {
                
                if (graph1.label(node_id1) == graph2.label(node_id2)) {
                    stack.emplace_back(node_id1, node_id2, 1);
                    backpointer[std::make_pair(node_id1, node_id2)] = std::pair<uint64_t, uint64_t>(-1, -1);
                }
            }
        }
        
        // DFS through matches in the alignment graph
        while (!stack.empty()) {
            
            uint64_t node_id1, node_id2;
            size_t path_len;
            std::tie(node_id1, node_id2, path_len) = stack.back();
            stack.pop_back();
            
            if (debug) {
                std::cerr << "fwd direction? " << forward << ", node1 " << node_id1 << ", node2 " << node_id2 << ", path length " << path_len << '\n';
            }
            
            if (path_len > max_path_len) {
                max_path_len = path_len;
                path_end = std::make_pair(node_id1, node_id2);
            }
            
            for (auto next_id1 : (forward ? graph1.next(node_id1) : graph1.previous(node_id1))) {
                for (auto next_id2 : (forward ? graph2.next(node_id2) : graph2.previous(node_id2))) {
                    
                    if (graph1.label(next_id1) == graph2.label(next_id2) &&
                        !backpointer.count(std::make_pair(next_id1, next_id2))) {
                        
                        if (debug) {
                            std::cerr << "\tqueue nodes " << next_id1 << ", " << next_id2 << '\n';
                        }
                        
                        backpointer[std::make_pair(next_id1, next_id2)] = std::make_pair(node_id1, node_id2);
                        stack.emplace_back(next_id1, next_id2, path_len + 1);
                    }
                }
            }
        }
        
        // traceback
        Alignment& aln = forward ? aln_fwd : aln_rev;
        while (path_end != std::pair<uint64_t, uint64_t>(-1, -1)) {
            aln.emplace_back(path_end.first, path_end.second);
            path_end = backpointer.at(path_end);
        }
        
        if (forward) {
            // the forward alignment is build backwards
            std::reverse(aln.begin(), aln.end());
        }
    }
    
    if (debug) {
        std::cerr << "forward greedy match:\n";
        for (const auto& aln_pair : aln_fwd) {
            std::cerr << '\t' << (int64_t) aln_pair.node_id1 << '\t' << (int64_t) aln_pair.node_id2 << '\n';
        }
        std::cerr << "reverse greedy match:\n";
        for (const auto& aln_pair : aln_rev) {
            std::cerr << '\t' << (int64_t) aln_pair.node_id1 << '\t' << (int64_t) aln_pair.node_id2 << '\n';
        }
    }
    
    // the distance oracles are pretty memory intensive, so we'll try to avoid them if possible
    // both code paths will fill out these variables
    size_t left_trim = 0;
    size_t right_trim = 0;
    std::vector<uint64_t> shortest_path1;
    std::vector<uint64_t> shortest_path2;
    
    bool found_path = false;
    if (aln_fwd.empty() || aln_rev.empty() ||
        (aln_fwd.back().node_id1 != aln_rev.front().node_id1 && aln_fwd.back().node_id2 != aln_rev.front().node_id2)) {
        // note: these conditions avoid an edge case where a shortest path will be returned despite the fact
        // that alignments overlap
        
        // try a shortest path from the end of the alignment on graph1
        std::vector<uint64_t> shortest_path_start1, shortest_path_end1;
        if (aln_fwd.empty()) {
            shortest_path_start1 = sources1;
        }
        else {
            shortest_path_start1.push_back(aln_fwd.back().node_id1);
        }
        if (aln_rev.empty()) {
            shortest_path_end1 = sinks1;
        }
        else {
            shortest_path_end1.push_back(aln_rev.front().node_id1);
        }
        if (debug) {
            std::cerr << "checking for simple path from graph 1 starts:\n";
            for (auto n : shortest_path_start1) {
                std::cerr << '\t' << n << '\n';
            }
            std::cerr << "and graph 1 ends:\n";
            for (auto n : shortest_path_end1) {
                std::cerr << '\t' << n << '\n';
            }
        }
        if (!shortest_path_start1.empty() && !shortest_path_end1.empty()) {
            shortest_path1 = std::move(shortest_path(graph1, shortest_path_start1, shortest_path_end1));
            if (debug) {
                std::cerr << "got graph 1 path of length " << shortest_path1.size() << "\n";
//                for (auto n : shortest_path1) {
//                    std::cerr << ' ' << n;
//                }
//                std::cerr << '\n';
            }
        }
        if (!shortest_path1.empty()) {
            // try a shortest path from the end of the alignment on graph1
            std::vector<uint64_t> shortest_path_start2, shortest_path_end2;
            if (aln_fwd.empty()) {
                shortest_path_start2 = sources2;
            }
            else {
                shortest_path_start2.push_back(aln_fwd.back().node_id2);
            }
            if (aln_rev.empty()) {
                shortest_path_end2 = sinks2;
            }
            else {
                shortest_path_end2.push_back(aln_rev.front().node_id2);
            }
            
            if (debug) {
                std::cerr << "checking for simple path from graph 2 starts:\n";
                for (auto n : shortest_path_start2) {
                    std::cerr << '\t' << n << '\n';
                }
                std::cerr << "and graph 2 ends:\n";
                for (auto n : shortest_path_end2) {
                    std::cerr << '\t' << n << '\n';
                }
            }
            
            if (!shortest_path_start2.empty() && !shortest_path_end2.empty()) {
                shortest_path2 = std::move(shortest_path(graph2, shortest_path_start2, shortest_path_end2));
            }
            
            if (!shortest_path2.empty()) {
                if (debug) {
                    std::cerr << "got graph 2 path of length " << shortest_path2.size() << "\n";
//                    for (auto n : shortest_path2) {
//                        std::cerr << ' ' << n;
//                    }
//                    std::cerr << '\n';
                }
                
                found_path = true;
                if (!aln_fwd.empty()) {
                    shortest_path1.erase(shortest_path1.begin());
                    shortest_path2.erase(shortest_path2.begin());
                }
                if (!aln_rev.empty()) {
                    shortest_path1.pop_back();
                    shortest_path2.pop_back();
                }
            }
        }
        else {
            if (debug) {
                std::cerr << "could not find graph 1 path\n";
            }
        }
    }
    
    if (!found_path) {
        
        // for reachability testing
        std::unique_ptr<SuperbubbleDistanceOracle> dist_oracle1(nullptr);
        std::unique_ptr<SuperbubbleDistanceOracle> dist_oracle2(nullptr);
        
        // the number of simple algorithmic tests we'll do before indexing for reachability
        // (this is to avoid the memory investment for indexing on simple cases)
        int unindexed_measurements_remaining = 8;
        
        auto test_reachability = [&](size_t trim_left, size_t trim_right) -> bool {
            // choose which nodes' reachability we'll be determining
            
            bool allow_equal = false; // if we're actually at before-the-first or past-the-last
            std::vector<AlignedPair> left_ends;
            if (trim_left == aln_fwd.size()) {
                for (auto node_id1 : sources1) {
                    for (auto node_id2 : sources2) {
                        left_ends.emplace_back(node_id1, node_id2);
                    }
                }
                allow_equal = true;
            }
            else {
                left_ends.emplace_back(aln_fwd[aln_fwd.size() - 1 - trim_left]);
            }
            
            std::vector<AlignedPair> right_ends;
            if (trim_right == aln_rev.size()) {
                for (auto node_id1 : sinks1) {
                    for (auto node_id2 : sinks2) {
                        right_ends.emplace_back(node_id1, node_id2);
                    }
                }
                allow_equal = true;
            }
            else {
                right_ends.emplace_back(aln_rev[trim_right]);
            }
            
            for (const auto& aln_pair_left : left_ends) {
                for (const auto& aln_pair_right : right_ends) {
                    
                    if (debug) {
                        std::cerr << "checking reachability between graph 1: " << aln_pair_left.node_id1 << " -> " <<  aln_pair_right.node_id1;
                        if (dist_oracle1.get()) {
                            std::cerr << " (" << (int64_t) dist_oracle1->min_distance(aln_pair_left.node_id1, aln_pair_right.node_id1) << ") and graph 2: " << aln_pair_left.node_id2 << " -> " <<  aln_pair_right.node_id2 << " (" << (int64_t) dist_oracle2->min_distance(aln_pair_left.node_id2, aln_pair_right.node_id2) << ")";
                        }
                        std::cerr << '\n';
                    }
                    
                    if (!allow_equal &&
                        (aln_pair_left.node_id1 == aln_pair_right.node_id1 ||
                         aln_pair_left.node_id2 == aln_pair_right.node_id2)) {
                        continue;
                    }
                    
                    if (unindexed_measurements_remaining > 0) {
                        --unindexed_measurements_remaining;
                        // TODO: use "is_reachable" instead? doesn't currently have a multi-source multi-target implementation
                        if (!shortest_path(graph1, aln_pair_left.node_id1, aln_pair_right.node_id1).empty() &&
                            !shortest_path(graph2, aln_pair_left.node_id2, aln_pair_right.node_id2).empty()) {
                            return true;
                        }
                    }
                    else {
                        if (dist_oracle1.get() == nullptr) {
                            if (debug) {
                                std::cerr << "switching to indexed reachability\n";
                            }
                            dist_oracle1.reset(new SuperbubbleDistanceOracle(graph1));
                            dist_oracle2.reset(new SuperbubbleDistanceOracle(graph2));
                        }
                        if (dist_oracle1->min_distance(aln_pair_left.node_id1, aln_pair_right.node_id1) != -1 &&
                            dist_oracle2->min_distance(aln_pair_left.node_id2, aln_pair_right.node_id2) != -1) {
                            // this combo can reach each other
                            return true;
                        }
                    }
                }
            }
            return false;
        };
        
        // bisect search to find the longest portion of the paths that we can include
        int64_t lo = 1; // we already checked 0 with the special case
        int64_t hi = aln_fwd.size() + aln_rev.size();
        while (lo <= hi) {
            int64_t total_trim = (lo + hi) / 2;
            bool success = false;
            
            if (debug) {
                std::cerr << "attempting trim total " << total_trim << " from range " << lo << ", " << hi << '\n';
            }
            size_t l_min = std::max<int64_t>(0, total_trim - aln_rev.size());
            size_t l_max = std::min<size_t>(total_trim, aln_fwd.size());
            for (size_t l = l_min; l <= l_max; ++l) {
                if (test_reachability(l, total_trim -  l)) {
                    // the aligments can reach each other if we trim this amount
                    left_trim = l;
                    right_trim = total_trim - l;
                    success = true;
                    if (debug) {
                        std::cerr << "found a reachable trim combo of " << left_trim << ", " << right_trim << '\n';
                    }
                    break;
                }
            }
            
            if (success) {
                // try to trim less
                hi = total_trim - 1;
            }
            else {
                // we'll have to trim more
                lo = total_trim + 1;
            }
        }
        
        // find the shortest paths between the reachable portions, which we'll use to make a
        // double deletion
        std::vector<uint64_t> search_sources1, search_sources2, search_sinks1, search_sinks2;
        if (left_trim == aln_fwd.size()) {
            search_sources1 = sources1;
            search_sources2 = sources2;
        }
        else {
            auto aln_pair = aln_fwd[aln_fwd.size() - left_trim - 1];
            search_sources1.push_back(aln_pair.node_id1);
            search_sources2.push_back(aln_pair.node_id2);
        }
        if (right_trim == aln_rev.size()) {
            search_sinks1 = sinks1;
            search_sinks2 = sinks2;
        }
        else {
            auto aln_pair = aln_rev[right_trim];
            search_sinks1.push_back(aln_pair.node_id1);
            search_sinks2.push_back(aln_pair.node_id2);
        }
        
        shortest_path1 = std::move(shortest_path(graph1, search_sources1, search_sinks1));
        shortest_path2 = std::move(shortest_path(graph2, search_sources2, search_sinks2));
        
        // if we used the graph sources/sinks in the shortest path, we should include them
        // in the alignment, but not if we used the final position of the fwd/rev alignment
        // TODO: ugly
        if (left_trim != aln_fwd.size()) {
            shortest_path1.erase(shortest_path1.begin());
            shortest_path2.erase(shortest_path2.begin());
        }
        if (right_trim != aln_rev.size()) {
            shortest_path1.pop_back();
            shortest_path2.pop_back();
        }
    }
    
    Alignment alignment;
    alignment.reserve(aln_fwd.size() + shortest_path1.size() + shortest_path2.size() + aln_rev.size() - left_trim - right_trim);
    
    // add the untrimmed part of the forward alignment
    for (size_t i = 0, end = aln_fwd.size() - left_trim; i < end; ++i) {
        alignment.emplace_back(aln_fwd[i]);
    }
    
    // add the shortest paths as a double deletion
    for (size_t i = 0; i < shortest_path1.size(); ++i) {
        alignment.emplace_back(shortest_path1[i], AlignedPair::gap);
    }
    for (size_t i = 0; i < shortest_path2.size(); ++i) {
        alignment.emplace_back(AlignedPair::gap, shortest_path2[i]);
    }
    
    // add the untrimmed part of the reverse alignment
    for (size_t i = right_trim; i < aln_rev.size(); ++i) {
        alignment.emplace_back(aln_rev[i]);
    }
    
    if (score_out) {
        *score_out = rescore(alignment, graph1, graph2, params, false);
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
template<bool Forward, class BackingMap, int NumPW, class Graph, class PruneFunc, class UpdateFunc,
         class NextFunc1, class NextFunc2, class StopFunc, class GreedyFunc>
inline std::pair<uint64_t, uint64_t>
wfa_iteration(std::deque<std::queue<std::tuple<uint64_t, uint64_t, int, uint64_t, uint64_t, int>>>& queue,
              int64_t& queue_min_score, BackingMap& backpointer,
              const Graph& graph1, const Graph& graph2, const AlignmentParameters<NumPW>& wfa_params,
              const PruneFunc& prune_function, const UpdateFunc& update_function,
              const NextFunc1& next_function1, const NextFunc2& next_function2,
              const StopFunc& stop_function, const GreedyFunc& greedy_function) {
    
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
    
    
    if (Forward) {
        if (here_comp == 0) {
            // match/mismatch
            
            if (greedy_function(here_id1, here_id2)) {
                // we don't need to branch out into scored edits here
                if (debug) {
                    std::cerr << "doing greedy extension\n";
                }
                enqueue(here_id1, here_id2, here_comp, next_function1(here_id1).front(), next_function2(here_id2).front(), 0, 0);
            }
            else {
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
                // extend deletion
                for (auto next_id2 : next_function2(here_id2)) {
                    enqueue(here_id1, here_id2, here_comp, here_id1, next_id2, here_comp,
                            wfa_params.gap_extend[-here_comp - 1]);
                }
            }
        }
    }
    else {
        // reverse iteration
        if (here_comp == 0) {
            
            if (here_id1 < graph1.node_size() && here_id2 < graph2.node_size()) {
                // match/mismatch
                uint64_t penalty = graph1.label(here_id1) == graph2.label(here_id2) ? 0 : wfa_params.mismatch;
                for (auto next_id1 : next_function1(here_id1)) {
                    for (auto next_id2 : next_function2(here_id2)) {
                        enqueue(here_id1, here_id2, here_comp, next_id1, next_id2, 0, penalty);
                    }
                }
            }
            // gap close
            for (int i = 0; i < NumPW; ++i) {
                enqueue(here_id1, here_id2, here_comp, here_id1, here_id2, i + 1, 0);
                enqueue(here_id1, here_id2, here_comp, here_id1, here_id2, -i - 1, 0);
            }
            
        }
        else if (here_comp > 0) {
            if (here_id1 < graph1.node_size()) {
                for (auto next_id1 : next_function1(here_id1)) {
                    // extend insertion
                    enqueue(here_id1, here_id2, here_comp, next_id1, here_id2, here_comp,
                            wfa_params.gap_extend[here_comp - 1]);
                    // open insertion
                    enqueue(here_id1, here_id2, here_comp, next_id1, here_id2, 0,
                            wfa_params.gap_open[here_comp - 1] + wfa_params.gap_extend[here_comp - 1]);
                }
            }
        }
        else {
            if (here_id2 < graph2.node_size()) {
                for (auto next_id2 : next_function2(here_id2)) {
                    // extend deletion
                    enqueue(here_id1, here_id2, here_comp, here_id1, next_id2, here_comp,
                            wfa_params.gap_extend[-here_comp - 1]);
                    // open deletion
                    enqueue(here_id1, here_id2, here_comp, here_id1, next_id2, 0,
                            wfa_params.gap_open[-here_comp - 1] + wfa_params.gap_extend[-here_comp - 1]);
                }
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

template<class BackingMap, class Graph>
Alignment wfa_traceback_rev(BackingMap& backpointer, uint64_t tb_node1, uint64_t tb_node2,
                            const Graph& graph1, const Graph& graph2) {
        
    Alignment alignment;
    
    // traceback using backpointers
    int tb_comp = 0;
    uint64_t next_id1, next_id2;
    int next_comp;
    std::tie(next_id1, next_id2, next_comp) = backpointer[std::make_tuple(tb_node1, tb_node2, tb_comp)];
    while (next_id1 != -1 && next_id2 != -1) {

        if (next_id1 != tb_node1 && next_id2 != tb_node2) {
            alignment.emplace_back(next_id1, next_id2);
        }
        else if (next_id1 != tb_node1) {
            alignment.emplace_back(next_id1, AlignedPair::gap);
        }
        else if (next_id2 != tb_node2) {
            alignment.emplace_back(AlignedPair::gap, next_id2);
        }
        
        tb_node1 = next_id1;
        tb_node2 = next_id2;
        tb_comp = next_comp;
        std::tie(next_id1, next_id2, next_comp) = backpointer[std::make_tuple(tb_node1, tb_node2, tb_comp)];
    }
    
    // note: don't reverse the traceback since it was already reversed
    
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
    queue.back().emplace(-1, -1, 0, graph1.node_size(), graph2.node_size(), 0);
    
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
    
    auto greedy = [&](uint64_t node_id1, uint64_t node_id2) -> bool {
        if (get_next1(node_id1).size() == 1 && get_next2(node_id2).size() == 1 && !sink_set1.count(node_id1) && !sink_set2.count(node_id2)) {
            return (graph1.label(get_next1(node_id1).front()) == graph2.label(get_next2(node_id2).front()));
        }
        return false;
    };
    
    std::pair<uint64_t, uint64_t> end(-1, -1);
    while (end == std::pair<uint64_t, uint64_t>(-1, -1)) {
        end = wfa_iteration<true>(queue, queue_min_score, backpointer, graph1, graph2, wfa_params,
                                  prune_function, update_function, get_next1, get_next2, stop, greedy);
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
    
    static const bool debug = false;
    
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
    queue_fwd.back().emplace(-1, -1, 0, short_graph.node_size(), long_graph.node_size(), 0);
    queue_rev.emplace_back();
    for (auto short_sink_id : sinks_short) {
        for (auto long_sink_id : sinks_long) {
            queue_rev.back().emplace(-1, -1, 0, short_sink_id, long_sink_id, 0);
        }
    }
    
    // for O(1) membership queries
    std::unordered_set<uint64_t> source_set_short(sources_short.begin(), sources_short.end());
    std::unordered_set<uint64_t> source_set_long(sources_long.begin(), sources_long.end());
    
    // transition functions for each direction
    auto get_next_short = [&](uint64_t node_id) -> const std::vector<uint64_t>& {
        return node_id == short_graph.node_size() ? sources_short : short_graph.next(node_id);
    };
    auto get_next_long = [&](uint64_t node_id) -> const std::vector<uint64_t>& {
        return node_id == long_graph.node_size() ? sources_long : long_graph.next(node_id);
    };
    // TODO: does it make more sense to do the edge condition checks in wfa_iteration here?
    auto get_prev_short = [&](uint64_t node_id) -> const std::vector<uint64_t> {
        auto prev = short_graph.previous(node_id);
        if (source_set_short.count(node_id)) {
            prev.push_back(short_graph.node_size());
        }
        return prev;
    };
    auto get_prev_long = [&](uint64_t node_id) -> const std::vector<uint64_t> {
        auto prev = long_graph.previous(node_id);
        if (source_set_long.count(node_id)) {
            prev.push_back(long_graph.node_size());
        }
        return prev;
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
        if (stop_score == std::numeric_limits<int64_t>::max()) {
            auto it = rev_score.find(std::get<0>(pos));
            if (it != rev_score.end()) {
                // the wavefronts have met on the short graph
                for (const auto& rev_pos : it->second) {
                    if (std::get<1>(pos) == rev_pos.first || // this covers the case that they're both in the boundary
                        (std::get<1>(pos) != long_graph.node_size() && rev_pos.first != long_graph.node_size() &&
                         dist_oracle.min_distance(std::get<1>(pos), rev_pos.first) != -1)) {
                        // the wavefronts are also reachable on the long graph
                        if (debug) {
                            std::cerr << "found fwd iteration meet on node " << std::get<0>(pos) << " with score " << s << '\n';
                        }
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
        if (stop_score == std::numeric_limits<int64_t>::max()) {
            auto it = fwd_score.find(std::get<0>(pos));
            if (it != fwd_score.end()) {
                // the wavefronts have met on the short graph
                for (const auto& fwd_pos : it->second) {
                    if (std::get<1>(pos) == fwd_pos.first || // this covers the case that they're both in the boundary
                        (std::get<1>(pos) != long_graph.node_size() && fwd_pos.first != long_graph.node_size() &&
                         dist_oracle.min_distance(fwd_pos.first, std::get<1>(pos)) != -1)) {
                        // the wavefronts are also reachable on the long graph
                        if (debug) {
                            std::cerr << "found rev iteration meet on node " << std::get<0>(pos) << " with score " << s << '\n';
                        }
                        stop_score = s + scope;
                    }
                }
            }
        }
    };

    auto stop = [&](uint64_t node_id1, uint64_t node_id2, int comp) {
        return queue_min_score_fwd >= stop_score && queue_min_score_rev >= stop_score;
    };
    
    auto no_greedy = [](uint64_t node_id1, uint64_t node_id2) -> bool {
        return false;
    };
    
    // do the core WFA iterations
    std::pair<uint64_t, uint64_t> end_fwd(-1, -1), end_rev(-1, -1);
    while (end_fwd == std::pair<uint64_t, uint64_t>(-1, -1) &&
           end_rev == std::pair<uint64_t, uint64_t>(-1, -1)) {
        if (queue_min_score_fwd <= queue_min_score_rev) {
            if (debug) {
                std::cerr << "do forward iteration\n";
            }
            end_fwd = wfa_iteration<true>(queue_fwd, queue_min_score_fwd, backpointer_fwd,
                                          short_graph, long_graph, wfa_params,
                                          no_prune, update_fwd, get_next_short, get_next_long, stop, no_greedy);
        }
        else {
            if (debug) {
                std::cerr << "do reverse iteration\n";
            }
            end_rev = wfa_iteration<false>(queue_rev, queue_min_score_rev, backpointer_rev,
                                           short_graph, long_graph, wfa_params,
                                           no_prune, update_rev, get_prev_short, get_prev_long, stop, no_greedy);
        }
    }
    
    if (debug) {
        std::cerr << "ending WFA iterations with queued scores " << queue_min_score_fwd << " and " << queue_min_score_rev << '\n';
    }
    
    int64_t opt_score = std::numeric_limits<int64_t>::max();
    uint64_t opt_short_id = -1, opt_long_fwd_id = -1, opt_long_rev_id = -1;
    
    // TODO: should i make it possible to transition into the deletion extend component after the big gap?
    
    // find the pair of locations with the best score
    for (const auto& fwd_rec : fwd_score) {
        auto it = rev_score.find(fwd_rec.first);
        if (it != rev_score.end()) {
            // TODO: can i avoid the N^2 all pairs behavior?
            for (const auto& fwd_pos : fwd_rec.second) {
                if (fwd_pos.first == long_graph.node_size()) {
                    continue;
                }
                for (const auto& rev_pos : it->second) {
                    if (rev_pos.first == long_graph.node_size()) {
                        continue;
                    }
                    size_t dist = dist_oracle.min_distance(fwd_pos.first, rev_pos.first);
                    if (dist != -1) {
                        // the positions are reachable in the long graph
                        int64_t score = wfa_params.gap_open[0] + wfa_params.gap_extend[0] * dist;
                        for (int i = 1; i < NumPW; ++i) {
                            score = std::min<int64_t>(score, wfa_params.gap_open[i] + wfa_params.gap_extend[i] * dist);
                        }
                        score += fwd_pos.second + rev_pos.second;
                        if (debug) {
                            std::cerr << "checking meet combo on short graph node " << fwd_rec.first << " between long graph nodes " << fwd_pos.first << " and " << rev_pos.first << " separated by " << dist << " with score " << score << "\n";
                        }
                        if (score < opt_score) {
                            // this pair is the new opt
                            opt_score = score;
                            opt_short_id = fwd_rec.first;
                            opt_long_fwd_id = fwd_pos.first;
                            opt_long_rev_id = rev_pos.first;
                            if (debug) {
                                std::cerr << "this is the new opt\n";
                            }
                        }
                    }
                }
            }
        }
    }
    
    if (debug) {
        std::cerr << "final opt is short graph node " << opt_short_id << " and long graph path between " << opt_long_fwd_id << " and " << opt_long_rev_id << '\n';
    }
    
    Alignment fwd_traceback = wfa_traceback(backpointer_fwd, opt_short_id, opt_long_fwd_id,
                                            short_graph, long_graph);
    Alignment rev_traceback = wfa_traceback_rev(backpointer_rev, opt_short_id, opt_long_rev_id,
                                                short_graph, long_graph);
    // find the path that allows the shortest gap between the forward and reverse alignments
    auto path_between = shortest_path(long_graph, opt_long_fwd_id, opt_long_rev_id);
    if (debug) {
        std::cerr << "forward traceback:\n";
        for (auto aln_pair : fwd_traceback) {
            std::cerr << '\t' << (int64_t) aln_pair.node_id1 << '\t' << (int64_t) aln_pair.node_id2 << '\n';
        }
        std::cerr << "middle path\n";
        for (size_t i = 0; i < path_between.size(); ++i) {
            std::cerr << '\t' << -1 << '\t' << path_between[i] << '\n';
        }
        std::cerr << "reverse traceback:\n";
        for (auto aln_pair : ReverseForEachAdapter<Alignment>(rev_traceback)) {
            std::cerr << '\t' << (int64_t) aln_pair.node_id1 << '\t' << (int64_t) aln_pair.node_id2 << '\n';
        }
    }
    for (size_t i = 1; i < path_between.size(); ++i) {
        fwd_traceback.emplace_back(AlignedPair::gap, path_between[i]);
    }
    for (auto& aln_pair : rev_traceback) {
        fwd_traceback.push_back(aln_pair);
    }
    
    if (score_out) {
        *score_out = convert_wfa_score(fwd_traceback, opt_score, params, factor);
    }
    
    if (debug) {
        std::cerr << "combined alignment:\n";
        for (const auto& aln_pair : fwd_traceback) {
            std::cerr << '\t' << (int64_t) aln_pair.node_id1 << '\t' << aln_pair.node_id2 << '\n';
        }
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
                                             up.I[pw] - params.gap_extend[pw]);
                cell.D[pw] = std::max<IntDP>(left.M - params.gap_open[pw] - params.gap_extend[pw],
                                             left.D[pw] - params.gap_extend[pw]);
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


template<int NumPW>
int64_t rescore(const Alignment& aln, const BaseGraph& graph1, const BaseGraph& graph2,
                const AlignmentParameters<NumPW>& params, bool wfa_style) {
    
    auto local_params = params;
    
    if (wfa_style) {
        local_params.match = 0;
        local_params.mismatch = 2 * (params.match + params.mismatch);
        local_params.gap_open[0] = 2 * params.gap_open[0];
        local_params.gap_extend[0] = 2 * params.gap_extend[0] + params.match;
    }
    
    int64_t score = 0;
    for (size_t i = 0; i < aln.size(); ++i) {
        if (aln[i].node_id1 != AlignedPair::gap && aln[i].node_id2 != AlignedPair::gap) {
            if (graph1.label(aln[i].node_id1) == graph2.label(aln[i].node_id2)) {
                score += local_params.match;
            }
            else {
                score -= local_params.mismatch;
            }
        }
    }
    
    for (size_t i = 0; i < aln.size();) {
        if (aln[i].node_id1 == AlignedPair::gap) {
            size_t j = i + 1;
            while (j < aln.size() && aln[j].node_id1 == AlignedPair::gap) {
                ++j;
            }
            score -= local_params.gap_open[0] + (j - i) * local_params.gap_extend[0];
            i = j;
        }
        else {
            ++i;
        }
    }
    
    for (size_t i = 0; i < aln.size();) {
        if (aln[i].node_id2 == AlignedPair::gap) {
            size_t j = i + 1;
            while (j < aln.size() && aln[j].node_id2 == AlignedPair::gap) {
                ++j;
            }
            score -= local_params.gap_open[0] + (j - i) * local_params.gap_extend[0];
            i = j;
        }
        else {
            ++i;
        }
    }
    
    return score;
}

// Broken!! The search limit doesn't guarantee that *earlier* points in the sparse data structure also obey that search limit
/*
template<class StringLike>
std::tuple<Alignment, std::pair<size_t, size_t>, std::pair<size_t, size_t>>
longest_common_subsequence_nonrepeating(const StringLike& str1, const StringLike& str2) {
    
    static const bool debug = true;

    if (debug) {
        std::cerr << "starting longest common nonrepeating subsequence algorithm\n";
        std::cerr << "str1:\n";
        for (auto c : str1) {
            std::cerr << c << ' ';
        }
        std::cerr << '\n';
        std::cerr << "str2:\n";
        for (auto c : str2) {
            std::cerr << c << ' ';
        }
        std::cerr << '\n';
    }
    
    // find the previous occurrences and the active indexes
    std::vector<std::vector<size_t>> active_indexes(str1.size());
    std::vector<size_t> search_limit1(str1.size());
    std::vector<size_t> search_limit2(str2.size());
    {
        std::unordered_map<typename StringLike::value_type, std::vector<size_t>> occurrences1;
        occurrences1.reserve(str1.size());
        
        size_t max_lim = 0;
        for (size_t i = 0; i < str1.size(); ++i) {
            auto it = occurrences1.find(str1[i]);
            if (it != occurrences1.end()) {
                max_lim = std::max(max_lim, it->second.back() + 1);
                it->second.push_back(i);
            }
            else {
                occurrences1[str1[i]].push_back(i);
            }
            search_limit1[i] = max_lim;
        }

        max_lim = 0;
        std::unordered_map<typename StringLike::value_type, size_t> prev2;
        for (size_t j = 0; j < str2.size(); ++j) {
            auto it = prev2.find(str2[j]);
            if (it != prev2.end()) {
                max_lim = std::max(max_lim, it->second + 1);
            }
            prev2[str2[j]] = j;
            auto it2 = occurrences1.find(str2[j]);
            if (it2 != occurrences1.end()) {
                for (size_t i : it2->second) {
                    active_indexes[i].push_back(j);
                }
            }
            search_limit2[j] = max_lim;
        }
    }
    
    if (debug) {
        std::cerr << "got search limits:\n";
        std::cerr << "str1:\n";
        for (auto p : search_limit1) {
            std::cerr << '\t' << p << '\n';
        }
        std::cerr << "str2:\n";
        for (auto p : search_limit2) {
            std::cerr << '\t' << p << '\n';
        }
        std::cerr << "active indexes\n";
        for (size_t i = 0; i < active_indexes.size(); ++i) {
            for (auto j : active_indexes[i]) {
                std::cerr << '\t' << i << '\t' << j << '\n';
            }
        }
    }
    
    if (debug) {
        std::cerr << "initializing RMQ\n";
    }
    
    // initialize the RMQ
    std::vector<std::tuple<size_t, size_t, int64_t>> tree_data;
    for (size_t i = 0; i < active_indexes.size(); ++i) {
        for (size_t j : active_indexes[i]) {
            tree_data.emplace_back(i, j, -1);
        }
    }
    OrthogonalMaxSearchTree<size_t, size_t, int64_t> orthogonal_rmq(tree_data);
    {
        // we no longer need this
        auto dummy = std::move(tree_data);
    }
    
    if (debug) {
        std::cerr << "beginning sparse DP\n";
    }
        
    std::unordered_map<std::pair<size_t, size_t>, std::pair<size_t, size_t>> backpointer;
    int64_t opt = -1;
    size_t opt_i = -1;
    size_t opt_j = -1;
    for (size_t i = 0; i < active_indexes.size(); ++i) {
        for (size_t j : active_indexes[i]) {
            
            if (debug) {
                std::cerr << "active index " << i << ", " << j << " with limits " << search_limit1[i] << ", " << search_limit2[j] << '\n';
            };
            
            auto it_max = orthogonal_rmq.range_max(search_limit1[i], i, search_limit2[j], j);
            auto it = orthogonal_rmq.find(i, j);
            if (it_max == orthogonal_rmq.end() || std::get<2>(*it_max) == -1) {
                // this is the start of a new common subsequence
                orthogonal_rmq.update(it, 1);
                if (debug) {
                    std::cerr << "no extension possible, starting new subsequence\n";
                }
                if (opt == -1) {
                    if (debug) {
                        std::cerr << "is new opt\n";
                    }
                    opt = 1;
                    opt_i = i;
                    opt_j = j;
                }
            }
            else {
                // we can extend from a previous common subsequence
                orthogonal_rmq.update(it, std::get<2>(*it_max) + 1);
                backpointer[std::make_pair(i, j)] = std::make_pair(std::get<0>(*it_max), std::get<1>(*it_max));
                if (debug) {
                    std::cerr << "got extension from " << std::get<0>(*it_max) << ", " << std::get<1>(*it_max) << " length " << std::get<2>(*it_max) << "\n";
                }
                if (std::get<2>(*it) > opt) {
                    if (debug) {
                        std::cerr << "is new opt\n";
                    }
                    opt = std::get<2>(*it);
                    opt_i = i;
                    opt_j = j;
                }
            }
        }
    }
    
    if (debug) {
        std::cerr << "tracing back LCS\n";
    }
    
    std::tuple<Alignment, std::pair<size_t, size_t>, std::pair<size_t, size_t>> lcs;
    
    if (opt_i != -1) {
        // traceback to find alignment
        auto& aln = std::get<0>(lcs);
        size_t tb_i = opt_i, tb_j = opt_j;
        std::get<2>(lcs) = std::make_pair(tb_i, tb_j);
        aln.emplace_back(str1[tb_i], str2[tb_j]);
        while (backpointer.count(std::make_pair(tb_i, tb_j))) {
            size_t next_i, next_j;
            std::tie(next_i, next_j) = backpointer[std::make_pair(tb_i, tb_j)];
            for (size_t i = tb_i - 1; i > next_i; --i) {
                aln.emplace_back(str1[i], AlignedPair::gap);
            }
            for (size_t j = tb_j - 1; j > next_j; --j) {
                aln.emplace_back(AlignedPair::gap, str2[j]);
            }
        }
        std::get<1>(lcs) = std::make_pair(tb_i, tb_j);
        std::reverse(aln.begin(), aln.end());
    }
    else {
        // there was nothing in common
        std::get<1>(lcs) = std::pair<size_t, size_t>(-1, -1);
        std::get<2>(lcs) = std::pair<size_t, size_t>(-1, -1);
    }
    
    return lcs;
}
*/

template<class StringLike>
Alignment long_common_subsequence_nonrepeating(const StringLike& str1, const StringLike& str2) {
    
    // get the LCS without any non-repeating constraint
    auto lcs_aln = align_hs(str1, str2);
    
    // make prefix sum vector over number of matched positions
    std::vector<uint64_t> matched_prefix_sum(lcs_aln.size() + 1, 0);
    for (size_t i = 0; i < lcs_aln.size(); ++i) {
        matched_prefix_sum[i + 1] = matched_prefix_sum[i] + uint64_t(lcs_aln[i].node_id1 != AlignedPair::gap && lcs_aln[i].node_id2 != AlignedPair::gap);
    }
    
    // make an inverse map from string characters to their location in the alignment
    std::vector<size_t> aln_idx1(str1.size()), aln_idx2(str2.size());
    for (size_t i = 0, idx1 = 0, idx2 = 0; i < lcs_aln.size(); ++i) {
        const auto& ap = lcs_aln[i];
        if (ap.node_id1 != AlignedPair::gap) {
            aln_idx1[idx1++] = i;
        }
        if (ap.node_id2 != AlignedPair::gap) {
            aln_idx2[idx2++] = i;
        }
    }
    
    // figure out how far back you can look from each position without a repeat
    std::vector<size_t> search_limit1(str1.size()), search_limit2(str2.size());
    for (bool do_str1 : {true, false}) {
        const auto& str = do_str1 ? str1 : str2;
        auto& search_limit = do_str1 ? search_limit1 : search_limit2;
        
        std::unordered_map<typename StringLike::value_type, size_t> prev;
        prev.reserve(str.size());
        
        size_t max_lim = 0;
        for (size_t i = 0; i < str.size(); ++i) {
            auto it = prev.find(str1[i]);
            if (it != prev.end()) {
                max_lim = std::max(max_lim, it->second + 1);
                it->second = i;
            }
            else {
                prev[str[i]] = i;
            }
            search_limit[i] = max_lim;
        }
    }
    
    // get the non-repeating interval with the most matches
    size_t opt_begin = 0;
    size_t opt_end = 1;
    for (size_t i = 1; i < lcs_aln.size(); ++i) {
        size_t begin = std::max(aln_idx1[search_limit1[i]], aln_idx2[search_limit2[i]]);
        if (matched_prefix_sum[i + 1] - matched_prefix_sum[begin] > matched_prefix_sum[opt_end] - matched_prefix_sum[opt_begin]) {
            opt_begin = begin;
            opt_end = i + 1;
        }
    }
    
    // trim gaps off the end
    while (opt_begin < opt_end &&
           (lcs_aln[opt_begin].node_id1 == AlignedPair::gap || lcs_aln[opt_begin].node_id2 == AlignedPair::gap)) {
        ++opt_begin;
    }
    while (opt_end > opt_begin &&
           (lcs_aln[opt_end - 1].node_id1 == AlignedPair::gap || lcs_aln[opt_end - 1].node_id2 == AlignedPair::gap)) {
        --opt_end;
    }
    for (size_t i = 0, j = opt_begin; j < opt_end; ++i, ++j) {
        lcs_aln[i] = lcs_aln[j];
    }
    lcs_aln.resize(opt_end - opt_begin);
    
    return lcs_aln;
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
