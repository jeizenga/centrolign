#ifndef centrolign_stitcher_hpp
#define centrolign_stitcher_hpp

#include <vector>
#include <cstdint>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <chrono>

#include "centrolign/chain_merge.hpp"
#include "centrolign/alignment.hpp"
#include "centrolign/topological_order.hpp"
#include "centrolign/step_index.hpp"
#include "centrolign/anchorer.hpp"
#include "centrolign/subgraph_extraction.hpp"
#include "centrolign/logging.hpp"

namespace centrolign {

/*
 * Object that connects an anchor chain into a base-level alignment
 */
class Stitcher : public Extractor {
public:
    
    Stitcher();
    ~Stitcher() = default;
    
    template<class BGraph1, class BGraph2, class XMerge1, class XMerge2>
    Alignment stitch(const std::vector<std::vector<anchor_t>>& anchor_segments,
                     const BGraph1& graph1, const BGraph2& graph2,
                     const SentinelTableau& tableau1, const SentinelTableau& tableau2,
                     XMerge1& chain_merge1, XMerge2& chain_merge2) const;
    
    AlignmentParameters<3> alignment_params;
    // the maximum size matrix that will always do PO-POA in spite of overrides (sized to be ~1x1 monomer)
    size_t max_trivial_size = 30000;
    // the minimum matrix size that will trigger WFA instead of PO-POA
    size_t min_wfa_size = 10000000;
    // the maximum size matrix that will still trigger WFA
    size_t max_wfa_size = 50000000;
    // only use WFA if the distances suggest their could be a gapless alignment, up to this ratio
    double max_wfa_ratio = 1.05;
    // prune WFA positions that are this far behind the opt
    size_t wfa_pruning_dist = 25;
    // use approximate deletion alignment when one graph is this much larger than the other
    size_t deletion_alignment_ratio = 4;
    // the longest the shorter of the two graphs can be
    size_t deletion_alignment_short_max_size = 4000;
    // the shortest the longer of the two graph can be
    size_t deletion_alignment_long_min_size = 2000;
    // if true, destroy the XMerge data structures after extracting stitch graphs
    bool clean_merge_structs = false;
    
private:
    
    static const bool debug;
    static const bool instrument;
    
    void subalign(const SubGraphInfo& extraction1, const SubGraphInfo& extraction2,
                  Alignment& stitched, bool only_deletion_alns) const;
    
    template<int NumPW>
    Alignment do_alignment(const SubGraphInfo& extraction1, const SubGraphInfo& extraction2,
                           bool only_deletion_alns, const AlignmentParameters<NumPW>& params) const;

    
    void do_instrument(const SubGraphInfo& extraction1, const SubGraphInfo& extraction2,
                       int64_t score1, int64_t score2, double dur1, double dur2) const;
    
    template<class BGraph>
    void log_subpath_info(const BGraph& graph1, const BGraph& graph2,
                          const std::vector<std::pair<SubGraphInfo, SubGraphInfo>>& extractions) const;
    
};

/*
 * Template implementations
 */


template<class BGraph1, class BGraph2, class XMerge1, class XMerge2>
Alignment Stitcher::stitch(const std::vector<std::vector<anchor_t>>& anchor_segments,
                           const BGraph1& graph1, const BGraph2& graph2,
                           const SentinelTableau& tableau1, const SentinelTableau& tableau2,
                           XMerge1& xmerge1, XMerge2& xmerge2) const {

    size_t next_log_idx = 0;
    std::vector<size_t> logging_indexes;
    size_t size = 0;
    if (logging::level >= logging::Debug) {
        size_t total_anchor_len = 0;
        for (const auto& segment : anchor_segments) {
            size += segment.size() - 1;
            for (const auto& anchor : segment) {
                total_anchor_len += anchor.walk1.size();
            }
        }
        logging_indexes = get_logging_indexes(size);
        
        logging::log(logging::Debug, "Stitching alignment between " + std::to_string(anchor_segments.size()) + " anchor segments containing " + std::to_string(size + anchor_segments.size()) + " anchors with total length " + std::to_string(total_anchor_len));
    }
    
    // get the subgraphs
    std::vector<std::vector<std::pair<SubGraphInfo, SubGraphInfo>>> within_segment_graphs;
    std::vector<std::pair<SubGraphInfo, SubGraphInfo>> between_segment_graphs;
    std::tie(within_segment_graphs, between_segment_graphs) = extract_graphs_between(anchor_segments,
                                                                                     graph1, graph2,
                                                                                     tableau1, tableau2,
                                                                                     xmerge1, xmerge2);
    
    assert(within_segment_graphs.size() + 1 == between_segment_graphs.size());
    
    if (clean_merge_structs) {
        // this is the last time we need these structs, so we can clear them out here
        {
            XMerge1 dummy = std::move(xmerge1);
        }
        {
            XMerge2 dummy = std::move(xmerge2);
        }
    }
    
    // TODO: break out cases for between/within segments
    
    if (instrument) {
        // record the locations of the subgraphs
        std::cerr << "between segment graphs:\n";
        log_subpath_info(graph1, graph2, between_segment_graphs);
    }
    
    Alignment stitched;
    
    size_t idx = 0;
    for (size_t i = 0; i < between_segment_graphs.size(); ++i) {
        
        if (i != 0) {
            // add the within-segment alignments
            const auto& segment_graphs = within_segment_graphs[i - 1];
            const auto& segment = anchor_segments[i - 1];
            
            assert(segment_graphs.size() + 1 == segment.size());
            if (instrument) {
                // record the locations of the subgraphs
                std::cerr << "segment " << (i - 1) << " graphs:\n";
                log_subpath_info(graph1, graph2, segment_graphs);
            }
            
            for (size_t j = 0; j < segment.size(); ++j) {
                
                if (j != 0) {
                    // align graph between anchors, possibly with expensive alignments
                    
                    if (instrument) {
                        std::cerr << "subalign " << (j - 1) << " within segment " << (i - 1) << '\n';
                    }
                    const auto& stitch_pair = segment_graphs[j - 1];
                    subalign(stitch_pair.first, stitch_pair.second, stitched, false);
                    
                    if (next_log_idx < logging_indexes.size() && idx == logging_indexes[next_log_idx]) {
                        logging::log(logging::Debug, "Stitching iteration " + std::to_string(idx) + " of " + std::to_string(size));
                        ++next_log_idx;
                    }
                    ++idx;
                }
                
                // copy the anchor
                const auto& anchor = segment[j];
                for (size_t k = 0; k < anchor.walk1.size(); ++k) {
                    stitched.emplace_back(anchor.walk1[k], anchor.walk2[k]);
                }
            }
        }
        
        // align between segments, but only if it's a near perfect deletion
        if (instrument) {
            std::cerr << "subalign before segment " << i << '\n';
        }
        const auto& between_stitch_pair = between_segment_graphs[i];
        subalign(between_stitch_pair.first, between_stitch_pair.second, stitched, true);
    }
//
//    for (size_t i = 0; i < stitch_graphs.size(); ++i) {
//
//        const auto& stitch_pair = stitch_graphs[i];
//
//        if (next_log_idx < logging_indexes.size() && i == logging_indexes[next_log_idx]) {
//            logging::log(logging::Debug, "Stitching iteration " + std::to_string(i + 1) + " of " + std::to_string(size));
//            ++next_log_idx;
//        }
//
//        // make an intervening alignment
//        if (instrument) {
//            std::cerr << "subalign " << i << '\n';
//        }
//        subalign(stitch_pair.first, stitch_pair.second, stitched);
//
//        if (i < anchor_chain.size()) {
//            // copy the anchor
//            const auto& anchor = anchor_chain[i];
//            for (size_t j = 0; j < anchor.walk1.size(); ++j) {
//                stitched.emplace_back(anchor.walk1[j], anchor.walk2[j]);
//            }
//        }
//    }
    
    return stitched;
}


template<int NumPW>
Alignment Stitcher::do_alignment(const SubGraphInfo& extraction1, const SubGraphInfo& extraction2,
                                 bool only_deletion_alns, const AlignmentParameters<NumPW>& params) const {
    
    size_t mat_size = (extraction1.subgraph.node_size() + 1) * (extraction2.subgraph.node_size() + 1);
    if (instrument) {
        std::cerr << '%' << '\t' << extraction1.subgraph.node_size() << '\t' << extraction2.subgraph.node_size() << '\t' << mat_size;
        auto minmax1 = std::minmax_element(extraction1.back_translation.begin(), extraction1.back_translation.end());
        auto minmax2 = std::minmax_element(extraction2.back_translation.begin(), extraction2.back_translation.end());
        if (minmax1.first != extraction1.back_translation.end()) {
            std::cerr << '\t' << *minmax1.first << '\t' << *minmax1.second;
        }
        else {
            std::cerr << "\t-1\t-1";
        }
        if (minmax2.first != extraction2.back_translation.end()) {
            std::cerr << '\t' << *minmax2.first << '\t' << *minmax2.second;
        }
        else {
            std::cerr << "\t-1\t-1";
        }
        std::cerr << '\t';
    }

    auto begin = std::chrono::high_resolution_clock::now();
    Alignment inter_aln;
        
    if (extraction2.subgraph.node_size() == 0) {
        // align with shortest path
        inter_aln = std::move(pure_deletion_alignment(extraction1.subgraph,
                                                      extraction1.sources,
                                                      extraction1.sinks,
                                                      params));
        if (instrument) {
            std::cerr << "pd1";
        }
    }
    else if (extraction1.subgraph.node_size() == 0) {
        // align with shortest path
        inter_aln = std::move(pure_deletion_alignment(extraction2.subgraph,
                                                      extraction2.sources,
                                                      extraction2.sinks,
                                                      params));
        swap_graphs(inter_aln);
        if (instrument) {
            std::cerr << "pd2";
        }
    }
    else {
        // find the range of possible distances through this connecting graph
        size_t min1, max1, min2, max2;
        std::tie(min1, max1) = source_sink_minmax(extraction1);
        std::tie(min2, max2) = source_sink_minmax(extraction2);
        
        if (max1 * deletion_alignment_ratio <= min2 &&
            max1 <= deletion_alignment_short_max_size &&
            min2 >= deletion_alignment_long_min_size) {
            // graph1 is probably mostly a deletion of graph2
            inter_aln = std::move(deletion_wfa_po_poa(extraction1.subgraph, extraction2.subgraph,
                                                      extraction1.sources, extraction2.sources,
                                                      extraction1.sinks, extraction2.sinks, params));
            if (instrument) {
                std::cerr << "ad1";
            }
        }
        else if (max2 * deletion_alignment_ratio <= min1 &&
                 max2 <= deletion_alignment_short_max_size &&
                 min1 >= deletion_alignment_long_min_size) {
            // graph2 is probably mostly a deletion of graph1
            inter_aln = std::move(deletion_wfa_po_poa(extraction2.subgraph, extraction1.subgraph,
                                                      extraction2.sources, extraction1.sources,
                                                      extraction2.sinks, extraction1.sinks, params));
            swap_graphs(inter_aln);
            if (instrument) {
                std::cerr << "ad2";
            }
        }
        else if (mat_size <= min_wfa_size && (!only_deletion_alns || mat_size <= max_trivial_size)) {
            // the graph is small enough that we'll do full PO-POA
            inter_aln = std::move(po_poa(extraction1.subgraph, extraction2.subgraph,
                                         extraction1.sources, extraction2.sources,
                                         extraction1.sinks, extraction2.sinks, params));
            if (instrument) {
                std::cerr << "po";
            }
        }
        else if (mat_size < max_wfa_size &&
                 ((min2 * max_wfa_ratio >= min1 && min2 <= max1 * max_wfa_ratio) ||
                  (max2 * max_wfa_ratio >= min1 && max2 <= max1 * max_wfa_ratio) ||
                  (min1 * max_wfa_ratio >= min2 && min1 <= max2 * max_wfa_ratio) ||
                  (max1 * max_wfa_ratio >= min2 && max1 <= max2 * max_wfa_ratio))
                 && !only_deletion_alns) {
            // one of the endpoints is in the others (expanded) interval, which is necessary
            // and sufficient for them to overlap. attempt this aligment with WFA
            // TODO: a bail out condition?
            inter_aln = std::move(pwfa_po_poa(extraction1.subgraph, extraction2.subgraph,
                                              extraction1.sources, extraction2.sources,
                                              extraction1.sinks, extraction2.sinks, params,
                                              2 * wfa_pruning_dist));
            if (instrument) {
                std::cerr << "w";
            }
        }
        else {
            // this looks like an unalignable gap
            // TODO: is it more consistent with the confident-first alignment principal
            inter_aln = std::move(greedy_partial_alignment(extraction1.subgraph, extraction2.subgraph,
                                                           extraction1.sources, extraction2.sources,
                                                           extraction1.sinks, extraction2.sinks, params));
            if (instrument) {
                // include this as a reminder to myself that this procedure can fill them in
                uint64_t num_aligned = 0;
                for (const auto& aln_pair : inter_aln) {
                    if (aln_pair.node_id1 != AlignedPair::gap && aln_pair.node_id2 != AlignedPair::gap) {
                        ++num_aligned;
                    }
                }
                std::cerr << "u(" << num_aligned << ")";
            }
        }
    }
    
    if (instrument) {
        auto end = std::chrono::high_resolution_clock::now();
        double dur = std::chrono::duration<double, std::nano>(end - begin).count();
        std::cerr << '\t' << dur << '\n';
    }
    
    return inter_aln;
}

template<class BGraph>
void Stitcher::log_subpath_info(const BGraph& graph1, const BGraph& graph2,
                                const std::vector<std::pair<SubGraphInfo, SubGraphInfo>>& extractions) const {
        
    StepIndex steps1(graph1);
    StepIndex steps2(graph2);
    
    for (size_t i = 0; i < extractions.size(); ++i) {
        
        std::cerr << '&' << '\t' << i;
        
        for (auto val : {std::make_tuple(&steps1, &extractions[i].first, &graph1), std::make_tuple(&steps2, &extractions[i].second, &graph2)}) {
                        
            const auto& steps = *std::get<0>(val);
            const auto& extraction = *std::get<1>(val);
            const auto& graph = *std::get<2>(val);
            
            std::unordered_map<uint64_t, std::pair<size_t, size_t>> path_intervals;
            
            // scan paths by iterating in topological order
            for (auto node_id : topological_order(extraction.subgraph)) {
                for (auto step : steps.path_steps(extraction.back_translation[node_id])) {
                    
                    if (path_intervals.count(step.first)) {
                        path_intervals[step.first].second = step.second;
                    }
                    else {
                        path_intervals[step.first] = std::make_pair(step.second, step.second);
                    }
                }
            }
            
            // put them in alphabetical order
            std::vector<std::tuple<std::string, size_t, size_t>> intervals;
            for (const auto& rec : path_intervals) {
                intervals.emplace_back(graph.path_name(rec.first), rec.second.first, rec.second.second);
            }
            std::sort(intervals.begin(), intervals.end());
            
            std::cerr << '\t';
            for (size_t i = 0; i < intervals.size(); ++i) {
                if (i) {
                    std::cerr << ',';
                }
                std::cerr << std::get<0>(intervals[i]) << ':' << std::get<1>(intervals[i]) << '-' << std::get<2>(intervals[i]);
            }
        }
        std::cerr << '\n';
        
        
    }
}

}

#endif /* centrolign_stitcher_hpp */
