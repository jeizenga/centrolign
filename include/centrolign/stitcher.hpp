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
class Stitcher {
public:
    
    Stitcher();
    ~Stitcher() = default;
    
    // TODO: maybe push this to an Extractor parent class?
    template<class BGraph1, class BGraph2, class XMerge1, class XMerge2>
    std::vector<std::pair<SubGraphInfo, SubGraphInfo>>
    extract_stitch_graphs(const std::vector<anchor_t>& anchor_chain,
                          const BGraph1& graph1, const BGraph2& graph2,
                          const SentinelTableau& tableau1, const SentinelTableau& tableau2,
                          const XMerge1& chain_merge1, const XMerge2& chain_merge2) const;
    
    template<class BGraph1, class BGraph2, class XMerge1, class XMerge2>
    Alignment stitch(const std::vector<anchor_t>& anchor_chain,
                     const BGraph1& graph1, const BGraph2& graph2,
                     const SentinelTableau& tableau1, const SentinelTableau& tableau2,
                     const XMerge1& chain_merge1, const XMerge2& chain_merge2) const;
    
    // split up matches
    template<class BGraph1, class BGraph2>
    std::vector<std::vector<match_set_t>>
    divvy_matches(const std::vector<match_set_t>& matches,
                  const BGraph1& graph1, const BGraph2& graph2,
                  const std::vector<std::pair<SubGraphInfo, SubGraphInfo>>& stitch_graphs) const;
    
    // add paths onto stitch graphs
    template<class BGraph1, class BGraph2>
    void project_paths(const BGraph1& graph1, const BGraph2& graph2,
                       std::vector<std::pair<SubGraphInfo, SubGraphInfo>>& stitch_graphs) const;
    
    void merge_stitch_chains(std::vector<anchor_t>& anchor_chain,
                             const std::vector<std::vector<anchor_t>> stitch_chains,
                             const std::vector<std::pair<SubGraphInfo, SubGraphInfo>>& stitch_graphs) const;
    
    AlignmentParameters<3> alignment_params;
    // the minimum matrix size that will trigger WFA instead of PO-POA
    size_t min_wfa_size = 10000000;
    // prune WFA positions that are this far behind the opt
    size_t wfa_pruning_dist = 25;
    // use approximate deletion alignment when one graph is this much larger than the other
    size_t deletion_alignment_ratio = 100;
    
private:
    
    static const bool debug;
    static const bool instrument;
    
    void subalign(const SubGraphInfo& extraction1, const SubGraphInfo& extraction2,
                  Alignment& stitched) const;
    
    template<int NumPW>
    Alignment do_alignment(const SubGraphInfo& extraction1, const SubGraphInfo& extraction2,
                           const AlignmentParameters<NumPW>& params) const;
    
    
    template<class BGraph>
    void do_project(const BGraph& graph, SubGraphInfo& subgraph,
                    const StepIndex& step_index) const;
    
    void do_instrument(const SubGraphInfo& extraction1, const SubGraphInfo& extraction2,
                       int64_t score1, int64_t score2, double dur1, double dur2) const;
    
    template<class BGraph>
    void log_subpath_info(const BGraph& graph1, const BGraph& graph2,
                          const std::vector<std::pair<SubGraphInfo, SubGraphInfo>>& extractions) const;
    
    static std::vector<size_t> get_logging_indexes(const std::vector<anchor_t>& anchor_chain);
    
    static std::pair<int64_t, int64_t> source_sink_minmax(const SubGraphInfo& extraction);
    
};

/*
 * Template implementations
 */

template<class BGraph1, class BGraph2, class XMerge1, class XMerge2>
std::vector<std::pair<SubGraphInfo, SubGraphInfo>>
Stitcher::extract_stitch_graphs(const std::vector<anchor_t>& anchor_chain,
                                const BGraph1& graph1, const BGraph2& graph2,
                                const SentinelTableau& tableau1, const SentinelTableau& tableau2,
                                const XMerge1& chain_merge1, const XMerge2& chain_merge2) const {
    
    size_t next_log_idx = 0;
    std::vector<size_t> logging_indexes;
    if (logging::level >= logging::Debug) {
        logging_indexes = get_logging_indexes(anchor_chain);
        
        logging::log(logging::Debug, "Extracting graphs in between chain of " + std::to_string(anchor_chain.size()) + " anchors");
    }
    
    std::vector<std::pair<SubGraphInfo, SubGraphInfo>> stitch_pairs;
    stitch_pairs.reserve(anchor_chain.size() + 1);
    
    stitch_pairs.emplace_back(extract_connecting_graph(graph1, tableau1.src_id,
                                                       anchor_chain.front().walk1.front(),
                                                       chain_merge1),
                              extract_connecting_graph(graph2, tableau2.src_id,
                                                       anchor_chain.front().walk2.front(),
                                                       chain_merge2));
    
    for (size_t i = 1; i < anchor_chain.size(); ++i) {
        
        if (next_log_idx < logging_indexes.size() && i == logging_indexes[next_log_idx]) {
            logging::log(logging::Debug, "Graph extraction iteration " + std::to_string(i + 1) + " of " + std::to_string(anchor_chain.size()));
            ++next_log_idx;
        }
        
        const auto& prev_anchor = anchor_chain[i - 1];
        const auto& anchor = anchor_chain[i];
        
        stitch_pairs.emplace_back(extract_connecting_graph(graph1,
                                                           prev_anchor.walk1.back(),
                                                           anchor.walk1.front(),
                                                           chain_merge1),
                                  extract_connecting_graph(graph2,
                                                           prev_anchor.walk2.back(),
                                                           anchor.walk2.front(),
                                                           chain_merge2));
    }
    
    
    stitch_pairs.emplace_back(extract_connecting_graph(graph1, anchor_chain.back().walk1.back(),
                                                       tableau1.snk_id, chain_merge1),
                              extract_connecting_graph(graph2, anchor_chain.back().walk2.back(),
                                                       tableau2.snk_id, chain_merge2));
    
    return stitch_pairs;
}


template<class BGraph1, class BGraph2, class XMerge1, class XMerge2>
Alignment Stitcher::stitch(const std::vector<anchor_t>& anchor_chain,
                           const BGraph1& graph1, const BGraph2& graph2,
                           const SentinelTableau& tableau1, const SentinelTableau& tableau2,
                           const XMerge1& xmerge1, const XMerge2& xmerge2) const {

    size_t next_log_idx = 0;
    std::vector<size_t> logging_indexes;
    if (logging::level >= logging::Debug) {
        logging_indexes = get_logging_indexes(anchor_chain);
        
        logging::log(logging::Debug, "Stitching alignment between chain of " + std::to_string(anchor_chain.size()) + " anchors");
    }
    
    Alignment stitched;
    
    auto stitch_graphs = extract_stitch_graphs(anchor_chain, graph1, graph2,
                                               tableau1, tableau2,
                                               xmerge1, xmerge2);
    
    if (instrument) {
        // record the locations of the subgraphs
        log_subpath_info(graph1, graph2, stitch_graphs);
    }
    
    for (size_t i = 0; i < stitch_graphs.size(); ++i) {
        
        const auto& stitch_pair = stitch_graphs[i];
        
        if (next_log_idx < logging_indexes.size() && i == logging_indexes[next_log_idx]) {
            logging::log(logging::Debug, "Stitching iteration " + std::to_string(i + 1) + " of " + std::to_string(anchor_chain.size()));
            ++next_log_idx;
        }
        
        // make an intervening alignment
        if (instrument) {
            std::cerr << "subalign " << i << '\n';
        }
        subalign(stitch_pair.first, stitch_pair.second, stitched);
        
        if (i < anchor_chain.size()) {
            // copy the anchor
            const auto& anchor = anchor_chain[i];
            for (size_t j = 0; j < anchor.walk1.size(); ++j) {
                stitched.emplace_back(anchor.walk1[j], anchor.walk2[j]);
            }
        }
    }
    
    return stitched;
}

template<class BGraph1, class BGraph2>
std::vector<std::vector<match_set_t>>
Stitcher::divvy_matches(const std::vector<match_set_t>& matches,
                        const BGraph1& graph1, const BGraph2& graph2,
                        const std::vector<std::pair<SubGraphInfo, SubGraphInfo>>& stitch_graphs) const {
    
    static const bool debug = false;
    
    // identify which nodes are in the stitch graphs
    std::vector<std::pair<size_t, uint64_t>> forward_trans1(graph1.node_size(), std::pair<size_t, uint64_t>(-1, -1));
    std::vector<std::pair<size_t, uint64_t>> forward_trans2(graph2.node_size(), std::pair<size_t, uint64_t>(-1, -1));
    for (size_t i = 0; i < stitch_graphs.size(); ++i) {
        const auto& back_trans1 = stitch_graphs[i].first.back_translation;
        for (uint64_t fwd_id1 = 0; fwd_id1 < back_trans1.size(); ++fwd_id1) {
            forward_trans1[back_trans1[fwd_id1]] = std::make_pair(i, fwd_id1);
        }
        const auto& back_trans2 = stitch_graphs[i].second.back_translation;
        for (uint64_t fwd_id2 = 0; fwd_id2 < back_trans2.size(); ++fwd_id2) {
            forward_trans2[back_trans2[fwd_id2]] = std::make_pair(i, fwd_id2);
        }
    }
    
    if (debug) {
        std::cerr << "forward translations to stitch graphs:\n";
        std::cerr << "graph1:\n";
        for (uint64_t n = 0; n < graph1.node_size(); ++n) {
            std::cerr << n << ":\t" << forward_trans1[n].first << '\t' << forward_trans1[n].second << '\n';
        }
        std::cerr << "graph2:\n";
        for (uint64_t n = 0; n < graph2.node_size(); ++n) {
            std::cerr << n << ":\t" << forward_trans2[n].first << '\t' << forward_trans2[n].second << '\n';
        }
    }
    
    std::vector<std::vector<match_set_t>> divvied(stitch_graphs.size());
    
    for (size_t i = 0; i < matches.size(); ++i) {
        const auto& match_set = matches[i];
        std::unordered_set<size_t> stitch_sets_initialized;
        for (const auto& walk1 : match_set.walks1) {
            size_t stitch_idx = forward_trans1[walk1.front()].first;
            if (stitch_idx != -1 && stitch_idx == forward_trans1[walk1.back()].first) {
                // both ends of the path are in the same stitch graph
                
                if (!stitch_sets_initialized.count(stitch_idx)) {
                    // we haven't initialized a corresponding match set for the stitch graph yet
                    divvied[stitch_idx].emplace_back();
                    // the counts get retained from the original match sets
                    divvied[stitch_idx].back().count1 = match_set.count1;
                    divvied[stitch_idx].back().count2 = match_set.count2;
                    stitch_sets_initialized.insert(stitch_idx);
                }
                
                // copy and translate the walk
                divvied[stitch_idx].back().walks1.emplace_back();
                auto& stitch_walk = divvied[stitch_idx].back().walks1.back();
                for (auto node_id : walk1) {
                    stitch_walk.push_back(forward_trans1[node_id].second);
                }
            }
        }
        
        for (const auto& walk2 : match_set.walks2) {
            size_t stitch_idx = forward_trans2[walk2.front()].first;
            if (stitch_sets_initialized.count(stitch_idx) && stitch_idx == forward_trans2[walk2.back()].first) {
                // both ends of the path are in the same stitch graph
                // note: checking whether it's initialized also checks whether it's -1 because
                // of the previous loop
                
                // copy and translate the walk
                divvied[stitch_idx].back().walks2.emplace_back();
                auto& stitch_walk = divvied[stitch_idx].back().walks2.back();
                for (auto node_id : walk2) {
                    stitch_walk.push_back(forward_trans2[node_id].second);
                }
            }
        }
        
        // clear out any initialized sets that didn't end up getting a walk from graph2
        for (auto stitch_idx : stitch_sets_initialized) {
            if (divvied[stitch_idx].back().walks2.empty()) {
                divvied[stitch_idx].pop_back();
            }
        }
    }
    
    return divvied;
}

template<class BGraph>
void Stitcher::do_project(const BGraph& graph, SubGraphInfo& subgraph,
                          const StepIndex& step_index) const {
    
    std::unordered_map<uint64_t, uint64_t> path_ids;
    for (auto node_id : topological_order(subgraph.subgraph)) {
        for (const auto& step : step_index.path_steps(subgraph.back_translation[node_id])) {
            auto it = path_ids.find(step.first);
            if (it == path_ids.end()) {
                auto new_path_id = subgraph.subgraph.add_path(graph.path_name(step.first));
                it = path_ids.insert(std::make_pair(step.first, new_path_id)).first;
            }
            subgraph.subgraph.extend_path(it->second, node_id);
        }
    }
}

template<class BGraph1, class BGraph2>
void Stitcher::project_paths(const BGraph1& graph1, const BGraph2& graph2,
                             std::vector<std::pair<SubGraphInfo, SubGraphInfo>>& stitch_graphs) const {
    
    StepIndex step_index1(graph1);
    StepIndex step_index2(graph2);
    for (auto& stitch_pair : stitch_graphs) {
        do_project(graph1, stitch_pair.first, step_index1);
        do_project(graph2, stitch_pair.second, step_index2);
    }
}

template<int NumPW>
Alignment Stitcher::do_alignment(const SubGraphInfo& extraction1, const SubGraphInfo& extraction2,
                                 const AlignmentParameters<NumPW>& params) const {
    
//    int pdist = 5;
//    if (params.gap_open.size() >= 2) {
//        // twice the distance where the final gap param becomes active
//        pdist = (2 * (params.gap_open.back() - params.gap_open[params.gap_open.size() - 2])
//                 / (params.gap_extend[params.gap_extend.size() - 2] - params.gap_extend.back()));
//    }
//
//    auto begin = std::chrono::high_resolution_clock::now();
//    auto inter_aln = pwfa_po_poa(extraction1.subgraph, extraction2.subgraph,
//                                 extraction1.sources, extraction2.sources,
//                                 extraction1.sinks, extraction2.sinks, params, 2 * wfa_pruning_dist);
//    auto middle = std::chrono::high_resolution_clock::now();
//    auto inter_aln = po_poa(extraction1.subgraph, extraction2.subgraph,
//                            extraction1.sources, extraction2.sources,
//                            extraction1.sinks, extraction2.sinks, params);
//    auto end = std::chrono::high_resolution_clock::now();
//
//    double po_poa_dur = std::chrono::duration<double, std::nano>(end - middle).count();
//    double wfa_dur = std::chrono::duration<double, std::nano>(middle - begin).count();
    
    size_t mat_size = (extraction1.subgraph.node_size() + 1) * (extraction2.subgraph.node_size() + 1);
    if (instrument) {
        std::cerr << '%' << '\t' << extraction1.subgraph.node_size() << '\t' << extraction2.subgraph.node_size() << '\t' << mat_size;
    }

    auto begin = std::chrono::high_resolution_clock::now();
    Alignment inter_aln;
    
    // TODO: should i use shortest/longest path instead of node count to identify deletions?
    
    if (extraction2.subgraph.node_size() == 0) {
        // align with shortest path
        inter_aln = std::move(pure_deletion_alignment(extraction1.subgraph,
                                                      extraction1.sources,
                                                      extraction1.sinks,
                                                      params));
    }
    else if (extraction1.subgraph.node_size() == 0) {
        // align with shortest path
        inter_aln = std::move(pure_deletion_alignment(extraction2.subgraph,
                                                      extraction2.sources,
                                                      extraction2.sinks,
                                                      params));
        swap_graphs(inter_aln);
    }
    else if (extraction1.subgraph.node_size() * deletion_alignment_ratio <= extraction2.subgraph.node_size()) {
        // graph1 is probably mostly a deletion of graph2
        inter_aln = std::move(deletion_wfa_po_poa(extraction1.subgraph, extraction2.subgraph,
                                                  extraction1.sources, extraction2.sources,
                                                  extraction1.sinks, extraction2.sinks, params));
    }
    else if (extraction2.subgraph.node_size() * deletion_alignment_ratio <= extraction1.subgraph.node_size()) {
        // graph2 is probably mostly a deletion of graph1
        inter_aln = std::move(deletion_wfa_po_poa(extraction2.subgraph, extraction1.subgraph,
                                                  extraction2.sources, extraction1.sources,
                                                  extraction2.sinks, extraction1.sinks, params));
        swap_graphs(inter_aln);
    }
    else if (mat_size > min_wfa_size) {
        inter_aln = std::move(pwfa_po_poa(extraction1.subgraph, extraction2.subgraph,
                                          extraction1.sources, extraction2.sources,
                                          extraction1.sinks, extraction2.sinks, params,
                                          2 * wfa_pruning_dist));
    }
    else {
        inter_aln = std::move(po_poa(extraction1.subgraph, extraction2.subgraph,
                                     extraction1.sources, extraction2.sources,
                                     extraction1.sinks, extraction2.sinks, params));
    }
    
    if (instrument) {
        auto end = std::chrono::high_resolution_clock::now();
        double dur = std::chrono::duration<double, std::nano>(end - begin).count();
        std::cerr << '\t' << dur << '\n';
        
//        int64_t wfa_score = score_alignment(extraction1.subgraph, extraction2.subgraph,
//                                            dummy, params);
//        int64_t po_poa_score = score_alignment(extraction1.subgraph, extraction2.subgraph,
//                                               inter_aln, params);
//        int64_t score = score_alignment(extraction1.subgraph, extraction2.subgraph, inter_aln, params);
//
//        // note: includes some dummy measurements
//        do_instrument(extraction1, extraction2, score, 0, dur, 0.0);
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
