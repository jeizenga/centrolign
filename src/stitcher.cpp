#include "centrolign/stitcher.hpp"

#include <chrono>
#include <stdexcept>

namespace centrolign {

using namespace std;

const bool Stitcher::debug = false;
const bool Stitcher::instrument = false;

Stitcher::Stitcher() {
    alignment_params.match = 20;
    alignment_params.mismatch = 40;
    alignment_params.gap_open[0] = 30;
    alignment_params.gap_extend[0] = 20;
    alignment_params.gap_open[1] = 800;
    alignment_params.gap_extend[1] = 5;
    alignment_params.gap_open[2] = 2500;
    alignment_params.gap_extend[2] = 1;
}

std::vector<size_t> Stitcher::get_logging_indexes(const std::vector<anchor_t>& anchor_chain) {
    std::vector<size_t> logging_indexes;
    for (size_t i = 1; i < 10; ++i) {
        logging_indexes.push_back((anchor_chain.size() * i) / 10);
    }
    auto end = std::unique(logging_indexes.begin(), logging_indexes.end());
    logging_indexes.resize(end - logging_indexes.begin());
    return logging_indexes;
}

pair<int64_t, int64_t> Stitcher::source_sink_minmax(const SubGraphInfo& extraction) {
    
    auto dists = minmax_distance(extraction.subgraph, &extraction.sources);
    pair<int64_t, int64_t> mm(numeric_limits<int64_t>::max(), -1);
    for (auto node_id : extraction.sinks) {
        mm.first = min(mm.first, dists[node_id].first);
        mm.second = max(mm.second, dists[node_id].second);
    }
    return mm;
}

void Stitcher::subalign(const SubGraphInfo& extraction1,
                        const SubGraphInfo& extraction2,
                        Alignment& stitched) const {
    
    // TODO: it's wasteful doing this every time, but i expect it not to dominate
    
    // find out what gap size makes you need the larger gap parameters
    std::vector<size_t> cutoffs(alignment_params.gap_extend.size() - 1);
    for (size_t i = 1; i < alignment_params.gap_open.size(); ++i) {
        // the algebra is much easier if we get to assume a fixed order
        if (alignment_params.gap_open[i - 1] > alignment_params.gap_open[i] ||
            alignment_params.gap_extend[i - 1] < alignment_params.gap_extend[i]) {
            throw std::runtime_error("Affine gap parameters must be increasing in gap open penalty and decreasing in gap extend penalty");
        }
        
        auto diff_open = alignment_params.gap_open[i] - alignment_params.gap_open[i - 1];
        auto diff_extend = alignment_params.gap_extend[i - 1] - alignment_params.gap_extend[i];
        
        // round down to the nearest integer
        cutoffs[i - 1] = (diff_open + diff_extend - 1) / diff_extend;
    }
    
    Alignment inter_aln;
    if (extraction1.subgraph.node_size() < cutoffs.front() &&
        extraction2.subgraph.node_size() < cutoffs.front()) {
        // we can tell without doing min/max distance that these won't use the
        auto trunc_params = truncate_parameters<3, 1>(alignment_params);
        inter_aln = std::move(do_alignment(extraction1, extraction2, trunc_params));
    }
    else {
        // TODO: is this worth doing, or should i always just use the node size?
        int64_t max_gap = 0;
        if (extraction1.subgraph.node_size() != 0) {
            max_gap = source_sink_minmax(extraction1).second;
        }
        if (extraction2.subgraph.node_size() != 0) {
            max_gap = max(max_gap, source_sink_minmax(extraction1).second);
        }
        size_t c = 0;
        while (c < cutoffs.size() && max_gap >= cutoffs[c]) {
            ++c;
        }
        
        if (c == 0) {
            auto trunc_params = truncate_parameters<3, 1>(alignment_params);
            inter_aln = std::move(do_alignment(extraction1, extraction2, trunc_params));
        }
        else if (c == 1) {
            auto trunc_params = truncate_parameters<3, 2>(alignment_params);
            inter_aln = std::move(do_alignment(extraction1, extraction2, trunc_params));
        }
        else {
            inter_aln = std::move(do_alignment(extraction1, extraction2, alignment_params));
        }
    }
    
    translate(inter_aln, extraction1.back_translation, extraction2.back_translation);
    
    if (debug) {
        std::cerr << "inter-anchor alignment:\n";
        for (auto& ap : inter_aln) {
            std::cerr << '\t' << (int64_t) ap.node_id1 << '\t' << (int64_t) ap.node_id2 << '\n';
        }
    }
    
    for (auto aln_pair : inter_aln) {
        stitched.push_back(aln_pair);
    }
}

void Stitcher::do_instrument(const SubGraphInfo& extraction1,
                             const SubGraphInfo& extraction2,
                             int64_t score1, int64_t score2,
                             double dur1, double dur2) const {
    
    int64_t min1 = -1;
    int64_t max1 = -1;
    int64_t min2 = -1;
    int64_t max2 = -1;
    if (!extraction1.back_translation.empty()) {
        min1 = *std::min_element(extraction1.back_translation.begin(),
                                 extraction1.back_translation.end());
        max1 = *std::max_element(extraction1.back_translation.begin(),
                                 extraction1.back_translation.end());
    }
    if (!extraction2.back_translation.empty()) {
        min2 = *std::min_element(extraction2.back_translation.begin(),
                                 extraction2.back_translation.end());
        max2 = *std::max_element(extraction2.back_translation.begin(),
                                 extraction2.back_translation.end());
    }
    int64_t min_penalty = std::numeric_limits<int64_t>::max();
    int64_t approx_gap = std::abs((int64_t) (extraction1.subgraph.node_size() - extraction2.subgraph.node_size()));
    for (size_t i = 0; i < alignment_params.gap_open.size(); ++i) {
        min_penalty = min<int64_t>(min_penalty,
                                   alignment_params.gap_extend[i] * approx_gap + alignment_params.gap_open[i]);
    }
    int64_t approx_max_match = alignment_params.match * min(extraction1.subgraph.node_size(),
                                                            extraction2.subgraph.node_size());
    int64_t approx_max_score = approx_max_match - min_penalty;
    std::cerr << '#' << '\t' << extraction1.subgraph.node_size() << '\t' << extraction2.subgraph.node_size() << '\t' << ((extraction1.subgraph.node_size() + 1) * (extraction2.subgraph.node_size() + 1)) << '\t' << score1 << '\t' << score2 << '\t' << approx_max_match << '\t' << approx_max_score << '\t' << min1 << '\t' << max1 << '\t' << min2 << '\t' << max2 << '\t' << dur1 << '\t' << dur2 << '\n';
}

}
