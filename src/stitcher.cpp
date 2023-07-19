#include "centrolign/stitcher.hpp"

namespace centrolign {

using namespace std;

const bool Stitcher::debug = false;
const bool Stitcher::instrument = true;

Stitcher::Stitcher() {
    alignment_params.match = 4;
    alignment_params.mismatch = 8;
    alignment_params.gap_open[0] = 6;
    alignment_params.gap_extend[0] = 4;
    alignment_params.gap_open[1] = 150;
    alignment_params.gap_extend[1] = 1;
}

void Stitcher::do_instrument(const SubGraphInfo& extraction1,
                             const SubGraphInfo& extraction2,
                             int64_t score) const {
    
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
    std::cerr << '#' << '\t' << extraction1.subgraph.node_size() << '\t' << extraction2.subgraph.node_size() << '\t' << (extraction1.subgraph.node_size() * extraction2.subgraph.node_size()) << '\t' << score << '\t' << approx_max_match << '\t' << approx_max_score << '\t' << min1 << '\t' << max1 << '\t' << min2 << '\t' << max2 << '\n';
}

}
