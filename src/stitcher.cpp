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
    auto min1 = *std::min_element(extraction1.back_translation.begin(),
                                  extraction1.back_translation.end());
    auto max1 = *std::max_element(extraction1.back_translation.begin(),
                                  extraction1.back_translation.end());
    auto min2 = *std::min_element(extraction2.back_translation.begin(),
                                  extraction2.back_translation.end());
    auto max2 = *std::max_element(extraction2.back_translation.begin(),
                                  extraction2.back_translation.end());
    std::cerr << '#' << '\t' << extraction1.subgraph.node_size() << '\t' << extraction2.subgraph.node_size() << '\t' << (extraction1.subgraph.node_size() * extraction2.subgraph.node_size()) << '\t' << score << '\t' << min1 << '\t' << max1 << '\t' << min2 << '\t' << max2 << '\n';
}

}
