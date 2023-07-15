#include "centrolign/stitcher.hpp"

namespace centrolign {

using namespace std;

const bool Stitcher::debug = false;
const bool Stitcher::instrument = false;

Stitcher::Stitcher() {
    alignment_params.match = 4;
    alignment_params.mismatch = 8;
    alignment_params.gap_open[0] = 6;
    alignment_params.gap_extend[0] = 4;
    alignment_params.gap_open[1] = 150;
    alignment_params.gap_extend[1] = 1;
}

}
