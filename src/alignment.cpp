#include "centrolign/alignment.hpp"

namespace centrolign {

using namespace std;

const uint64_t AlignedPair::gap = -1;

AlignedPair::AlignedPair(uint64_t node_id1, uint64_t node_id2) :
    node_id1(node_id1),
    node_id2(node_id2) {
    // nothing else to do
}

bool AlignedPair::operator==(const AlignedPair& other) const {
    return node_id1 == other.node_id1 && node_id2 == other.node_id2;
}

}
