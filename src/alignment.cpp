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

void translate(Alignment& alignment,
               const vector<uint64_t>& back_translation1,
               const vector<uint64_t>& back_translation2) {
    
    for (size_t i = 0; i < alignment.size(); ++i) {
        auto& aln_pair = alignment[i];
        if (aln_pair.node_id1 != AlignedPair::gap) {
            aln_pair.node_id1 = back_translation1[aln_pair.node_id1];
        }
        if (aln_pair.node_id2 != AlignedPair::gap) {
            aln_pair.node_id2 = back_translation2[aln_pair.node_id2];
        }
    }
}

std::string cigar(const Alignment& alignment) {
    // TODO: copypasta from explicit version
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
        else {
            op = 'M';
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
