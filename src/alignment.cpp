#include "centrolign/alignment.hpp"

#include <unordered_map>

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

void swap_graphs(Alignment& alignment) {
    for (auto& aln_pair : alignment) {
        std::swap(aln_pair.node_id1, aln_pair.node_id2);
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

std::string explicit_cigar(const Alignment& alignment,
                           const std::string& seq1, const std::string& seq2) {
    // TODO: copypasta from graph version
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
        else if (seq1.at(aln_pair.node_id1) == seq2.at(aln_pair.node_id2)) {
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

Alignment induced_pairwise_alignment(const BaseGraph& graph, uint64_t path_id1, uint64_t path_id2) {
        
    unordered_map<uint64_t, uint64_t> index_in_path1;
    
    const auto& path1 = graph.path(path_id1);
    const auto& path2 = graph.path(path_id2);
    
    for (uint64_t i = 0; i < path1.size(); ++i) {
        index_in_path1[path1[i]] = i;
    }
    
    Alignment alignment;
    
    // scan path 2
    uint64_t j = 0;
    for (uint64_t i = 0; i < path2.size(); ++i) {
        auto node_id = path2[i];
        auto it = index_in_path1.find(node_id);
        if (it == index_in_path1.end()) {
            // no aligned base, must be a gap
            alignment.emplace_back(AlignedPair::gap, i);
        }
        else {
            // clear out earlier gaps in path1
            while (j < it->second) {
                alignment.emplace_back(j++, AlignedPair::gap);
            }
            // record the aligned pair
            alignment.emplace_back(j++, i);
        }
    }
    
    // finish off path 1, if necessary
    while (j < path1.size()) {
        alignment.emplace_back(j++, AlignedPair::gap);
    }
    
    // consolidate mismatches equal length gap runs runs as a mismatch
    size_t removed = 0;
    for (size_t i = 0; i < alignment.size();) {
        auto& aln_pair = alignment[i];
        if (aln_pair.node_id1 != AlignedPair::gap && aln_pair.node_id2 != AlignedPair::gap) {
            alignment[i - removed] = aln_pair;
            ++i;
        }
        else {
            // find the extent of the gap run
            size_t j = i;
            size_t gaps1 = 0;
            size_t gaps2 = 0;
            while (j < alignment.size() && (alignment[j].node_id1 == AlignedPair::gap ||
                                            alignment[j].node_id2 == AlignedPair::gap)) {
                gaps1 += (alignment[j].node_id1 == AlignedPair::gap);
                gaps2 += (alignment[j].node_id2 == AlignedPair::gap);
                ++j;
            }
            
            
            // find the final index for each sequence in this run
            uint64_t last1 = -1;
            uint64_t last2 = -1;
            if (i != 0) {
                // this node is guaranteed to be a match or it would have
                // ended up in this run
                auto& last_pair = alignment[i - removed - 1];
                last1 = last_pair.node_id1;
                last2 = last_pair.node_id2;
            }
            
            if (gaps1 == gaps2) {
                // they have the same length, so we'll call this a mismatch
                for (uint64_t n = 0; n < gaps1; ++n) {
                    alignment[i - removed + n] = AlignedPair(last1 + n + 1, last2 + n + 1);
                }
                removed += gaps1;
            }
            else {
                // since they're not the same length, we interpret this as a separate insertion
                // and deletion
                
                // consolidate into a single gap for each sequence
                for (uint64_t n = 0; n < gaps2; ++n) {
                    alignment[i - removed + n] = AlignedPair(last1 + n + 1, AlignedPair::gap);
                }
                for (uint64_t n = 0; n < gaps1; ++n) {
                    alignment[i - removed + n + gaps2] = AlignedPair(AlignedPair::gap, last2 + n + 1);
                }
            }
            i = j;
        }
    }
    alignment.resize(alignment.size() - removed);
    
    return alignment;
}

}
