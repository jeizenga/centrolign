#include "centrolign/determinize.hpp"

namespace centrolign {

using namespace std;


SentinelTableau translate_tableau(const BaseGraph& determinized,
                                  const SentinelTableau& original_tableau) {
    
    SentinelTableau translated;
    translated.src_sentinel = original_tableau.src_sentinel;
    translated.snk_sentinel = original_tableau.snk_sentinel;
    for (uint64_t node_id = 0; node_id < determinized.node_size(); ++node_id) {
        if (determinized.label(node_id) == translated.src_sentinel) {
            translated.src_id = node_id;
        }
        if (determinized.label(node_id) == translated.snk_sentinel) {
            translated.snk_id = node_id;
        }
    }
    return translated;
}


}
