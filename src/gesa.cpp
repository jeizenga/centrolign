#include "centrolign/gesa.hpp"

namespace centrolign {

using namespace std;


size_t GESA::component_size() const {
    return component_sentinels.size();
}

char GESA::sentinel(size_t component, bool beginning) const {
    return beginning ? component_sentinels[component].first : component_sentinels[component].second;
}

GESANode GESA::root() const {
    return GESANode(0, lcp_array.size());
}

void GESA::construct_child_array() {
    
    child_array.resize(lcp_array.size(), -1);
    
    // get the next l-indexes
    vector<size_t> stack;
    for (size_t i = 0; i < lcp_array.size(); ++i) {
        // ensure invariant that LCPs are increasing over the stacked indexes, with intervening
        // LCPs being strictly longer
        while (!stack.empty() && lcp_array[stack.back()] > lcp_array[i]) {
            stack.pop_back();
        }
        if (!stack.empty() && lcp_array[i] == lcp_array[stack.back()]) {
            // we are at the next l-index
            child_array[stack.back()] = i;
            // clear this index out of the way now that we've found its next l-index
            stack.pop_back();
        }
        stack.push_back(i);
    }
    
    // get the up and down indexes
    stack.clear();
    for (size_t i = 0; i < lcp_array.size(); ++i) {
        size_t last_idx = -1;
        while (!stack.empty() && lcp_array[stack.back()] > lcp_array[i]) {
            last_idx = stack.back();
            stack.pop_back();
            if (!stack.empty() &&
                child_array[stack.back()] == -1 && // not already holding a next l-index
                lcp_array[i] == lcp_array[stack.back()] && // is latest, later blocked by equality at i
                lcp_array[last_idx] != lcp_array[stack.back()]) { // inequality from increasing stack invariant is strict
                // set childtab[stack.back()].down
                child_array[stack.back()] = last_idx;
            }
        }
        if (last_idx != -1) {
            // final item popped is earliest, must be strict inequality, all intervening have longer LCP
            // by the stack invariant
            // set childtab[i].up (but stored in position i - 1, which is guaranteed to be empty)
            child_array[i - 1] = last_idx;
        }
        stack.push_back(i);
    }
}

}
