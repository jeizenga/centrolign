#include "centrolign/superbubbles.hpp"

namespace centrolign {

using namespace std;

std::vector<std::pair<uint64_t, bool>> SuperbubbleTree::postorder() const {
    
    std::vector<std::pair<uint64_t, bool>> result;
    result.reserve(chain_size() + superbubble_size());
    
    for (uint64_t chain_id = 0; chain_id < chain_size(); ++chain_id) {
        if (superbubble_containing(chain_id) != -1) {
            // this is not a top-level chain
            continue;
        }
        
        // stack of (id, is chain, children have been added)
        std::vector<std::tuple<uint64_t, bool, bool>> stack;
        stack.emplace_back(chain_id, true, false);
        
        while (!stack.empty()) {
            if (std::get<2>(stack.back())) {
                result.emplace_back(std::get<0>(stack.back()),
                                    std::get<1>(stack.back()));
                stack.pop_back();
            }
            else {
                // add the children to the stack
                std::get<2>(stack.back()) = true;
                if (std::get<1>(stack.back())) {
                    // is a chain, add superbubble children
                    for (auto child_id : superbubbles_inside(std::get<0>(stack.back()))) {
                        stack.emplace_back(child_id, false, false);
                    }
                }
                else {
                    // is a superbubble, add chain children
                    for (auto child_id : chains_inside(std::get<0>(stack.back()))) {
                        stack.emplace_back(child_id, true, false);
                    }
                }
            }
        }
    }
    
    return result;
}

}
