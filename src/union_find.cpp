#include "centrolign/union_find.hpp"

namespace centrolign {

UnionFind::UnionFind(size_t size) noexcept : nodes(size) {
    for (size_t i = 0; i < nodes.size(); i++) {
        nodes[i].head = i;
    }
}

}
