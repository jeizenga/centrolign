#ifndef centrolign_union_find_hpp
#define centrolign_union_find_hpp

#include <vector>
#include <cstdint>
#include <algorithm>

namespace centrolign {

/**
 * A union-find / disjoint-set data structure
 */
class UnionFind {
public:

    UnionFind(size_t size) noexcept;
    UnionFind() noexcept = default;
    UnionFind(const UnionFind& other) noexcept = default;
    UnionFind(UnionFind&& other) noexcept = default;
    ~UnionFind() = default;
    UnionFind& operator=(const UnionFind& other) noexcept = default;
    UnionFind& operator=(UnionFind&& other) noexcept = default;
    
    // returns the group ID that index i belongs to (can change after calling union)
    inline size_t find_group(size_t i);
    
    // merges the group containing index i with the group containing index j
    inline void union_groups(size_t i, size_t j);
    
    // Returns all of the groups, each in a separate vector
    inline std::vector<std::vector<size_t>> get_groups();
    
private:
    
    struct Node {
        Node() = default;
        ~Node() = default;
        size_t rank = 0;
        size_t head = 0;
    };
    std::vector<Node> nodes;
};

inline size_t UnionFind::find_group(size_t i) {
    std::vector<size_t> path;
    // traverse tree upwards
    while (nodes[i].head != i) {
        path.push_back(i);
        i = nodes[i].head;
    }
    // compress path
    for (size_t p = 1; p < path.size(); p++) {
        size_t j = path[p - 1];
        nodes[j].head = i;
    }
    return i;
}

inline void UnionFind::union_groups(size_t i, size_t j) {
    size_t head_i = find_group(i);
    size_t head_j = find_group(j);
    if (head_i != head_j) {
        // use rank to determine which group to make the head
        auto& node_i = nodes[head_i];
        auto& node_j = nodes[head_j];
        if (node_i.rank > node_j.rank) {
            node_j.head = head_i;
        }
        else {
            node_i.head = head_j;
            
            if (node_j.rank == node_i.rank) {
                node_j.rank++;
            }
        }
    }
}

inline std::vector<std::vector<size_t>> UnionFind::get_groups() {
    std::vector<std::vector<size_t>> groups(nodes.size());
    for (size_t i = 0; i < nodes.size(); i++) {
        groups[find_group(i)].push_back(i);
    }
    auto new_end = std::remove_if(groups.begin(), groups.end(),
                                  [](const std::vector<size_t>& grp) { return grp.empty(); });
    groups.resize(new_end - groups.begin());
    return groups;
}


}

#endif /* centrolign_union_find_hpp */
