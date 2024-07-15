#ifndef centrolign_union_find_hpp
#define centrolign_union_find_hpp

#include <vector>
#include <cstdint>

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
    size_t find_group(size_t i);
    
    // merges the group containing index i with the group containing index j
    void union_groups(size_t i, size_t j);
    
    // Returns all of the groups, each in a separate vector
    std::vector<std::vector<size_t>> get_groups();
    
private:
    
    struct Node {
        Node() = default;
        ~Node() = default;
        size_t rank = 0;
        size_t head = 0;
    };
    std::vector<Node> nodes;
};

}

#endif /* centrolign_union_find_hpp */
