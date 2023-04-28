#ifndef centrolign_tree_hpp
#define centrolign_tree_hpp

#include <vector>
#include <limits>
#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_map>

namespace centrolign {

/*
 * Object that connects an anchor chain into a base-level alignment
 */
class Tree {
public:
    
    // construct from a newick string
    Tree(const std::string& newick);
    Tree() = default;
    ~Tree() = default;
    
    // query interface
    uint64_t get_id(const std::string& label) const;
    uint64_t get_root() const;
    uint64_t get_parent(uint64_t node_id) const;
    const std::vector<uint64_t>& get_children(uint64_t node_id) const;
    const std::string& label(uint64_t node_id) const;
    
    // convert all polytomies into an arbitrary sequence of binary nodes
    void binarize();
    
    // node IDs in post order
    std::vector<uint64_t> postorder() const;
    
protected:
    
    struct Node {
        std::string label;
        uint64_t parent = -1;
        double distance = std::numeric_limits<double>::max();
        std::vector<uint64_t> children;
    };
    
    const static bool debug;
    
    void debug_print(std::ostream& out) const;
    void debug_print_internal(uint64_t node, size_t level, std::ostream& out) const;
    
    uint64_t add_child(uint64_t parent);
    void parse_label(uint64_t node_id,
                     std::string::const_iterator begin,
                     std::string::const_iterator end);
    
    uint64_t root = -1;
    
    std::vector<Node> nodes;
    std::unordered_map<std::string, uint64_t> label_map;
};

}

#endif /* centrolign_tree_hpp */
