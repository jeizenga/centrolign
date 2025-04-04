#ifndef centrolign_tree_hpp
#define centrolign_tree_hpp

#include <vector>
#include <limits>
#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_map>

namespace centrolign {

// make a Newick sring that
std::string in_order_newick_string(const std::vector<std::string>& sequences);

/*
 * Object that connects an anchor chain into a base-level alignment
 */
class Tree {
public:
    
    // construct from a newick string
    Tree(const std::string& newick);
    Tree() = default;
    Tree(const Tree& other) = default;
    Tree(Tree&& tree) = default;
    ~Tree() = default;
    Tree& operator=(const Tree& other) = default;
    Tree& operator=(Tree&& other) = default;
    
    // query interface
    size_t node_size() const;
    bool has_label(const std::string& label) const;
    uint64_t get_id(const std::string& label) const;
    uint64_t get_root() const;
    uint64_t get_parent(uint64_t node_id) const;
    const std::vector<uint64_t>& get_children(uint64_t node_id) const;
    const std::string& label(uint64_t node_id) const;
    double distance(uint64_t node_id) const;
    bool is_leaf(uint64_t node_id) const;
    
    // convert all polytomies into an arbitrary sequence of binary nodes
    void binarize();
    
    // convert 0-distance children into polytomies
    void polytomize();
    
    // remove all nodes that are not in this set or one of their ancestors
    void prune(const std::vector<uint64_t>& node_ids);
    
    // compress non-branching paths into a single edge
    void compact();
    
    // node IDs in pre order
    std::vector<uint64_t> preorder() const;
    
    // node IDs in post order
    std::vector<uint64_t> postorder() const;
    
    // node IDs in a leaf-to-root order that places small subtrees first
    std::vector<uint64_t> small_first_postorder() const;
    
    // get Newick string encoding this tree
    std::string to_newick() const;
    
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
    
    std::string::const_iterator find_skipping_quotes(std::string::const_iterator begin,
                                                     std::string::const_iterator end,
                                                     const std::string& values) const;
    
    uint64_t add_child(uint64_t parent);
    void parse_label(uint64_t node_id,
                     std::string::const_iterator begin,
                     std::string::const_iterator end);
    void filter(const std::vector<bool>& keep);
    
    uint64_t root = -1;
    
    std::vector<Node> nodes;
    std::unordered_map<std::string, uint64_t> label_map;
};

}

#endif /* centrolign_tree_hpp */
