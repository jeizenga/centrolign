#ifndef centrolign_trie_hpp
#define centrolign_trie_hpp

#include <vector>
#include <cstdint>
#include <string>
#include <unordered_map>

namespace centrolign {

/*
 * Standard trie data structure
 */
class Trie {
public:
    Trie();
    Trie(Trie&& other) noexcept = default;
    Trie& operator=(Trie&& other) noexcept = default;
    ~Trie() = default;
    
    // returns the path ID
    uint64_t insert_sequence(const std::string& name,
                             const std::string& sequence);
    
    // query interface
    size_t node_size() const;
    uint64_t get_root() const;
    uint64_t get_parent(uint64_t node_id) const;
    std::vector<uint64_t> next(uint64_t node_id) const;
    size_t next_size(uint64_t node_id) const;
    std::vector<uint64_t> previous(uint64_t node_id) const;
    size_t previous_size(uint64_t node_id) const;
    char label(uint64_t node_id) const;
    uint64_t follow(uint64_t node_id, char label) const;
    
    size_t path_size() const;
    const std::string& path_name(uint64_t path_id) const;
    const std::vector<uint64_t>& path(uint64_t path_id) const;
    
protected:
    
    struct Node {
        char label = '\0';
        uint64_t parent = -1;
        std::unordered_map<char, uint64_t> children;
        
        Node() = default;
        ~Node() = default;
    };
    
    std::vector<Node> nodes;
    
    std::vector<std::pair<std::string, std::vector<uint64_t>>> paths;
};

}

#endif /* centrolign_trie_hpp */
