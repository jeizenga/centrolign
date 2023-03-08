#ifndef centrolign_graph_hpp
#define centrolign_graph_hpp

#include <vector>
#include <cstdint>
#include <string>

namespace centrolign {

class SequenceGraph {
public:
    
    SequenceGraph() = default;
    ~SequenceGraph() = default;
    
    uint64_t add_node(const std::string& sequence);
    void add_edge(uint64_t node_id_from, uint64_t node_id_to);
    
    size_t node_size() const;
    std::string sequence(uint64_t node_id) const;
    const std::vector<uint64_t>& next(uint64_t node_id) const;
    const std::vector<uint64_t>& previous(uint64_t node_id) const;
    size_t next_size(uint64_t node_id) const;
    size_t previous_size(uint64_t node_id) const;
    
private:
    
    struct SequenceGraphNode {
        
        SequenceGraphNode(const std::string& sequence) : seq(sequence) { }
        SequenceGraphNode() = default;
        ~SequenceGraphNode() = default;
        
        std::string seq;
        std::vector<uint64_t> next;
        std::vector<uint64_t> prev;
        
    };
    
    std::vector<SequenceGraphNode> nodes;
};

class BaseGraphOverlay {
public:
    BaseGraphOverlay(const SequenceGraph* sequence_graph);
    BaseGraphOverlay() = default;
    ~BaseGraphOverlay() = default;
    
    size_t node_size() const;
    char base(uint64_t node_id) const;
    std::vector<uint64_t> next(uint64_t node_id) const;
    std::vector<uint64_t> previous(uint64_t node_id) const;
    size_t next_size(uint64_t node_id) const;
    size_t previous_size(uint64_t node_id) const;
    
private:
    
    const SequenceGraph* seq_graph = nullptr;
    std::vector<size_t> cumul_len;
    std::vector<uint64_t> origin;
}

class BaseGraph {
public:
    BaseGraph() = default;
    ~BaseGraph() = default;
    
    uint64_t add_node(char base);
    void add_edge(uint64_t node_id_from, uint64_t node_id_to);
    
    size_t node_size() const;
    char base(uint64_t node_id) const;
    std::vector<uint64_t> next(uint64_t node_id) const;
    std::vector<uint64_t> previous(uint64_t node_id) const;
    size_t next_size(uint64_t node_id) const;
    size_t previous_size(uint64_t node_id) const;
    
private:
    
    struct BaseGraphNode {
        
        BaseGraphNode(char base) : base(base) { }
        BaseGraphNode() = default;
        ~BaseGraphNode() = default;
        
        char base = 0;
        std::vector<uint64_t> next;
        std::vector<uint64_t> prev;
        
    };
    
    std::vector<BaseGraphNode> nodes;
}

}

#endif /* centrolign_graph_hpp */
