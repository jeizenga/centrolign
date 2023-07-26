#ifndef centrolign_graph_hpp
#define centrolign_graph_hpp

#include <vector>
#include <cstdint>
#include <string>
#include <utility>
#include <unordered_map>

namespace centrolign {

class SequenceGraph {
public:
    
    SequenceGraph() noexcept = default;
    SequenceGraph(const SequenceGraph& other) noexcept = default;
    SequenceGraph(SequenceGraph&& other) noexcept = default;
    ~SequenceGraph() noexcept = default;
    
    SequenceGraph& operator=(const SequenceGraph& other) noexcept = default;
    SequenceGraph& operator=(SequenceGraph&& other) noexcept = default;
    
    uint64_t add_node(const std::string& sequence);
    void add_edge(uint64_t node_id_from, uint64_t node_id_to);
    void remove_edge(uint64_t node_id_from, uint64_t node_id_to);
    void relabel(uint64_t node_id, const std::string& sequence);
    
    size_t node_size() const;
    std::string label(uint64_t node_id) const;
    size_t label_size(uint64_t node_id) const;
    const std::vector<uint64_t>& next(uint64_t node_id) const;
    const std::vector<uint64_t>& previous(uint64_t node_id) const;
    size_t next_size(uint64_t node_id) const;
    size_t previous_size(uint64_t node_id) const;
    
    uint64_t add_path(const std::string& name);
    void extend_path(uint64_t path_id, uint64_t node_id);
    void pre_extend_path(uint64_t path_id, uint64_t node_id);
    
    size_t path_size() const;
    const std::string& path_name(uint64_t path_id) const;
    const std::vector<uint64_t>& path(uint64_t path_id) const;
    uint64_t path_id(const std::string& name) const;
    
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
    std::unordered_map<std::string, uint64_t> name_to_id;
    std::vector<std::pair<std::string, std::vector<uint64_t>>> paths;
};

class BaseGraphOverlay {
public:
    BaseGraphOverlay(const SequenceGraph* sequence_graph) noexcept;
    BaseGraphOverlay() = default;
    ~BaseGraphOverlay() = default;
    
    size_t node_size() const;
    char label(uint64_t node_id) const;
    size_t label_size(uint64_t node_id) const;
    std::vector<uint64_t> next(uint64_t node_id) const;
    std::vector<uint64_t> previous(uint64_t node_id) const;
    size_t next_size(uint64_t node_id) const;
    size_t previous_size(uint64_t node_id) const;
    
    size_t path_size() const;
    const std::string& path_name(uint64_t path_id) const;
    std::vector<uint64_t> path(uint64_t path_id) const;
    uint64_t path_id(const std::string& name) const;
    
private:
    
    const SequenceGraph* seq_graph = nullptr;
    std::vector<size_t> cumul_len;
    std::vector<uint64_t> origin;
};

class BaseGraph {
public:
    
    BaseGraph() noexcept = default;
    BaseGraph(const BaseGraph& other) noexcept = default;
    BaseGraph(BaseGraph&& other) noexcept = default;
    ~BaseGraph() noexcept = default;
    
    BaseGraph& operator=(const BaseGraph& other) noexcept = default;
    BaseGraph& operator=(BaseGraph&& other) noexcept = default;
    
    uint64_t add_node(char base);
    void add_edge(uint64_t node_id_from, uint64_t node_id_to);
    void remove_edge(uint64_t node_id_from, uint64_t node_id_to);
    void relabel(uint64_t node_id, char base);
    
    size_t node_size() const;
    char label(uint64_t node_id) const;
    size_t label_size(uint64_t node_id) const;
    const std::vector<uint64_t>& next(uint64_t node_id) const;
    const std::vector<uint64_t>& previous(uint64_t node_id) const;
    size_t next_size(uint64_t node_id) const;
    size_t previous_size(uint64_t node_id) const;
    
    uint64_t add_path(const std::string& name);
    void extend_path(uint64_t path_id, uint64_t node_id);
    void pre_extend_path(uint64_t path_id, uint64_t node_id);
    
    size_t path_size() const;
    const std::string& path_name(uint64_t path_id) const;
    const std::vector<uint64_t>& path(uint64_t path_id) const;
    uint64_t path_id(const std::string& name) const;
    
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
    std::unordered_map<std::string, uint64_t> name_to_id;
    std::vector<std::pair<std::string, std::vector<uint64_t>>> paths;
};

}

#endif /* centrolign_graph_hpp */
