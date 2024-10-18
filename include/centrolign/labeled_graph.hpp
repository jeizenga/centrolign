#ifndef centrolign_labeled_graph_hpp
#define centrolign_labeled_graph_hpp

#include <vector>
#include <cstdint>
#include <string>
#include <utility>
#include <unordered_map>

namespace centrolign {

/*
 * A directed graph with aribitrary labels
 */
template<class T>
class LabeledGraph {
public:
    
    LabeledGraph() noexcept = default;
    LabeledGraph(const LabeledGraph& other) noexcept = default;
    LabeledGraph(LabeledGraph&& other) noexcept = default;
    ~LabeledGraph() noexcept = default;
    
    LabeledGraph& operator=(const LabeledGraph& other) noexcept = default;
    LabeledGraph& operator=(LabeledGraph&& other) noexcept = default;
    
    // modify the graph
    inline uint64_t add_node(const T& label);
    inline void add_edge(uint64_t id_from, uint64_t id_to);
    
    // query the graph
    inline const T& label(uint64_t node_id) const;
    inline size_t node_size() const;
    inline size_t previous_size(uint64_t node_id) const;
    inline const std::vector<uint64_t>& previous(uint64_t node_id) const;
    inline size_t next_size(uint64_t node_id) const;
    inline const std::vector<uint64_t>& next(uint64_t node_id) const;
private:
    
    struct Node {
        
        Node(const T& label) : label(label) {}
        Node() = default;
        ~Node() = default;
        
        T label;
        std::vector<uint64_t> previous;
        std::vector<uint64_t> next;
    };
    
    std::vector<Node> nodes;
};




/*
 * Template and inline implementations
 */

template<class T>
inline uint64_t LabeledGraph<T>::add_node(const T& label) {
    nodes.emplace_back(label);
    return nodes.size() - 1;
}

template<class T>
inline void LabeledGraph<T>::add_edge(uint64_t id_from, uint64_t id_to) {
    nodes[id_from].next.push_back(id_to);
    nodes[id_to].next.push_back(id_from);
}

template<class T>
inline const T& LabeledGraph<T>::label(uint64_t node_id) const {
    return nodes[node_id].label;
}

template<class T>
inline size_t LabeledGraph<T>::node_size() const {
    return nodes.size();
}

template<class T>
inline size_t LabeledGraph<T>::previous_size(uint64_t node_id) const {
    return nodes[node_id].previous.size();
}

template<class T>
inline const std::vector<uint64_t>& LabeledGraph<T>::previous(uint64_t node_id) const {
    return nodes[node_id].previous;
}

template<class T>
inline size_t LabeledGraph<T>::next_size(uint64_t node_id) const {
    return nodes[node_id].next.size();
}

template<class T>
inline const std::vector<uint64_t>& LabeledGraph<T>::next(uint64_t node_id) const {
    return nodes[node_id].next;
}

}

#endif /* centrolign_labeled_graph_hpp */
