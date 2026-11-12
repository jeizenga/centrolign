#ifndef centrolign_compressible_graph_hpp
#define centrolign_compressible_graph_hpp

#include <vector>
#include <string>
#include <atomic>
#include <cstdint>

#include "centrolign/graph.hpp"
#include "centrolign/modify_graph.hpp"
#include "centrolign/packed_vector.hpp"

namespace centrolign {

/*
 * A static, memory-efficient representation of a BaseGraph, in which non-branching paths
 * are compacted into sequence nodes. Sequence nodes are maximal except that each path must
 * start at the beginning of a sequence node and end at the end of one. Node IDs are not
 * preserved, but the topology, the labels, the paths, and the order of each node's next
 * edges are.
 */
class PackedSequenceGraph {
public:
    PackedSequenceGraph() noexcept = default;
    // all paths must follow edges in the graph, else an error is thrown
    // if provided, the translation is filled with the decompressed node ID of each original node ID
    PackedSequenceGraph(const BaseGraph& graph, const SentinelTableau& tableau,
                        std::vector<uint64_t>* translation = nullptr);
    PackedSequenceGraph(const PackedSequenceGraph& other) noexcept = default;
    PackedSequenceGraph(PackedSequenceGraph&& other) noexcept = default;
    ~PackedSequenceGraph() noexcept = default;

    PackedSequenceGraph& operator=(const PackedSequenceGraph& other) noexcept = default;
    PackedSequenceGraph& operator=(PackedSequenceGraph&& other) noexcept = default;

    // reconstruct the BaseGraph and update the tableau's sentinel node IDs to match it
    BaseGraph decompress(SentinelTableau& tableau) const;

    // the number of nodes in the decompressed graph
    size_t node_size() const;

    size_t memory_size() const;

private:

    struct PackedPath {
        // the index of the first sequence node, or -1 if the path is empty
        uint64_t start = -1;
        // the index of the next edge taken at the end of each sequence node but the last
        PackedVector steps;
    };

    // sentinels are the final two sequence nodes, and their labels are stored separately
    // so that they do not inflate the width of the labels
    bool has_sentinels = false;
    char src_label = 0;
    char snk_label = 0;

    // the node labels, concatenated in decompressed node ID order
    PackedVector labels;
    // the first decompressed node ID of each sequence node, with a final past-the-end value
    PackedVector node_starts;
    // the number of next edges from each sequence node
    PackedVector out_degrees;
    // the sequence node index of each edge's target, grouped by source sequence node
    PackedVector edge_targets;

    std::vector<std::string> path_names;
    std::vector<PackedPath> paths;
};

/*
 * A wrapper for a BaseGraph and its sentinel tableau that can switch to a memory-efficient
 * internal representation while the graph is not in use. Node IDs are not preserved by
 * compression.
 */
class CompressibleGraph {
public:
    CompressibleGraph() noexcept = default;
    CompressibleGraph(BaseGraph&& graph, const SentinelTableau& tableau) noexcept;
    CompressibleGraph(const CompressibleGraph& other) noexcept;
    CompressibleGraph(CompressibleGraph&& other) noexcept;
    ~CompressibleGraph() noexcept = default;

    CompressibleGraph& operator=(const CompressibleGraph& other) noexcept;
    CompressibleGraph& operator=(CompressibleGraph&& other) noexcept;

    // ensure that the internal representation is memory-efficient (invalidates references
    // to the BaseGraph and the tableau)
    void compress();

    // ensure that the internal representation is a BaseGraph and return it
    BaseGraph& get_base_graph();
    // ensure that the internal representation is a BaseGraph and return its tableau
    SentinelTableau& get_tableau();

    // these require that the graph not be compressed, else an error is thrown
    const BaseGraph& get_base_graph() const;
    const SentinelTableau& get_tableau() const;

    bool is_compressed() const;

    // valid in either representation
    size_t node_size() const;

    // the memory size as of the last change of representation
    size_t memory_size() const;

private:

    void decompress();

    BaseGraph graph;
    SentinelTableau tableau;
    PackedSequenceGraph packed_graph;
    bool compressed = false;

    // can be queried while another thread changes the representation
    std::atomic<size_t> cached_memory_size{0};
};

}

#endif /* centrolign_compressible_graph_hpp */
