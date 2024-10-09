#ifndef centrolign_execution_hpp
#define centrolign_execution_hpp

#include <vector>
#include <string>
#include <tuple>
#include <functional>

#include "centrolign/graph.hpp"
#include "centrolign/modify_graph.hpp"
#include "centrolign/alignment.hpp"
#include "centrolign/tree.hpp"

namespace centrolign {

/*
 * An alignment subproblem in the process of making the full MSA
 */
struct Subproblem {
    Subproblem() noexcept = default;
    Subproblem(const Subproblem& other) noexcept = default;
    Subproblem(Subproblem&& other) noexcept = default;
    ~Subproblem() = default;
    Subproblem& operator=(const Subproblem& other) noexcept = default;
    Subproblem& operator=(Subproblem&& other) noexcept = default;
    
    BaseGraph graph;
    SentinelTableau tableau;
    Alignment alignment;
    std::string name;
    bool complete = false;
};

/*
 * Class that keeps track of the inputs, state, and execution order of the MSA
 */
class Execution {
public:
    
    Execution(std::vector<std::pair<std::string, std::string>>&& names_and_sequences,
              Tree&& tree);
    
    Execution() = default;
    Execution(const Execution& other) = default;
    Execution(Execution&& other) = default;
    ~Execution() = default;
    Execution& operator=(const Execution& other) = default;
    Execution& operator=(Execution&& other) = default;
    
    bool finished() const;
    // the next subproblem and its two children (parent first)
    std::tuple<Subproblem*, Subproblem*, Subproblem*> next();
    // the problem that is currently being executed (or nullptrs if not currently mid-execution)
    std::tuple<const Subproblem*, const Subproblem*, const Subproblem*> current() const;
    
    // restart from saved partial results
    void restart(std::function<std::string(const Subproblem&)>& file_location,
                 bool preserve_leaves, bool preserve_internal_nodes);
    
    // get the subproblem corresponding to one input sequence
    const Subproblem& leaf_subproblem(const std::string& name) const;
    
    // get the last subproblem, which contains the full alignment after execution
    Subproblem& final_subproblem();
    const Subproblem& final_subproblem() const;
    
    // get all of the subproblems that correspond to an input sequence
    std::vector<Subproblem*> leaf_subproblems();
    
    // get a hash identifier for a subproblem
    uint64_t subproblem_hash(const Subproblem& subproblem) const;
    
    // get the names of the sequences involved in a given subproblem
    std::vector<std::string> leaf_descendents(const Subproblem& subproblem) const;
    
    // TODO: i don't love exposing this, but it's currently needed to make the subpath
    // tree during cyclic graph polishing
    const Tree& get_tree() const;
    
private:
    
    uint64_t get_tree_id(const Subproblem& subproblem) const;
    
    // the guide tree
    Tree tree;
    
    // the individual alignment subproblems (including single-sequence leaves)
    std::vector<Subproblem> subproblems;
    
    // a postorder of the inner nodes of the tree
    std::vector<uint64_t> execution_order;
    
    // the next subproblem in the execution order that we need to do
    size_t next_subproblem = -1;
    
public:
    
    size_t memory_size() const;
};


}
#endif /* centrolign_execution_hpp */
