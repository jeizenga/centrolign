#ifndef centrolign_core_hpp
#define centrolign_core_hpp

#include <vector>
#include <string>

#include "centrolign/chain_merge.hpp"
#include "centrolign/modify_graph.hpp"
#include "centrolign/anchorer.hpp"
#include "centrolign/stitcher.hpp"
#include "centrolign/tree.hpp"
#include "centrolign/simplifier.hpp"
#include "centrolign/match_finder.hpp"

namespace centrolign {

/*
 * The core runner for centrolign
 */
class Core {
public:
    
    // parse files to construct core (either may be - for stdin)
    Core(const std::string& fasta_file, const std::string& tree_file);
    
    // construct core from already-parsed inputs (consumes the inputs)
    Core(std::vector<std::pair<std::string, std::string>>&& names_and_sequences,
         Tree&& tree);
    
    Core() = default;
    ~Core() = default;
    
    // trigger the MSA
    void execute();
    
    // configurable submodules
    Simplifier simplifier;
    MatchFinder match_finder;
    Anchorer anchorer;
    Stitcher stitcher;
    
    struct Subproblem {
        Subproblem() = default;
        ~Subproblem() = default;
        
        BaseGraph graph;
        SentinelTableau tableau;
        Alignment alignment;
        std::string name;
        bool complete = false;
    };
    
    // preserve subproblems whose parent problems have been completed
    bool preserve_subproblems = true;
    
    // if non-empty prefix to give GFA output for all suproblems
    std::string subproblems_prefix;
    
    // load alignments from the prefix and start where they left off
    void restart();
    
    // the root subproblem
    const Subproblem& root_subproblem() const;
    
    const Subproblem& leaf_subproblem(const std::string& name) const;
    
    // the narrowest subproblem that includes all these sequences
    const Subproblem& subproblem_covering(const std::vector<std::string>& names) const;
    
private:
    
    void init(std::vector<std::pair<std::string, std::string>>&& names_and_sequences,
              Tree&& tree_in);
    
    std::string subproblem_file_name(uint64_t tree_id) const;
    
    void emit_subproblem(uint64_t tree_id) const;
    
    std::vector<std::string> leaf_descendents(uint64_t tree_id) const;
    
    std::vector<match_set_t> find_matches(ExpandedGraph& expanded1,
                                          ExpandedGraph& expanded2) const;
    
    Tree tree;
    
    std::vector<Subproblem> subproblems;
    
};

}

#endif /* centrolign_core_hpp */
