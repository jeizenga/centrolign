#ifndef centrolign_core_hpp
#define centrolign_core_hpp

#include <vector>
#include <string>

#include "centrolign/chain_merge.hpp"
#include "centrolign/modify_graph.hpp"
#include "centrolign/anchorer.hpp"
#include "centrolign/stitcher.hpp"
#include "centrolign/tree.hpp"

namespace centrolign {

/*
 * The core runner for centrolign
 */
class Core {
public:
    
    Core(const std::string& fasta_file, const std::string& tree_file);
    Core() = default;
    ~Core() = default;
    
    Anchorer anchorer;
    Stitcher stitcher;
    
    // preserve subproblems whose parent problems have been completed
    bool preserve_subproblems = true;
    
private:
    
    struct Subproblem {
        Subproblem() = default;
        ~Subproblem() = default;
        
        BaseGraph graph;
        ChainMerge chain_merge;
        SentinelTableau tableau;
        std::string name;
    };
    
    Tree tree;
    
    std::vector<Subproblem> subproblems;
    
};

}

#endif /* centrolign_core_hpp */
