#include "centrolign/match_finder.hpp"



namespace centrolign {

std::vector<match_set_t> GESAMatchFinder::index_and_query(ExpandedGraph& expanded1,
                                                          ExpandedGraph& expanded2) const {
    
    std::vector<match_set_t> matches;
    try {
        std::vector<const BaseGraph*> graph_ptrs{&expanded1.graph, &expanded2.graph};
        std::vector<const std::vector<uint64_t>*> trans_ptrs{&expanded1.back_translation, &expanded2.back_translation};
        std::vector<const SentinelTableau*> tableau_ptrs{&expanded1.tableau, &expanded2.tableau};
        
        // FIXME: this doesn't play well with the graphs being simplified recursively
        size_t size_limit = size_limit_factor * (expanded1.graph.node_size() + expanded2.graph.node_size());
        
        GESA gesa(graph_ptrs, trans_ptrs, size_limit);
        matches = std::move(query_index(gesa));
    }
    catch (GESASizeException& ex) {
        
        logging::log(logging::Verbose, "Graph not simple enough to index, resimplifying.");
        
        auto targets = simplifier.identify_target_nodes(ex.from_counts());
        
        size_t simplify_dist = (1 << ex.doubling_step());
        
        size_t pre_simplify_size1 = expanded1.graph.node_size();
        size_t pre_simplify_size2 = expanded2.graph.node_size();
        
        auto expanded_more1 = simplifier.targeted_simplify(expanded1.graph, expanded1.tableau,
                                                           targets[0], simplify_dist);
        auto expanded_more2 = simplifier.targeted_simplify(expanded2.graph, expanded2.tableau,
                                                           targets[1], simplify_dist);
        
        for (auto& tr : expanded_more1.back_translation) {
            tr = expanded1.back_translation[tr];
        }
        for (auto& tr : expanded_more2.back_translation) {
            tr = expanded2.back_translation[tr];
        }
        
        expanded1 = std::move(expanded_more1);
        expanded2 = std::move(expanded_more2);
        
        if (pre_simplify_size1 == expanded1.graph.node_size() &&
            pre_simplify_size2 == expanded2.graph.node_size()) {
            
            throw std::runtime_error("Simplification algorithm failed to simplify graph");
        }
        
        // recursively try again with a more simplified graph
        matches = std::move(index_and_query(expanded1, expanded2));
    }
    
    return matches;
}


}
