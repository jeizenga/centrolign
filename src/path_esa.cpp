#include <numeric>
#include <limits>
#include <unordered_set>

#include "centrolign/logging.hpp"
#include "centrolign/path_esa.hpp"
#include "centrolign/range_unique_query.hpp"

namespace centrolign {

using namespace std;

void PathESA::minimal_rare_matches(size_t max_count) const {

    logging::log(logging::Debug, "Constructing Range-Unique-Query structures");
    
    // construct range unique queries to compute subtree counts
    vector<RUQ> ruqs;
    ruqs.reserve(component_ranked_ids.size());
    for (const auto& ranked_ids : component_ranked_ids) {
        ruqs.emplace_back(ranked_ids);
    }
    
}

}
