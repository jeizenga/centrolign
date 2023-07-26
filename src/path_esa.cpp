#include <numeric>
#include <limits>
#include <unordered_set>

#include "centrolign/logging.hpp"
#include "centrolign/path_esa.hpp"
#include "centrolign/range_unique_query.hpp"

namespace centrolign {

using namespace std;

vector<tuple<SANode, size_t, vector<uint64_t>>> PathESA::minimal_rare_matches(size_t max_count) const {
    auto label_getter = [&](const SANode& node) {
        return joined_seq[suffix_array[node.begin] + depth(node)];
    };
    return minimal_rare_matches_internal(max_count, label_getter);
}

}
