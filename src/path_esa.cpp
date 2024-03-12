#include "centrolign/path_esa.hpp"

namespace centrolign {

using namespace std;

vector<tuple<SANode, size_t, vector<uint64_t>>> PathESA::minimal_rare_matches(size_t max_count) const {
    auto label_getter = [&](const SANode& parent, const SANode& child) {
        return joined_seq[suffix_array[child.begin] + depth(parent)];
    };
    return minimal_rare_matches_internal(max_count, label_getter);
}

vector<pair<size_t, vector<uint64_t>>> PathESA::walk_matches(const SANode& node, size_t length) const {
    auto advance = [&](size_t i) -> size_t {
        return next(i);
    };
    return ESA::walk_matches_internal(node, length, advance);
}

size_t PathESA::memory_size() const {
    return (esa_struct_size()
            + sizeof(suffix_array) + suffix_array.capacity() * sizeof(decltype(suffix_array)::value_type)
            + sizeof(inverse_suffix_array) + inverse_suffix_array.capacity() * sizeof(decltype(inverse_suffix_array)::value_type)
            + sizeof(joined_seq) + joined_seq.capacity() * sizeof(decltype(joined_seq)::value_type));
}

void PathESA::print(ostream& out) const {
    out << "i" << '\t' << "SA" << '\t' << "Cmp" << '\t' << "Nd" << '\t' << "LCP" << '\t' << "Ch" << '\t' << "SL1" << '\t' << "SL2" << '\n';
    for (size_t i = 0; i < lcp_array.size(); ++i) {
        const auto& ranked_id = component_ranked_ids[leaf_to_comp[i]];
        const auto& nearest_rank = nearest_comp_rank[leaf_to_comp[i]];
        out << i << '\t' << suffix_array[i] << '\t' << leaf_to_comp[i] << '\t' << ranked_id[nearest_rank[i]] << '\t' << lcp_array[i] << '\t' << (i < child_array.size() ? (int) child_array[i] : -1);
        out << " (";
        if (child_array_is_down(i)) {
            out << 'D';
        }
        else if (child_array_is_l_index(i)) {
            out << 'N';
        }
        else if (child_array_is_up(i)) {
            out << 'U';
        }
        else {
            out << '.';
        }
        out << ')';
        if (suffix_links[i] != SANode()) {
            out << '\t' << suffix_links[i].begin << '\t' << suffix_links[i].end;
        }
        else {
            out << '\t' << -1 << '\t' << -1;
        }
        out << '\n';
    }
}

}
