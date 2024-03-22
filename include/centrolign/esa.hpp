#ifndef centrolign_esa_hpp
#define centrolign_esa_hpp

#include <vector>
#include <cstdint>
#include <algorithm>
#include <unordered_set>
#include <limits>

#include "centrolign/utility.hpp"
#include "centrolign/modify_graph.hpp"
#include "centrolign/logging.hpp"
#include "centrolign/range_unique_query.hpp"

namespace centrolign {


/*
 * A conceptual node in the suffix tree
 */
struct SANode {
    SANode(size_t begin, size_t end) : begin(begin), end(end) { }
    SANode() = default;
    ~SANode() = default;
    inline bool is_leaf() const;
    inline bool operator<(const SANode& other) const;
    inline bool operator==(const SANode& other) const;
    inline bool operator!=(const SANode& other) const;
    
    size_t begin = -1;
    size_t end = -1;
};

/*
 * Shared components for GESA and PathESA
 */
class ESA {
public:
    
    ESA() = default;
    ~ESA() = default;
    
    // the number of components used to construct
    inline size_t component_size() const;
    
    // the length of the shared prefix of an internal node, or for a leaf the length
    // of its minimal unique prefix
    inline size_t depth(const SANode& node) const;
    
protected:
    
    static const bool debug_esa;
    
    template<class LabelGetter>
    std::vector<std::tuple<SANode, size_t, std::vector<uint64_t>>>
    minimal_rare_matches_internal(size_t max_count, const LabelGetter& label_getter) const;
    
    template<class RUQType, class LabelGetter>
    std::vector<std::tuple<SANode, size_t, std::vector<uint64_t>>>
    minimal_rare_matches_internal_query(const std::vector<RUQType>& ruqs,
                                        size_t max_count, const LabelGetter& label_getter) const;
    
    template<class Advancer>
    std::vector<std::pair<size_t, std::vector<uint64_t>>> walk_matches_internal(const SANode& node, size_t length,
                                                                                const Advancer& advancer) const;
    
    void construct_child_array();
    
    template<class Advancer>
    void construct_suffix_links(const Advancer& advancer);
    
    inline SANode root() const;
    inline size_t first_l_index(const SANode& node) const;
    
    inline bool child_array_is_up(size_t i) const;
    inline bool child_array_is_down(size_t i) const;
    inline bool child_array_is_l_index(size_t i) const;
    
    std::vector<SANode> children(const SANode& parent) const;
    
    inline SANode link(const SANode& node) const;
    
    inline std::pair<size_t, size_t> st_node_annotation_idx(const SANode& node) const;
    
    inline size_t esa_struct_size() const;
    
    std::vector<size_t> child_array;
    std::vector<size_t> lcp_array;
    // component origin of each leaf node
    std::vector<uint16_t> leaf_to_comp;
    // suffix links of internal nodes
    std::vector<SANode> suffix_links;
    // for each component, the original node IDs of path nodes, in order of prefix rank
    std::vector<std::vector<uint64_t>> component_ranked_ids;
    // for each component, for each path node, the lowest within-component rank that is at
    // or before the path node
    std::vector<std::vector<size_t>> nearest_comp_rank;
    
};

/*
 * Template implementations
 */


bool SANode::is_leaf() const {
    return begin == end;
}

bool SANode::operator<(const SANode& other) const {
    return begin < other.begin;
}

bool SANode::operator==(const SANode& other) const {
    return begin == other.begin && end == other.end;
}

bool SANode::operator!=(const SANode& other) const {
    return !(*this == other);
}

size_t ESA::component_size() const {
    return component_ranked_ids.size();
}

inline SANode ESA::root() const {
    return SANode(0, lcp_array.size() - 1);
}

inline size_t ESA::first_l_index(const SANode& node) const {
    if (node == root()) {
        // a special case
        // TODO: which we could maybe get around by storing a special .up value at the end?
        return child_array.front();
    }
    else if (child_array_is_down(node.begin) && child_array[node.begin] <= node.end) {
        // TODO: i'm not really sure why we also need the <= end check, although this
        // agrees with the paper...
        
        // we prefer the down value, because the up value could belong to an ancestor
        // interval
        return child_array[node.begin];
    }
    else {
        return child_array[node.end];
    }
}

inline std::pair<size_t, size_t> ESA::st_node_annotation_idx(const SANode& node) const {
    if (node.is_leaf()) {
        return std::pair<size_t, size_t>(1, node.begin);
    }
    else {
        return std::pair<size_t, size_t>(0, first_l_index(node));
    }
}

inline bool ESA::child_array_is_down(size_t i) const {
    return i < child_array.size() ? (child_array[i] > i && lcp_array[child_array[i]] != lcp_array[i]) : false;
}

inline bool ESA::child_array_is_up(size_t i) const {
    // other pointers are all downward
    return i < child_array.size() ? child_array[i] <= i : false;
}

inline bool ESA::child_array_is_l_index(size_t i) const {
    return i < child_array.size() ? (child_array[i] > i && lcp_array[child_array[i]] == lcp_array[i]) : false;
}

inline size_t ESA::depth(const SANode& node) const {
    if (node.is_leaf()) {
        // longer of the two adjacent matches
        size_t l = lcp_array[node.begin];
        if (node.begin + 1 < lcp_array.size()) {
            l = std::max(l, lcp_array[node.begin + 1]);
        }
        // add 1 to make the prefix unique
        return l + 1;
    }
    else {
        return lcp_array[first_l_index(node)];
    }
}

inline SANode ESA::link(const SANode& node) const {
    size_t i, j;
    std::tie(i, j) = st_node_annotation_idx(node);
    return suffix_links[j];
}

template<class LabelGetter>
std::vector<std::tuple<SANode, size_t, std::vector<uint64_t>>>
ESA::minimal_rare_matches_internal(size_t max_count, const LabelGetter& label_getter) const {
    
    if (debug_esa) {
        std::cerr << "finding minimal rare matches with max count " << max_count << '\n';
    }
    
    logging::log(logging::Debug, "Constructing Range-Unique-Query structures");
    
    // the max size of the IDs vector across components
    size_t max_ids_size = 0;
    for (const auto& ranked_ids : component_ranked_ids) {
        max_ids_size = std::max(max_ids_size, ranked_ids.size());
    }
    
    std::vector<std::tuple<SANode, size_t, std::vector<uint64_t>>> matches;
    
    // macro to generate a templated code block
    #define _gen_match_query(UIntSize, NAry, Sampling) \
        /* construct range unique queries to compute subtree counts */ \
        size_t ruq_mem_size = 0; \
        std::vector<RUQ<UIntSize, NAry, Sampling>> ruqs; \
        ruqs.reserve(component_ranked_ids.size()); \
        for (const auto& ranked_ids : component_ranked_ids) { \
            ruqs.emplace_back(ranked_ids); \
            if (logging::level >= logging::Debug) { \
                ruq_mem_size += ruqs.back().memory_size(); \
            } \
        } \
        if (logging::level >= logging::Debug) { \
            logging::log(logging::Debug, std::string("Range-Unique-Query structures (") + #NAry + "-ary, "  + #Sampling + "-sampled) are occupying " + format_memory_usage(ruq_mem_size) + "."); \
            logging::log(logging::Debug, "Current memory usage is " + format_memory_usage(current_memory_usage()) + "."); \
        } \
        matches = std::move(minimal_rare_matches_internal_query(ruqs, max_count, label_getter))
    
    // trade off execution time for memory thrift with larger inputs
    if (max_ids_size < (1ul << 27)) {
        _gen_match_query(uint32_t, 3, 1);
    }
    else if (max_ids_size < (1ul << 30)) {
        _gen_match_query(uint32_t, 4, 2);
    }
    else if (max_ids_size < (1ul << 33)) {
        // 64-bit ints are now required 
        _gen_match_query(uint64_t, 6, 4);
    }
    else {
        _gen_match_query(uint64_t, 7, 6);
    }
    
    #undef _gen_match_query
    
    return matches;
}

template<class RUQType, class LabelGetter>
std::vector<std::tuple<SANode, size_t, std::vector<uint64_t>>>
ESA::minimal_rare_matches_internal_query(const std::vector<RUQType>& ruqs,
                                         size_t max_count, const LabelGetter& label_getter) const {
    
    std::vector<std::tuple<SANode, size_t, std::vector<uint64_t>>> matches;
    auto add_matches = [&](const SANode& parent, const std::vector<SANode>& children,
                           const std::vector<bool>& children_too_frequent) {
        
        // in order to have the same count as the children, we need one more character
        // after the parent's depth
        
        bool any_too_frequent = false;
        size_t unique_length = depth(parent) + 1;
        if (debug_esa) {
            std::cerr << "checking non-leaf children of " << parent.begin << ',' << parent.end << ", which has depth " << depth(parent) << '\n';
        }
        
        if (unique_length == 1) {
            if (debug_esa) {
                std::cerr << "parent is the root, handling as a special case\n";
            }
            // we only need to check for the max count on children of root
            for (size_t i = 0; i < children.size(); ++i) {
                if (children_too_frequent[i]) {
                    if (debug_esa) {
                        std::cerr << "skipping child " << children[i].begin << ',' << children[i].end << " as too frequent" << '\n';
                    }
                    any_too_frequent = true;
                    continue;
                }
                const auto& child = children[i];
                if (debug_esa) {
                    std::cerr << "considering match node " << child.begin << ',' << child.end << '\n';
                }
                std::vector<uint64_t> counts(component_size(), 0);
                for (size_t c = 0; c < component_size(); ++c) {
                    uint64_t count = ruqs[c].range_unique(nearest_comp_rank[c][child.begin],
                                                          nearest_comp_rank[c][child.end + 1]);
                    if (count == 0) {
                        break;
                    }
                    counts[c] = count;
                }
                uint64_t total_count = 1;
                for (auto c : counts) {
                    total_count = sat_mult(total_count, c);
                }
                
                if (total_count > 0 && total_count <= max_count) {
                    if (debug_esa) {
                        std::cerr << "is a minimal match with length " << unique_length << " and counts:";
                        for (auto cnt : counts) {
                            std::cerr << ' ' << cnt;
                        }
                        std::cerr << '\n';
                    }
                    matches.emplace_back(child, unique_length, std::move(counts));
                }
                else {
                    if (debug_esa) {
                        std::cerr << "not a minimal rare match\n";
                    }
                    any_too_frequent = true;
                }
            }
            return any_too_frequent;
        }
        
        // align the children of the parent's suffix link to the children we're testing
        auto suf_link = link(parent);
        std::vector<SANode> link_children = this->children(suf_link);

        for (size_t i = 0, j = 0; i < children.size(); ++j) {
            if (label_getter(parent, children[i]) == label_getter(suf_link, link_children[j])) {
                link_children[i] = link_children[j];
                ++i;
            }
        }
        link_children.resize(children.size());
        
        for (size_t k = 0; k < children.size(); ++k) {
            
            if (children_too_frequent[k]) {
                if (debug_esa) {
                    std::cerr << "skipping child " << children[k].begin << ',' << children[k].end << " as too frequent" << '\n';
                }
                any_too_frequent = true;
                continue;
            }

            if (debug_esa) {
                std::cerr << "considering match node " << children[k].begin << ',' << children[k].end << '\n';
            }
            
            auto& child = children[k];
            auto& link_child = link_children[k];
            
            std::vector<uint64_t> counts(component_size(), 0);
            bool link_more_frequent = false;
            for (size_t c = 0; c < component_size(); ++c) {
                uint64_t count = ruqs[c].range_unique(nearest_comp_rank[c][child.begin],
                                                      nearest_comp_rank[c][child.end + 1]);
                if (count == 0) {
                    break;
                }
                counts[c] = count;
                uint64_t link_count = ruqs[c].range_unique(nearest_comp_rank[c][link_child.begin],
                                                           nearest_comp_rank[c][link_child.end + 1]);
                link_more_frequent = (link_more_frequent || count < link_count);
            }
            uint64_t total_count = 1;
            for (auto c : counts) {
                total_count = sat_mult(total_count, c);
            }
            
            if (total_count > 0 && total_count <= max_count && link_more_frequent) {
                // occurs on more than one component and removing the first character
                // involves introducing more matches
                if (debug_esa) {
                    std::cerr << "is a minimal match with length " << unique_length << " and counts:";
                    for (auto cnt : counts) {
                        std::cerr << ' ' << cnt;
                    }
                    std::cerr << '\n';
                }
                matches.emplace_back(children[k], unique_length, std::move(counts));
            }
            else if (total_count > max_count) {
                
                any_too_frequent = true;
                
                if (debug_esa) {
                    std::cerr << "total count is too high\n";
                }
            }
            else if (debug_esa && total_count == 0) {
                std::cerr << "missing on at least one component\n";
            }
            else if (debug_esa && !link_more_frequent) {
                std::cerr << "suffix link has same occurrences\n";
            }
        }
        
        return any_too_frequent;
    };
    
    
    logging::log(logging::Debug, "Traversing the LCP tree");
    
    // records of (lcp, left, children, children too frequent)
    std::vector<std::tuple<size_t, size_t, std::vector<SANode>, std::vector<bool>>> stack;
    stack.emplace_back(0, 0, std::vector<SANode>(), false);
    for (size_t i = 1; i < lcp_array.size(); ++i) {
        
        SANode last_node(-1, -1);
        bool has_child_too_freq = false;
        
        // figure out which internal nodes we're leaving
        size_t left = i - 1;
        while (std::get<0>(stack.back()) > lcp_array[i]) {
            
            auto& top = stack.back();
            
            // emit an internal node
            last_node = SANode(std::get<1>(top), i - 1);
            has_child_too_freq = add_matches(last_node, std::get<2>(top), std::get<3>(top));
            
            left = std::get<1>(top);
            stack.pop_back();
            if (std::get<0>(stack.back()) >= lcp_array[i]) {
                // record this as a child of the parent
                std::get<2>(stack.back()).push_back(last_node);
                std::get<3>(stack.back()).push_back(has_child_too_freq);
                last_node = SANode(-1, -1);
                has_child_too_freq = false;
            }
        }
        if (std::get<0>(stack.back()) < lcp_array[i]) {
            stack.emplace_back(lcp_array[i], left, std::vector<SANode>(), false);
            if (last_node != SANode(-1, -1)) {
                std::get<2>(stack.back()).push_back(last_node);
                std::get<3>(stack.back()).push_back(has_child_too_freq);
            }
        }
    }
    // clear the stack
    while (!stack.empty()) {
        auto& top = stack.back();
        
        // emit an internal node
        SANode node(std::get<1>(top), lcp_array.size() - 1);
        bool has_child_too_freq = add_matches(node, std::get<2>(top), std::get<3>(top));
        stack.pop_back();
        
        // record this as a child of the parent
        if (!stack.empty()) {
            std::get<2>(stack.back()).push_back(node);
            std::get<3>(stack.back()).push_back(has_child_too_freq);
        }
    }
    
    return matches;
}


template<class Advancer>
void ESA::construct_suffix_links(const Advancer& advancer) {
    
    if (debug_esa) {
        std::cerr << "computing suffix links\n";
    }
    
    // TODO: decide if I want the leaf suffix links
    
    // records of (lcp, left, right)
    std::vector<std::tuple<size_t, size_t, size_t>> stack;
    
    // do a top down traversal to construct lists of l-intervals
    std::vector<std::vector<SANode>> l_interval_lists;
    stack.emplace_back(0, 0, -1);
    for (size_t i = 1; i < lcp_array.size(); ++i) {
        if (debug_esa) {
            std::cerr << "iteration " << i << ", stack state:\n";
            for (auto x : stack) {
                std::cerr << " (" << std::get<0>(x) << ',' << std::get<1>(x) << ',' << (int64_t) std::get<2>(x) << ')';
            }
            std::cerr << '\n';
        }
        
        // figure out which internal nodes we're leaving
        size_t left = i - 1;
        while (std::get<0>(stack.back()) > lcp_array[i]) {
            auto& top = stack.back();
            std::get<2>(top) = i - 1;
            if (debug_esa) {
                std::cerr << "emit internal [" << std::get<1>(top) << ", " << std::get<2>(top) << "], l = " << std::get<0>(top) << '\n';
            }
            while (l_interval_lists.size() <= std::get<0>(top)) {
                l_interval_lists.emplace_back();
            }
            l_interval_lists[std::get<0>(top)].emplace_back(std::get<1>(top), std::get<2>(top));
            left = std::get<1>(top);
            stack.pop_back();
        }
        if (lcp_array[i] > std::get<0>(stack.back())) {
            stack.emplace_back(lcp_array[i], left, -1);
        }
    }
    if (debug_esa) {
        std::cerr << "clearing the stack\n";
    }
    // clear the stack
    while (!stack.empty()) {
        auto& top = stack.back();
        std::get<2>(top) = lcp_array.size() - 1;
        if (debug_esa) {
            std::cerr << "emit internal [" << std::get<1>(top) << ", " << std::get<2>(top) << "], l = " << std::get<0>(top) << '\n';
        }
        while (l_interval_lists.size() <= std::get<0>(top)) {
            l_interval_lists.emplace_back();
        }
        l_interval_lists[std::get<0>(top)].emplace_back(std::get<1>(top), std::get<2>(top));
        stack.pop_back();
    }
    
    if (debug_esa) {
        std::cerr << "l interval lists:\n";
        for (size_t l = 0; l < l_interval_lists.size(); ++l) {
            std::cerr << l << ':';
            for (auto node : l_interval_lists[l]) {
                std::cerr << " [" << node.begin << ',' << node.end << "]";
            }
            std::cerr << '\n';
        }
    }
    
    suffix_links.resize(lcp_array.size());
    for (size_t l = 1; l < l_interval_lists.size(); ++l) {
        if (debug_esa) {
            std::cerr << "forming links for nodes at depth " << l << "\n";
        }
        auto& link_interval_list = l_interval_lists[l - 1];
        for (const SANode& node : l_interval_lists[l]) {
            
            size_t i, j;
            std::tie(i, j) = st_node_annotation_idx(node);
            
            // either SA^-1[SA[i] + 1] or use a PathGraph edge
            // use -1 as an end sentinel for the final value (TODO: is this necessary?)
            size_t next_rank = advancer(node.begin);
            if (next_rank == -1) {
                // special case, no successor
                if (debug_esa) {
                    std::cerr << "link [" << node.begin << ',' << node.end << "] -> [" << root().begin << ',' << root().end << "], storing at " << j << "\n";
                }
                suffix_links[j] = root();
            }
            else {
                // binary search to find the containing interval
                size_t lo = 0, hi = link_interval_list.size() - 1;
                while (lo != hi) {
                    size_t mid = (lo + hi) / 2;
                    if (next_rank < link_interval_list[mid].begin) {
                        hi = mid - 1;
                    }
                    else if (next_rank > link_interval_list[mid].end) {
                        lo = mid + 1;
                    }
                    else {
                        hi = lo = mid;
                    }
                }
                if (debug_esa) {
                    std::cerr << "link [" << node.begin << ',' << node.end << "] -> [" << link_interval_list[lo].begin << ',' << link_interval_list[lo].end << "], storing at " << j << "\n";
                }
                
                // record the suffix link
                suffix_links[j] = link_interval_list[lo];
            }

        }
    }
}

template<class Advancer>
std::vector<std::pair<size_t, std::vector<uint64_t>>> ESA::walk_matches_internal(const SANode& node, size_t length,
                                                                                 const Advancer& advancer) const {
    
    if (debug_esa) {
        std::cerr << "walking node " << node.begin << " " << node.end << '\n';
    }
    
    std::vector<std::pair<size_t, std::vector<uint64_t>>> matches;
    
    std::unordered_set<std::pair<size_t, uint64_t>> match_starts;
    for (size_t i = node.begin; i <= node.end; ++i) {
        
        // get the first node and the component information
        size_t idx = i;
        size_t comp = leaf_to_comp[idx];
        const auto& ranked_ids = component_ranked_ids[comp];
        const auto& nearest_rank = nearest_comp_rank[comp];
        uint64_t node_id = ranked_ids[nearest_rank[idx]];
        
        if (match_starts.count(std::make_pair(comp, node_id))) {
            // this is a duplicated node
            if (debug_esa) {
                std::cerr << "walk at " << i << " has duplicated start " << comp << " " << node_id << '\n';
            }
            continue;
        }
        
        match_starts.emplace(comp, node_id);
        
        matches.emplace_back();
        auto& match = matches.back();
        match.second.reserve(length);
        match.first = comp;
        match.second.push_back(node_id);
        
        // walk the rest of the match
        // FIXME: how should i choose from the combinatorially many identical matches
        // that there might be?
        for (size_t j = 1; j < length; ++j) {
            idx = advancer(idx);
            match.second.push_back(ranked_ids[nearest_rank[idx]]);
        }
        if (debug_esa) {
            std::cerr << "walk from " << i << ", comp " << match.first << ":\n";
            for (auto v : match.second) {
                std::cerr << '\t' << v << '\n';
            }
        }
    }
    
    return matches;
}

inline size_t ESA::esa_struct_size() const {
    return (sizeof(child_array) + child_array.capacity() * sizeof(decltype(child_array)::value_type)
            + sizeof(lcp_array) + lcp_array.capacity() * sizeof(decltype(lcp_array)::value_type)
            + sizeof(leaf_to_comp) + leaf_to_comp.capacity() * sizeof(decltype(leaf_to_comp)::value_type)
            + sizeof(suffix_links) + suffix_links.capacity() * sizeof(decltype(suffix_links)::value_type)
            + sizeof(component_ranked_ids) + component_ranked_ids.capacity() * sizeof(decltype(component_ranked_ids)::value_type)
            + component_ranked_ids.size() * (component_ranked_ids[0].capacity() * sizeof(decltype(component_ranked_ids)::value_type::value_type))
            + sizeof(nearest_comp_rank) + nearest_comp_rank.capacity() * sizeof(decltype(nearest_comp_rank)::value_type)
            + nearest_comp_rank.size() * (nearest_comp_rank[0].capacity() * sizeof(decltype(nearest_comp_rank)::value_type::value_type)));
}

}

#endif /* centrolign_esa_hpp */
