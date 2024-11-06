#ifndef centrolign_packed_match_bank_hpp
#define centrolign_packed_match_bank_hpp

#include <vector>
#include <tuple>
#include <limits>
#include <utility>

#include "centrolign/match_finder.hpp"
#include "centrolign/packed_vector.hpp"
#include "centrolign/paged_vector.hpp"

namespace centrolign {

/*
 * A support structure for doing dynamic programming over a match set using bit-packed
 * underlying data structures
 */
template<typename UIntAnchor>
class PackedMatchBank {
public:
    
    template<class BGraph>
    PackedMatchBank(const BGraph& graph1, const std::vector<match_set_t>& matches, size_t num_match_sets, bool suppress_verbose_logging,
              const std::unordered_set<std::tuple<size_t, size_t, size_t>>* masked_matches = nullptr) noexcept;
    PackedMatchBank() noexcept = default;
    PackedMatchBank(const PackedMatchBank<UIntAnchor>& other) noexcept = default;
    PackedMatchBank(PackedMatchBank<UIntAnchor>&& other) noexcept = default;
    ~PackedMatchBank() noexcept = default;
    
    PackedMatchBank& operator=(const PackedMatchBank<UIntAnchor>& other) noexcept = default;
    PackedMatchBank& operator=(PackedMatchBank<UIntAnchor>&& other) noexcept = default;
    
    using match_id_t = UIntAnchor;
    
    inline std::tuple<size_t, size_t, size_t> get_match_indexes(const match_id_t& match_id) const;
    
    inline const std::vector<uint64_t>& walk1(const match_id_t& match_id) const;
    inline const std::vector<uint64_t>& walk2(const match_id_t& match_id) const;
    inline const match_set_t& match_set(const match_id_t& match_id) const;
    
    inline double dp_value(const match_id_t& match_id) const;
    inline match_id_t backpointer(const match_id_t& match_id) const;
    inline void update_dp(const match_id_t& match_id, double value, const match_id_t& traceback_id);
    
    std::vector<match_id_t> starts_on(uint64_t node_id) const;
    std::vector<match_id_t> ends_on(uint64_t node_id) const;
    
    inline static match_id_t max();
    inline static match_id_t min();
    
    
    class iterator {
    public:
        iterator() = default;
        ~iterator() = default;
        
        iterator& operator++();
        const match_id_t& operator*() const;
        const match_id_t* operator->() const;
        bool operator==(const iterator& other) const;
        bool operator!=(const iterator& other) const;
        
        
    private:
        friend class PackedMatchBank<UIntAnchor>;
        // internal constructor
        iterator(const PackedMatchBank<UIntAnchor>& iteratee, UIntAnchor i);
        
        const PackedMatchBank<UIntAnchor>* iteratee = nullptr;
        match_id_t match_id;
    };
    
    iterator begin() const;
    iterator end() const;
    
private:
    
    friend class iterator;
    
    const std::vector<match_set_t>* matches = nullptr;
    
    PagedVector<1028, 1027> set_idx;
    PagedVector<128, 15> walk1_idx;
    PackedVector walk2_idx;
    
    PagedVector<256, 31> start_interval;
    PagedVector<256, 31> end_interval;
    PackedVector start_match_ids;
    PackedVector end_match_ids;
    
    std::vector<double> dp;
    PackedVector backpointer_id;
    
//    std::vector<std::vector<std::vector<std::pair<double, match_id_t>>>> dp;
//    std::vector<std::vector<std::pair<UIntSet, UIntMatch>>> starts;
//    std::vector<std::vector<std::pair<UIntSet, UIntMatch>>> ends;
    
};

/*
 * Template implementations
 */


template<typename UIntAnchor>
template<class BGraph>
PackedMatchBank<UIntAnchor>::PackedMatchBank(const BGraph& graph1, const std::vector<match_set_t>& matches, size_t num_match_sets, bool suppress_verbose_logging, const std::unordered_set<std::tuple<size_t, size_t, size_t>>* masked_matches) noexcept : matches(&matches) {
        
    static const bool debug = false;
    
    size_t num_match_ids = 0;
    std::vector<UIntAnchor> num_starts(graph1.node_size() + 1, 0);
    std::vector<UIntAnchor> num_ends(graph1.node_size() + 1, 0);
    for (size_t i = 0; i < num_match_sets; ++i) {
        const auto& match_set = matches[i];
        for (size_t j = 0; j < match_set.walks1.size(); ++j) {
            const auto& walk1 = match_set.walks1[j];
            for (size_t k = 0; k < match_set.walks2.size(); ++k) {
                if (masked_matches && masked_matches->count(std::make_tuple(i, j, k))) {
                    continue;
                }
                ++num_match_ids;
                ++num_starts[walk1.front() + 1];
                ++num_ends[walk1.back() + 1];
            }
        }
    }
    
    
    dp.resize(num_match_ids, std::numeric_limits<typename decltype(dp)::value_type>::lowest());
    backpointer_id = std::move(decltype(backpointer_id)(num_match_ids));
    
    set_idx = std::move(decltype(set_idx)(num_match_ids));
    walk1_idx = std::move(decltype(walk1_idx)(num_match_ids));
    walk2_idx = std::move(decltype(walk2_idx)(num_match_ids));
    
    start_match_ids = std::move(decltype(start_match_ids)(num_match_ids));
    start_interval = std::move(decltype(start_interval)(graph1.node_size() + 1));
    
    end_match_ids = std::move(decltype(end_match_ids)(num_match_ids));
    end_interval = std::move(decltype(start_interval)(graph1.node_size() + 1));
    for (size_t i = 0; i < graph1.node_size(); ++i) {
        start_interval.set(i, num_starts[i]);
        end_interval.set(i, num_ends[i]);
        num_starts[i + 1] += num_starts[i];
        num_ends[i + 1] += num_ends[i];
    }
    start_interval.set(start_interval.size() - 1, num_starts.back());
    end_interval.set(end_interval.size() - 1, num_ends.back());
    
    size_t match_id = 0;
    for (size_t i = 0; i < num_match_sets; ++i) {
        const auto& match_set = matches[i];
        for (size_t j = 0; j < match_set.walks1.size(); ++j) {
            const auto& walk1 = match_set.walks1[j];
            for (size_t k = 0; k < match_set.walks2.size(); ++k) {
                if (masked_matches && masked_matches->count(std::make_tuple(i, j, k))) {
                    continue;
                }
                set_idx.set(match_id, i);
                walk1_idx.set(match_id, j);
                walk2_idx.set(match_id, k);
                start_match_ids.set(num_starts[walk1.front()]++, match_id);
                end_match_ids.set(num_ends[walk1.back()]++, match_id);
                ++match_id;
            }
        }
    }
    
    if (debug) {
        
        std::cerr << "ID translator\n";
        for (size_t i = 0; i < set_idx.size(); ++i ) {
            std::cerr << i << "\t" << set_idx.at(i) << '\t' << walk1_idx.at(i) << '\t' << walk2_idx.at(i) << '\n';
        }
        std::cerr << "start/end intervals\n";
        for (size_t i = 0; i < start_interval.size(); ++i) {
            std::cerr << i << '\t' << start_interval.at(i) << '\t' << end_interval.at(i) << '\n';
        }
        std::cerr << "start/end IDs\n";
        for (size_t i = 0; i < start_match_ids.size(); ++i) {
            std::cerr << i << '\t' << start_match_ids.at(i) << '\t' << end_match_ids.at(i) << '\n';
        }
    }
    
    if (!suppress_verbose_logging) {
        // dp table
        size_t dp_size = sizeof(dp) + dp.capacity() * sizeof(typename decltype(dp)::value_type) + backpointer_id.memory_size();
        
        // start and end locations
        size_t start_end_size = start_interval.memory_size() + start_match_ids.memory_size() + end_interval.memory_size() + end_match_ids.memory_size();
        
        // ID to match translation
        size_t id_size = set_idx.memory_size() + walk1_idx.memory_size() + walk2_idx.memory_size();
        
        logging::log(logging::Debug, "Dynamic programming table is occupying " + format_memory_usage(dp_size) + ".");
        logging::log(logging::Debug, "Match start and end locations are occupying " + format_memory_usage(start_end_size) + ".");
        logging::log(logging::Debug, "Match IDs are occupying " + format_memory_usage(id_size) + ".");
    }
//    std::cerr << "finish making packed bank\n";
}


template<typename UIntAnchor>
inline std::tuple<size_t, size_t, size_t> PackedMatchBank<UIntAnchor>::get_match_indexes(const match_id_t& match_id) const {
    return std::tuple<size_t, size_t, size_t>(set_idx.at(match_id), walk1_idx.at(match_id), walk2_idx.at(match_id));
}

template<typename UIntAnchor>
inline const std::vector<uint64_t>& PackedMatchBank<UIntAnchor>::walk1(const match_id_t& match_id) const {
    return (*matches)[set_idx.at(match_id)].walks1[walk1_idx.at(match_id)];
}

template<typename UIntAnchor>
inline const std::vector<uint64_t>& PackedMatchBank<UIntAnchor>::walk2(const match_id_t& match_id) const {
    return (*matches)[set_idx.at(match_id)].walks2[walk2_idx.at(match_id)];
}

template<typename UIntAnchor>
inline const match_set_t& PackedMatchBank<UIntAnchor>::match_set(const match_id_t& match_id) const {
    return (*matches)[set_idx.at(match_id)];
}

template<typename UIntAnchor>
inline double PackedMatchBank<UIntAnchor>::dp_value(const match_id_t& match_id) const {
    return dp[match_id];
}

template<typename UIntAnchor>
inline typename PackedMatchBank<UIntAnchor>::match_id_t PackedMatchBank<UIntAnchor>::backpointer(const match_id_t& match_id) const {
    auto bp = backpointer_id.at(match_id);
    return bp ? bp - 1 : max();
}

template<typename UIntAnchor>
inline void PackedMatchBank<UIntAnchor>::update_dp(const match_id_t& match_id, double value, const match_id_t& traceback_id) {
    static const bool debug = false;
    auto& dp_entry = dp[match_id];
    if (debug) {
        std::cerr << "checking match id " << match_id << " (" << set_idx.at(match_id) << ", " << walk1_idx.at(match_id) << ", " << walk2_idx.at(match_id) << ") with current value " << dp_entry << " against value " << value << '\n';
    }
    if (value > dp_entry) {
        dp_entry = value;
        // note: counting on -1's to overflow to 0
        backpointer_id.set(match_id, traceback_id + 1);
        if (debug) {
            std::cerr << "is new best\n";
        }
    }
}


template<typename UIntAnchor>
std::vector<typename PackedMatchBank<UIntAnchor>::match_id_t> PackedMatchBank<UIntAnchor>::starts_on(uint64_t node_id) const {
    std::vector<match_id_t> start_ids;
    size_t i = start_interval.at(node_id), n = start_interval.at(node_id + 1);
    start_ids.reserve(n - i);
    for (; i < n; ++i) {
        start_ids.push_back(start_match_ids.at(i));
    }
    return start_ids;
}

template<typename UIntAnchor>
std::vector<typename PackedMatchBank<UIntAnchor>::match_id_t> PackedMatchBank<UIntAnchor>::ends_on(uint64_t node_id) const {
    std::vector<match_id_t> end_ids;
    size_t i = end_interval.at(node_id), n = end_interval.at(node_id + 1);
    end_ids.reserve(n - i);
    for (; i < n; ++i) {
        end_ids.push_back(end_match_ids.at(i));
    }
    return end_ids;
}




template<typename UIntAnchor>
typename PackedMatchBank<UIntAnchor>::iterator PackedMatchBank<UIntAnchor>::begin() const {
    return iterator(*this, 0);
}

template<typename UIntAnchor>
typename PackedMatchBank<UIntAnchor>::iterator PackedMatchBank<UIntAnchor>::end() const {
    return iterator(*this, set_idx.size());
}


template<typename UIntAnchor>
inline typename PackedMatchBank<UIntAnchor>::match_id_t PackedMatchBank<UIntAnchor>::max() {
    return match_id_t(-1);
}

template<typename UIntAnchor>
inline typename PackedMatchBank<UIntAnchor>::match_id_t PackedMatchBank<UIntAnchor>::min() {
    return match_id_t(0);
}

template<typename UIntAnchor>
PackedMatchBank<UIntAnchor>::iterator::iterator(const PackedMatchBank<UIntAnchor>& iteratee, UIntAnchor i) : iteratee(&iteratee), match_id(i) {

}

template<typename UIntAnchor>
typename PackedMatchBank<UIntAnchor>::iterator& PackedMatchBank<UIntAnchor>::iterator::operator++() {
    ++match_id;
    return *this;
}

template<typename UIntAnchor>
const typename PackedMatchBank<UIntAnchor>::match_id_t& PackedMatchBank<UIntAnchor>::iterator::operator*() const {
    return match_id;
}

template<typename UIntAnchor>
const typename PackedMatchBank<UIntAnchor>::match_id_t* PackedMatchBank<UIntAnchor>::iterator::operator->() const {
    return &match_id;
}

template<typename UIntAnchor>
bool PackedMatchBank<UIntAnchor>::iterator::operator==(const iterator& other) const {
    return iteratee == other.iteratee && match_id == other.match_id;
}

template<typename UIntAnchor>
bool PackedMatchBank<UIntAnchor>::iterator::operator!=(const iterator& other) const {
    return !(*this == other);
}


}

#endif /* centrolign_packed_match_bank_hpp */
