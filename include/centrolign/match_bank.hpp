#ifndef centrolign_match_bank_hpp
#define centrolign_match_bank_hpp

#include <vector>
#include <tuple>
#include <limits>
#include <utility>

#include "centrolign/match_finder.hpp"

namespace centrolign {

/*
 * A support structure for doing dynamic programming over a match set
 */
template<typename UIntSet = size_t, typename UIntMatch = size_t, typename ScoreFloat = double>
class MatchBank {
public:
    
    template<class BGraph>
    MatchBank(const BGraph& graph1, const std::vector<match_set_t>& matches, size_t num_match_sets, bool suppress_verbose_logging,
              const std::unordered_set<std::tuple<size_t, size_t, size_t>>* masked_matches = nullptr) noexcept;
    MatchBank() noexcept = default;
    MatchBank(const MatchBank<UIntSet, UIntMatch, ScoreFloat>& other) noexcept = default;
    MatchBank(MatchBank<UIntSet, UIntMatch, ScoreFloat>&& other) noexcept = default;
    ~MatchBank() noexcept = default;
    
    MatchBank& operator=(const MatchBank<UIntSet, UIntMatch, ScoreFloat>& other) noexcept = default;
    MatchBank& operator=(MatchBank<UIntSet, UIntMatch, ScoreFloat>&& other) noexcept = default;
    
    using match_id_t = std::tuple<UIntSet, UIntMatch, UIntMatch>;
    
    inline std::tuple<UIntSet, UIntMatch, UIntMatch> get_match_indexes(const match_id_t& match_id) const;
    
    inline const std::vector<uint64_t>& walk1(const match_id_t& match_id) const;
    inline const std::vector<uint64_t>& walk2(const match_id_t& match_id) const;
    inline const match_set_t& match_set(const match_id_t& match_id) const;
    
    inline ScoreFloat dp_value(const match_id_t& match_id) const;
    inline const match_id_t& backpointer(const match_id_t& match_id) const;
    inline void update_dp(const match_id_t& match_id, ScoreFloat value, const match_id_t& traceback_id);
    
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
        friend class MatchBank<UIntSet, UIntMatch, ScoreFloat>;
        // internal constructor
        iterator(const MatchBank<UIntSet, UIntMatch, ScoreFloat>& iteratee, UIntSet i, UIntMatch j, UIntMatch k);
        
        const MatchBank<UIntSet, UIntMatch, ScoreFloat>* iteratee = nullptr;
        match_id_t match_id;
    };
    
    iterator begin() const;
    iterator end() const;
    
private:
    
    friend class iterator;
    
    const std::vector<match_set_t>* matches = nullptr;
    const std::unordered_set<std::tuple<size_t, size_t, size_t>>* masked_matches = nullptr;
    size_t num_match_sets = 0;
    
    std::vector<std::vector<std::vector<std::pair<ScoreFloat, match_id_t>>>> dp;
    std::vector<std::vector<std::pair<UIntSet, UIntMatch>>> starts;
    std::vector<std::vector<std::pair<UIntSet, UIntMatch>>> ends;
    
};

/*
 * Template implementations
 */


template<typename UIntSet, typename UIntMatch, typename ScoreFloat>
template<class BGraph>
MatchBank<UIntSet, UIntMatch, ScoreFloat>::MatchBank(const BGraph& graph1, const std::vector<match_set_t>& matches, size_t num_match_sets, bool suppress_verbose_logging, const std::unordered_set<std::tuple<size_t, size_t, size_t>>* masked_matches) noexcept : matches(&matches), starts(graph1.node_size()), ends(graph1.node_size()), dp(matches.size()), masked_matches(masked_matches), num_match_sets(num_match_sets) {
    
    static const bool debug = false;
    
    for (UIntSet i = 0; i < num_match_sets; ++i) {
        const auto& match_set = matches[i];
        dp[i].resize(match_set.walks1.size(),
                     std::vector<std::pair<ScoreFloat, match_id_t>>(match_set.walks2.size(),
                                                                    std::make_pair(std::numeric_limits<ScoreFloat>::lowest(), this->max())));
        for (size_t j = 0; j < match_set.walks1.size(); ++j) {
            starts[match_set.walks1[j].front()].emplace_back(i, j);
            ends[match_set.walks1[j].back()].emplace_back(i, j);
        }
    }
    
    for (auto& start_list : starts) {
        start_list.shrink_to_fit();
    }
    for (auto& end_list : ends) {
        end_list.shrink_to_fit();
    }
    
    if (!suppress_verbose_logging) {
        // dp table
        size_t dp_data_size = sizeof(dp) + dp.capacity() * sizeof(typename decltype(dp)::value_type);
        for (const auto& group : dp) {
            dp_data_size += group.capacity() * sizeof(typename decltype(dp)::value_type::value_type);
            for (const auto& row : group) {
                dp_data_size += row.capacity() * sizeof(typename decltype(dp)::value_type::value_type::value_type);
            }
        }
        
        // start and end locations
        size_t start_end_size = sizeof(starts) + sizeof(ends) + (starts.capacity() + ends.capacity()) * sizeof(typename decltype(starts)::value_type);
        assert(starts.size() == ends.size());
        for (size_t i = 0; i < starts.size(); ++i) {
            start_end_size += (starts[i].capacity() + ends[i].capacity()) * sizeof(typename decltype(starts)::value_type::value_type);
        }
        
        logging::log(logging::Debug, "Dynamic programming table is occupying " + format_memory_usage(dp_data_size) + ".");
        logging::log(logging::Debug, "Match start and end locations are occupying " + format_memory_usage(start_end_size) + ".");
    }
}


template<typename UIntSet, typename UIntMatch, typename ScoreFloat>
inline std::tuple<UIntSet, UIntMatch, UIntMatch> MatchBank<UIntSet, UIntMatch, ScoreFloat>::get_match_indexes(const match_id_t& match_id) const {
    return match_id;
}

template<typename UIntSet, typename UIntMatch, typename ScoreFloat>
inline const std::vector<uint64_t>& MatchBank<UIntSet, UIntMatch, ScoreFloat>::walk1(const match_id_t& match_id) const {
    return (*matches)[std::get<0>(match_id)].walks1[std::get<1>(match_id)];
}

template<typename UIntSet, typename UIntMatch, typename ScoreFloat>
inline const std::vector<uint64_t>& MatchBank<UIntSet, UIntMatch, ScoreFloat>::walk2(const match_id_t& match_id) const {
    return (*matches)[std::get<0>(match_id)].walks2[std::get<2>(match_id)];
}

template<typename UIntSet, typename UIntMatch, typename ScoreFloat>
inline const match_set_t& MatchBank<UIntSet, UIntMatch, ScoreFloat>::match_set(const match_id_t& match_id) const {
    return (*matches)[std::get<0>(match_id)];
}

template<typename UIntSet, typename UIntMatch, typename ScoreFloat>
inline ScoreFloat MatchBank<UIntSet, UIntMatch, ScoreFloat>::dp_value(const match_id_t& match_id) const {
    return dp[std::get<0>(match_id)][std::get<1>(match_id)][std::get<2>(match_id)].first;
}

template<typename UIntSet, typename UIntMatch, typename ScoreFloat>
inline const typename MatchBank<UIntSet, UIntMatch, ScoreFloat>::match_id_t& MatchBank<UIntSet, UIntMatch, ScoreFloat>::backpointer(const match_id_t& match_id) const {
    return dp[std::get<0>(match_id)][std::get<1>(match_id)][std::get<2>(match_id)].second;
}

template<typename UIntSet, typename UIntMatch, typename ScoreFloat>
inline void MatchBank<UIntSet, UIntMatch, ScoreFloat>::update_dp(const match_id_t& match_id, ScoreFloat value, const match_id_t& traceback_id) {
    static const bool debug = false;
    auto& dp_entry = dp[std::get<0>(match_id)][std::get<1>(match_id)][std::get<2>(match_id)];
    if (debug) {
        std::cerr << "check dp value at " << std::get<0>(match_id) << ", " << std::get<1>(match_id) << ", " << std::get<2>(match_id) << ", current " << dp_entry.first << " new value " << value << '\n';
    }
    if (value > dp_entry.first) {
        if (debug) {
            std::cerr << "dp value at " << std::get<0>(match_id) << ", " << std::get<1>(match_id) << ", " << std::get<2>(match_id) << " increases from " << dp_entry.first << " to " << value << '\n';
        }
        dp_entry.first = value;
        dp_entry.second = traceback_id;
    }
}


template<typename UIntSet, typename UIntMatch, typename ScoreFloat>
std::vector<typename MatchBank<UIntSet, UIntMatch, ScoreFloat>::match_id_t> MatchBank<UIntSet, UIntMatch, ScoreFloat>::starts_on(uint64_t node_id) const {
    std::vector<match_id_t> start_ids;
    for (const auto& start : starts[node_id]) {
        const auto& match_set = (*matches)[start.first];
        for (UIntMatch k = 0; k < match_set.walks2.size(); ++k) {
            if (masked_matches && masked_matches->count(std::tuple<size_t, size_t, size_t>(start.first, start.second, k))) {
                continue;
            }
            start_ids.emplace_back(start.first, start.second, k);
        }
    }
    return start_ids;
}

template<typename UIntSet, typename UIntMatch, typename ScoreFloat>
std::vector<typename MatchBank<UIntSet, UIntMatch, ScoreFloat>::match_id_t> MatchBank<UIntSet, UIntMatch, ScoreFloat>::ends_on(uint64_t node_id) const {
    std::vector<match_id_t> end_ids;
    for (const auto& end : ends[node_id]) {
        const auto& match_set = (*matches)[end.first];
        for (UIntMatch k = 0; k < match_set.walks2.size(); ++k) {
            if (masked_matches && masked_matches->count(std::tuple<size_t, size_t, size_t>(end.first, end.second, k))) {
                continue;
            }
            end_ids.emplace_back(end.first, end.second, k);
        }
    }
    return end_ids;
}




template<typename UIntSet, typename UIntMatch, typename ScoreFloat>
typename MatchBank<UIntSet, UIntMatch, ScoreFloat>::iterator MatchBank<UIntSet, UIntMatch, ScoreFloat>::begin() const {
    iterator it(*this, 0, 0, 0);
    if (masked_matches && it != end() && masked_matches->count(std::tuple<size_t, size_t, size_t>(0, 0, 0))) {
        ++it;
    }
    return it;
}

template<typename UIntSet, typename UIntMatch, typename ScoreFloat>
typename MatchBank<UIntSet, UIntMatch, ScoreFloat>::iterator MatchBank<UIntSet, UIntMatch, ScoreFloat>::end() const {
    return iterator(*this, num_match_sets, 0, 0);
}


template<typename UIntSet, typename UIntMatch, typename ScoreFloat>
inline typename MatchBank<UIntSet, UIntMatch, ScoreFloat>::match_id_t MatchBank<UIntSet, UIntMatch, ScoreFloat>::max() {
    return match_id_t(-1, -1, -1);
}

template<typename UIntSet, typename UIntMatch, typename ScoreFloat>
inline typename MatchBank<UIntSet, UIntMatch, ScoreFloat>::match_id_t MatchBank<UIntSet, UIntMatch, ScoreFloat>::min() {
    return match_id_t(0, 0, 0);
}

template<typename UIntSet, typename UIntMatch, typename ScoreFloat>
MatchBank<UIntSet, UIntMatch, ScoreFloat>::iterator::iterator(const MatchBank<UIntSet, UIntMatch, ScoreFloat>& iteratee, UIntSet i, UIntMatch j, UIntMatch k) : iteratee(&iteratee), match_id(i, j, k) {

}

template<typename UIntSet, typename UIntMatch, typename ScoreFloat>
typename MatchBank<UIntSet, UIntMatch, ScoreFloat>::iterator& MatchBank<UIntSet, UIntMatch, ScoreFloat>::iterator::operator++() {
    do {
        const auto& match_set = (*iteratee->matches)[std::get<0>(match_id)];
        std::get<2>(match_id)++;
        if (std::get<2>(match_id) == match_set.walks2.size()) {
            std::get<2>(match_id) = 0;
            ++std::get<1>(match_id);
            if (std::get<1>(match_id) == match_set.walks1.size()) {
                std::get<1>(match_id) = 0;
                ++std::get<0>(match_id);
            }
        }
    } while (*this != iteratee->end() && iteratee->masked_matches &&
             iteratee->masked_matches->count(std::tuple<size_t, size_t, size_t>(std::get<0>(match_id), std::get<1>(match_id), std::get<2>(match_id))));
    
    return *this;
}

template<typename UIntSet, typename UIntMatch, typename ScoreFloat>
const typename MatchBank<UIntSet, UIntMatch, ScoreFloat>::match_id_t& MatchBank<UIntSet, UIntMatch, ScoreFloat>::iterator::operator*() const {
    return match_id;
}

template<typename UIntSet, typename UIntMatch, typename ScoreFloat>
const typename MatchBank<UIntSet, UIntMatch, ScoreFloat>::match_id_t* MatchBank<UIntSet, UIntMatch, ScoreFloat>::iterator::operator->() const {
    return &match_id;
}

template<typename UIntSet, typename UIntMatch, typename ScoreFloat>
bool MatchBank<UIntSet, UIntMatch, ScoreFloat>::iterator::operator==(const iterator& other) const {
    return iteratee == other.iteratee && match_id == other.match_id;
}

template<typename UIntSet, typename UIntMatch, typename ScoreFloat>
bool MatchBank<UIntSet, UIntMatch, ScoreFloat>::iterator::operator!=(const iterator& other) const {
    return !(*this == other);
}


}

#endif /* centrolign_match_bank_hpp */
