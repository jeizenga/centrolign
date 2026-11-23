#ifndef centrolign_multi_orthogonal_max_search_tree_hpp
#define centrolign_multi_orthogonal_max_search_tree_hpp

#include <vector>
#include <array>
#include <cstdint>
#include <limits>
#include <algorithm>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include "centrolign/packed_vector.hpp"
#include "centrolign/utility.hpp"

namespace centrolign {

/*
 * The cross trees at every depth of a MultiOrthogonalMaxSearchTree, stored as interleaved records so that a
 * step of a cross tree search reads one record: the entry's key2, its index in the outer tree, and for each
 * channel the position of the max within the cross subtree rooted at the record. Each depth has a record for
 * every entry, and all of the depths share one buffer, which is large enough to use huge pages more often
 * than separate buffers would be.
 */
template<typename StoredK2, typename StoredIndex, size_t NumChannels>
class UnpackedCrossLevels {
public:
    UnpackedCrossLevels() = default;
    // the max positions at each depth are bounded by the corresponding entry of max_max_pos
    UnpackedCrossLevels(size_t size, uint64_t max_key2, uint64_t max_index, const std::vector<uint64_t>& max_max_pos);
    
    inline uint64_t key2(size_t depth, size_t pos) const;
    inline uint64_t index(size_t depth, size_t pos) const;
    inline uint64_t max_pos(size_t depth, size_t pos, size_t channel) const;
    
    inline void set_key2(size_t depth, size_t pos, uint64_t key2);
    inline void set_index(size_t depth, size_t pos, uint64_t index);
    inline void set_max_pos(size_t depth, size_t pos, size_t channel, uint64_t max_pos);
    
    // hint that the records in [pos_begin, pos_end) at this depth will be read soon
    inline void prefetch(size_t depth, size_t pos_begin, size_t pos_end) const;
    
    // the number of bytes for which huge pages were requested
    size_t huge_page_bytes() const;
    size_t memory_size() const;
    
private:
    
    struct Record {
        StoredK2 key2;
        StoredIndex index;
        std::array<StoredIndex, NumChannels> max_pos;
    };
    
    size_t size = 0;
    size_t advised_bytes = 0;
    std::vector<Record> records;
};

/*
 * Bit-packed version of UnpackedCrossLevels, in which each field of the records has the minimum width for its
 * maximum value (or the full width of the key type if PackKey2 is false), so it uses the same number of bits
 * as separate packed vectors would
 */
template<size_t NumChannels, bool PackKey2, size_t Key2Bits>
class PackedCrossLevels {
public:
    PackedCrossLevels() = default;
    // the max positions at each depth are bounded by the corresponding entry of max_max_pos
    PackedCrossLevels(size_t size, uint64_t max_key2, uint64_t max_index, const std::vector<uint64_t>& max_max_pos);
    
    inline uint64_t key2(size_t depth, size_t pos) const;
    inline uint64_t index(size_t depth, size_t pos) const;
    inline uint64_t max_pos(size_t depth, size_t pos, size_t channel) const;
    
    inline void set_key2(size_t depth, size_t pos, uint64_t key2);
    inline void set_index(size_t depth, size_t pos, uint64_t index);
    inline void set_max_pos(size_t depth, size_t pos, size_t channel, uint64_t max_pos);
    
    // hint that the records in [pos_begin, pos_end) at this depth will be read soon
    inline void prefetch(size_t depth, size_t pos_begin, size_t pos_end) const;
    
    // the number of bytes for which huge pages were requested
    size_t huge_page_bytes() const;
    size_t memory_size() const;
    
private:
    
    static inline uint8_t width_for(uint64_t max_value);
    
    // the bit where a record begins
    inline size_t record_bit(size_t depth, size_t pos) const;
    
    inline uint64_t get(size_t bit, uint8_t width) const;
    inline void set(size_t bit, uint8_t width, uint64_t value);
    
    uint8_t key2_width = 0;
    uint8_t index_width = 0;
    // the bit offset of the max positions within a record
    size_t max_pos_offset = 0;
    // for each depth, the width of the max positions, the width of a record, and the bit where the depth begins
    std::vector<uint8_t> max_pos_width;
    std::vector<size_t> record_width;
    std::vector<size_t> depth_begin;
    
    size_t advised_bytes = 0;
    std::vector<uint64_t> words;
};

/*
 * Chooses the bit-packed cross levels when the tree's index vectors are bit-packed
 */
template<typename K2, typename UIntSize, size_t NumChannels, class KeyVector2, class IndexVector>
struct CrossLevelsSelector {
    using type = UnpackedCrossLevels<typename KeyVector2::value_type, typename IndexVector::value_type, NumChannels>;
};

template<typename K2, typename UIntSize, size_t NumChannels, class KeyVector2>
struct CrossLevelsSelector<K2, UIntSize, NumChannels, KeyVector2, PackedVector> {
    using type = PackedCrossLevels<NumChannels, std::is_same<KeyVector2, PackedVector>::value, 8 * sizeof(K2)>;
};

/*
 * Static-topology search tree that supports 2D orthogonal range max queries over
 * [lo1, hi1) x [0, hi2) and value updates. Each entry has several channels of values
 * that share the same keys. The channels are organized into groups, and each query computes
 * the max for every channel in one group. Values can only be increased.
 *
 * Internally, this is a range tree in which the cross trees of all outer nodes at the
 * same depth are stored in shared arrays (each subtree's cross tree occupies the range of
 * sorted positions covered by the subtree), and the cross trees refer back to the values
 * in the outer tree rather than copying them.
 */
template<typename K1, typename K2, typename V, typename UIntSize, size_t NumGroups, size_t GroupSize,
         class KeyVector1 = std::vector<K1>, class KeyVector2 = std::vector<K2>, class IndexVector = std::vector<UIntSize>>
class MultiOrthogonalMaxSearchTree {
public:

    static const size_t NumChannels = NumGroups * GroupSize;

    // vector will be sorted, key pairs must be unique, every channel is initialized to the initial value
    MultiOrthogonalMaxSearchTree(std::vector<std::pair<K1, K2>>& data, const V& initial_value);
    MultiOrthogonalMaxSearchTree() = default;
    MultiOrthogonalMaxSearchTree(const MultiOrthogonalMaxSearchTree& other) = default;
    MultiOrthogonalMaxSearchTree(MultiOrthogonalMaxSearchTree&& other) = default;
    ~MultiOrthogonalMaxSearchTree() = default;
    MultiOrthogonalMaxSearchTree& operator=(const MultiOrthogonalMaxSearchTree& other) = default;
    MultiOrthogonalMaxSearchTree& operator=(MultiOrthogonalMaxSearchTree&& other) = default;

    // returns the index of the entry with these keys, or size() if there is none
    UIntSize find(const K1& key1, const K2& key2) const;
    // increase the values of an entry, where the value for channel c of group g is at index g * GroupSize + c,
    // any value that does not exceed the entry's current value for the channel is ignored
    void update(UIntSize i, const std::array<V, NumChannels>& values);
    // for each channel in the group, returns the (value, entry index) of the max value in [lo1, hi1) x [0, hi2),
    // or (lowest, size()) if there are no entries in the range
    std::array<std::pair<V, UIntSize>, GroupSize> range_max(const K1& lo1, const K1& hi1, const K2& hi2, size_t group) const;

    // the keys and values of an entry
    K1 key1_at(UIntSize i) const;
    K2 key2_at(UIntSize i) const;
    V value_at(UIntSize i, size_t group, size_t channel) const;

    size_t size() const;
    bool empty() const;

    size_t memory_size() const;
    // the number of bytes for which huge pages were requested
    size_t huge_page_bytes() const;

private:

    // subtrees smaller than this are scanned directly rather than having cross trees
    static const size_t min_cross_size = 8;
    // searches in the cross trees prefetch the records this many levels below the current node
    static const size_t prefetch_depth = 3;

    // an outer tree node, along with the range of sorted positions covered by its subtree
    struct Subtree {
        size_t node;
        size_t lo;
        size_t size;
        size_t depth;
    };

    inline Subtree root() const;
    static inline Subtree left_subtree(const Subtree& subtree);
    static inline Subtree right_subtree(const Subtree& subtree);

    // the size of the left subtree of a complete binary tree of this size
    static inline size_t left_size(size_t size);
    static inline size_t depth(size_t x);
    static inline size_t left(size_t x);
    static inline size_t right(size_t x);
    static inline size_t parent(size_t x);

    // calls func(rank, x) for the nodes x of a heap-ordered binary tree of the given size in in-order
    template<class Func>
    static void in_order(size_t size, const Func& func);

    inline V value_of(size_t i, size_t group, size_t channel) const;
    // hint that the values of an entry's group will be read soon
    inline void prefetch_value(size_t i, size_t group) const;
    // hint that the records of a cross tree that are a few levels below a node will be read soon (they are
    // contiguous in heap order)
    inline void prefetch_descendants(const Subtree& subtree, size_t x) const;
    // replace the current maxes of the group's channels with this entry's values if they are higher
    inline void consider(size_t i, size_t group, std::array<std::pair<V, UIntSize>, GroupSize>& maxes) const;
    inline void consider(size_t i, size_t group, size_t channel, std::array<std::pair<V, UIntSize>, GroupSize>& maxes) const;
    // consider all of the entries in this subtree with key2 in [0, hi2)
    inline void subtree_max(const Subtree& subtree, const K2& hi2, size_t group,
                            std::array<std::pair<V, UIntSize>, GroupSize>& maxes) const;

    // the keys of the outer tree
    KeyVector1 key1;
    KeyVector2 key2;
    // the values of all channels for each entry, stored contiguously
    std::vector<V> value;
    // for each depth, the cross trees for all subtrees at that depth, each occupying its range of sorted positions
    using CrossLevels = typename CrossLevelsSelector<K2, UIntSize, NumChannels, KeyVector2, IndexVector>::type;
    CrossLevels cross_levels;

    static_assert(std::is_integral<UIntSize>::value, "UIntSize must be an integer");
    static_assert(std::is_unsigned<UIntSize>::value, "UIntSize must be unsigned");
    static_assert(std::is_integral<K2>::value, "K2 must be an integer");
    static_assert(std::is_unsigned<K2>::value, "K2 must be unsigned");
};




/*
 * Template implementations
 */

template<typename StoredK2, typename StoredIndex, size_t NumChannels>
UnpackedCrossLevels<StoredK2, StoredIndex, NumChannels>::UnpackedCrossLevels(size_t size, uint64_t max_key2, uint64_t max_index,
                                                                             const std::vector<uint64_t>& max_max_pos) : size(size) {
    size_t num_records = size * max_max_pos.size();
    // request huge pages before the memory is first touched
    records.reserve(num_records);
    advised_bytes = advise_huge_pages(records.data(), num_records * sizeof(Record));
    records.resize(num_records);
}

template<typename StoredK2, typename StoredIndex, size_t NumChannels>
inline uint64_t UnpackedCrossLevels<StoredK2, StoredIndex, NumChannels>::key2(size_t depth, size_t pos) const {
    return records[depth * size + pos].key2;
}

template<typename StoredK2, typename StoredIndex, size_t NumChannels>
inline uint64_t UnpackedCrossLevels<StoredK2, StoredIndex, NumChannels>::index(size_t depth, size_t pos) const {
    return records[depth * size + pos].index;
}

template<typename StoredK2, typename StoredIndex, size_t NumChannels>
inline uint64_t UnpackedCrossLevels<StoredK2, StoredIndex, NumChannels>::max_pos(size_t depth, size_t pos, size_t channel) const {
    return records[depth * size + pos].max_pos[channel];
}

template<typename StoredK2, typename StoredIndex, size_t NumChannels>
inline void UnpackedCrossLevels<StoredK2, StoredIndex, NumChannels>::set_key2(size_t depth, size_t pos, uint64_t key2) {
    records[depth * size + pos].key2 = key2;
}

template<typename StoredK2, typename StoredIndex, size_t NumChannels>
inline void UnpackedCrossLevels<StoredK2, StoredIndex, NumChannels>::set_index(size_t depth, size_t pos, uint64_t index) {
    records[depth * size + pos].index = index;
}

template<typename StoredK2, typename StoredIndex, size_t NumChannels>
inline void UnpackedCrossLevels<StoredK2, StoredIndex, NumChannels>::set_max_pos(size_t depth, size_t pos, size_t channel, uint64_t max_pos) {
    records[depth * size + pos].max_pos[channel] = max_pos;
}

template<typename StoredK2, typename StoredIndex, size_t NumChannels>
inline void UnpackedCrossLevels<StoredK2, StoredIndex, NumChannels>::prefetch(size_t depth, size_t pos_begin, size_t pos_end) const {
    // intentionally a no-op: the unpacked records are large enough that prefetching several of them per
    // search step cost more in memory traffic than it saved in latency when benchmarked
}

template<typename StoredK2, typename StoredIndex, size_t NumChannels>
size_t UnpackedCrossLevels<StoredK2, StoredIndex, NumChannels>::huge_page_bytes() const {
    return advised_bytes;
}

template<typename StoredK2, typename StoredIndex, size_t NumChannels>
size_t UnpackedCrossLevels<StoredK2, StoredIndex, NumChannels>::memory_size() const {
    return sizeof(*this) - sizeof(records) + get_vec_memory_size(records);
}

template<size_t NumChannels, bool PackKey2, size_t Key2Bits>
PackedCrossLevels<NumChannels, PackKey2, Key2Bits>::PackedCrossLevels(size_t size, uint64_t max_key2, uint64_t max_index,
                                                                      const std::vector<uint64_t>& max_max_pos) :
    key2_width(PackKey2 ? width_for(max_key2) : Key2Bits), index_width(width_for(max_index)),
    max_pos_width(max_max_pos.size()), record_width(max_max_pos.size()), depth_begin(max_max_pos.size()) {
    
    max_pos_offset = key2_width + index_width;
    size_t num_bits = 0;
    for (size_t d = 0; d < max_max_pos.size(); ++d) {
        max_pos_width[d] = width_for(max_max_pos[d]);
        record_width[d] = max_pos_offset + NumChannels * max_pos_width[d];
        depth_begin[d] = num_bits;
        num_bits += size * record_width[d];
    }
    size_t num_words = (num_bits + 63) / 64;
    // request huge pages before the memory is first touched
    words.reserve(num_words);
    advised_bytes = advise_huge_pages(words.data(), num_words * sizeof(uint64_t));
    words.resize(num_words, 0);
}

template<size_t NumChannels, bool PackKey2, size_t Key2Bits>
inline uint8_t PackedCrossLevels<NumChannels, PackKey2, Key2Bits>::width_for(uint64_t max_value) {
    // matches the widths chosen by SizedVectorFactory<PackedVector>
    uint8_t width = 1;
    while (width < 64 && (max_value >> width) != 0) {
        ++width;
    }
    return width;
}

template<size_t NumChannels, bool PackKey2, size_t Key2Bits>
inline size_t PackedCrossLevels<NumChannels, PackKey2, Key2Bits>::record_bit(size_t depth, size_t pos) const {
    return depth_begin[depth] + pos * record_width[depth];
}

template<size_t NumChannels, bool PackKey2, size_t Key2Bits>
inline uint64_t PackedCrossLevels<NumChannels, PackKey2, Key2Bits>::get(size_t bit, uint8_t width) const {
    size_t k = bit / 64;
    size_t s = bit % 64;
    uint64_t value = words[k] >> s;
    if (s + width > 64) {
        // the field continues into the next integer
        value |= words[k + 1] << (64 - s);
    }
    return width == 64 ? value : value & ((uint64_t(1) << width) - 1);
}

template<size_t NumChannels, bool PackKey2, size_t Key2Bits>
inline void PackedCrossLevels<NumChannels, PackKey2, Key2Bits>::set(size_t bit, uint8_t width, uint64_t value) {
    size_t k = bit / 64;
    size_t s = bit % 64;
    uint64_t mask = width == 64 ? uint64_t(-1) : (uint64_t(1) << width) - 1;
    value &= mask;
    words[k] = (words[k] & ~(mask << s)) | (value << s);
    if (s + width > 64) {
        // the field continues into the next integer
        words[k + 1] = (words[k + 1] & ~(mask >> (64 - s))) | (value >> (64 - s));
    }
}

template<size_t NumChannels, bool PackKey2, size_t Key2Bits>
inline uint64_t PackedCrossLevels<NumChannels, PackKey2, Key2Bits>::key2(size_t depth, size_t pos) const {
    return get(record_bit(depth, pos), key2_width);
}

template<size_t NumChannels, bool PackKey2, size_t Key2Bits>
inline uint64_t PackedCrossLevels<NumChannels, PackKey2, Key2Bits>::index(size_t depth, size_t pos) const {
    return get(record_bit(depth, pos) + key2_width, index_width);
}

template<size_t NumChannels, bool PackKey2, size_t Key2Bits>
inline uint64_t PackedCrossLevels<NumChannels, PackKey2, Key2Bits>::max_pos(size_t depth, size_t pos, size_t channel) const {
    return get(record_bit(depth, pos) + max_pos_offset + channel * max_pos_width[depth], max_pos_width[depth]);
}

template<size_t NumChannels, bool PackKey2, size_t Key2Bits>
inline void PackedCrossLevels<NumChannels, PackKey2, Key2Bits>::set_key2(size_t depth, size_t pos, uint64_t key2) {
    set(record_bit(depth, pos), key2_width, key2);
}

template<size_t NumChannels, bool PackKey2, size_t Key2Bits>
inline void PackedCrossLevels<NumChannels, PackKey2, Key2Bits>::set_index(size_t depth, size_t pos, uint64_t index) {
    set(record_bit(depth, pos) + key2_width, index_width, index);
}

template<size_t NumChannels, bool PackKey2, size_t Key2Bits>
inline void PackedCrossLevels<NumChannels, PackKey2, Key2Bits>::set_max_pos(size_t depth, size_t pos, size_t channel, uint64_t max_pos) {
    set(record_bit(depth, pos) + max_pos_offset + channel * max_pos_width[depth], max_pos_width[depth], max_pos);
}

template<size_t NumChannels, bool PackKey2, size_t Key2Bits>
inline void PackedCrossLevels<NumChannels, PackKey2, Key2Bits>::prefetch(size_t depth, size_t pos_begin, size_t pos_end) const {
    prefetch_range(words.data() + record_bit(depth, pos_begin) / 64, words.data() + (record_bit(depth, pos_end) + 63) / 64);
}

template<size_t NumChannels, bool PackKey2, size_t Key2Bits>
size_t PackedCrossLevels<NumChannels, PackKey2, Key2Bits>::huge_page_bytes() const {
    return advised_bytes;
}

template<size_t NumChannels, bool PackKey2, size_t Key2Bits>
size_t PackedCrossLevels<NumChannels, PackKey2, Key2Bits>::memory_size() const {
    return (sizeof(*this) + words.capacity() * sizeof(uint64_t) + max_pos_width.capacity() * sizeof(uint8_t)
            + record_width.capacity() * sizeof(size_t) + depth_begin.capacity() * sizeof(size_t));
}

template<typename K1, typename K2, typename V, typename UIntSize, size_t NumGroups, size_t GroupSize, class KeyVector1, class KeyVector2, class IndexVector>
MultiOrthogonalMaxSearchTree<K1, K2, V, UIntSize, NumGroups, GroupSize, KeyVector1, KeyVector2, IndexVector>::MultiOrthogonalMaxSearchTree(std::vector<std::pair<K1, K2>>& data,
                                                                                                                                           const V& initial_value) : key1(data.size()), value(data.size() * NumChannels, initial_value) {

    if (size_t(std::numeric_limits<UIntSize>::max()) < data.size()) {
        throw std::runtime_error("MultiOrthogonalMaxSearchTree integer size is too small to represent its indexes.");
    }

    if (data.empty()) {
        return;
    }

    if (!std::is_sorted(data.begin(), data.end())) {
        std::sort(data.begin(), data.end());
    }

    size_t n = data.size();
    uint64_t max_key2 = 0;
    for (const auto& item : data) {
        max_key2 = std::max<uint64_t>(max_key2, item.second);
    }

    // assign the keys to the outer tree in sorted order
    key2 = std::move(SizedVectorFactory<KeyVector2>::make(n, max_key2));
    in_order(n, [&](size_t r, size_t x) {
        key1[x] = data[r].first;
        key2[x] = data[r].second;
    });
    const auto& const_key2 = key2;

    // record the range of sorted positions that each subtree covers
    std::vector<UIntSize> subtree_lo(n), subtree_size(n);
    subtree_lo[0] = 0;
    subtree_size[0] = n;
    for (size_t x = 0; x < n; ++x) {
        Subtree subtree{x, subtree_lo[x], subtree_size[x], 0};
        auto l = left_subtree(subtree);
        if (l.node < n) {
            subtree_lo[l.node] = l.lo;
            subtree_size[l.node] = l.size;
        }
        auto r = right_subtree(subtree);
        if (r.node < n) {
            subtree_lo[r.node] = r.lo;
            subtree_size[r.node] = r.size;
        }
    }

    // the outer indexes, sorted by key2 with ties broken by index
    std::vector<UIntSize> order(n);
    for (size_t i = 0; i < n; ++i) {
        order[i] = i;
    }
    std::sort(order.begin(), order.end(), [&](UIntSize a, UIntSize b) {
        K2 key_a = const_key2[a];
        K2 key_b = const_key2[b];
        return key_a < key_b || (key_a == key_b && a < b);
    });

    // the leftmost subtree at each depth is the largest, so it determines which depths have any cross trees
    size_t num_levels = 0;
    while ((size_t(1) << num_levels) - 1 < n && subtree_size[(size_t(1) << num_levels) - 1] >= min_cross_size) {
        ++num_levels;
    }

    std::vector<uint64_t> max_max_pos(num_levels);
    for (size_t d = 0; d < num_levels; ++d) {
        max_max_pos[d] = subtree_size[(size_t(1) << d) - 1] - 1;
    }
    cross_levels = std::move(CrossLevels(n, max_key2, n - 1, max_max_pos));

    std::vector<UIntSize> level_order(n);
    std::vector<UIntSize> num_filled;
    for (size_t d = 0; d < num_levels; ++d) {

        size_t level_begin = (size_t(1) << d) - 1;
        size_t level_end = std::min(n, 2 * level_begin + 1);

        // divvy up the sorted outer indexes into the cross trees for the subtrees at this depth
        num_filled.assign(level_end - level_begin, 0);
        for (UIntSize i : order) {
            size_t depth_i = depth(i);
            if (depth_i < d) {
                continue;
            }
            size_t ancestor = ((size_t(i) + 1) >> (depth_i - d)) - 1;
            if (subtree_size[ancestor] < min_cross_size) {
                continue;
            }
            level_order[subtree_lo[ancestor] + num_filled[ancestor - level_begin]++] = i;
        }


        // lay out each cross tree in heap order
        for (size_t v = level_begin; v < level_end; ++v) {
            size_t lo = subtree_lo[v];
            if (subtree_size[v] < min_cross_size) {
                continue;
            }
            in_order(subtree_size[v], [&](size_t r, size_t x) {
                UIntSize i = level_order[lo + r];
                cross_levels.set_key2(d, lo + x, K2(const_key2[i]));
                cross_levels.set_index(d, lo + x, i);
                // all values are equal, so each node can initially be its own max
                for (size_t ch = 0; ch < NumChannels; ++ch) {
                    cross_levels.set_max_pos(d, lo + x, ch, x);
                }
            });
        }
    }
}

template<typename K1, typename K2, typename V, typename UIntSize, size_t NumGroups, size_t GroupSize, class KeyVector1, class KeyVector2, class IndexVector>
inline typename MultiOrthogonalMaxSearchTree<K1, K2, V, UIntSize, NumGroups, GroupSize, KeyVector1, KeyVector2, IndexVector>::Subtree
MultiOrthogonalMaxSearchTree<K1, K2, V, UIntSize, NumGroups, GroupSize, KeyVector1, KeyVector2, IndexVector>::root() const {
    return Subtree{0, 0, size(), 0};
}

template<typename K1, typename K2, typename V, typename UIntSize, size_t NumGroups, size_t GroupSize, class KeyVector1, class KeyVector2, class IndexVector>
inline typename MultiOrthogonalMaxSearchTree<K1, K2, V, UIntSize, NumGroups, GroupSize, KeyVector1, KeyVector2, IndexVector>::Subtree
MultiOrthogonalMaxSearchTree<K1, K2, V, UIntSize, NumGroups, GroupSize, KeyVector1, KeyVector2, IndexVector>::left_subtree(const Subtree& subtree) {
    return Subtree{left(subtree.node), subtree.lo, left_size(subtree.size), subtree.depth + 1};
}

template<typename K1, typename K2, typename V, typename UIntSize, size_t NumGroups, size_t GroupSize, class KeyVector1, class KeyVector2, class IndexVector>
inline typename MultiOrthogonalMaxSearchTree<K1, K2, V, UIntSize, NumGroups, GroupSize, KeyVector1, KeyVector2, IndexVector>::Subtree
MultiOrthogonalMaxSearchTree<K1, K2, V, UIntSize, NumGroups, GroupSize, KeyVector1, KeyVector2, IndexVector>::right_subtree(const Subtree& subtree) {
    if (subtree.size == 0) {
        return Subtree{right(subtree.node), subtree.lo, 0, subtree.depth + 1};
    }
    size_t l = left_size(subtree.size);
    return Subtree{right(subtree.node), subtree.lo + l + 1, subtree.size - l - 1, subtree.depth + 1};
}

template<typename K1, typename K2, typename V, typename UIntSize, size_t NumGroups, size_t GroupSize, class KeyVector1, class KeyVector2, class IndexVector>
inline size_t MultiOrthogonalMaxSearchTree<K1, K2, V, UIntSize, NumGroups, GroupSize, KeyVector1, KeyVector2, IndexVector>::left_size(size_t size) {
    if (size <= 1) {
        return 0;
    }
    // the last row of a complete tree is filled from the left
    size_t height = 63 - __builtin_clzll(size);
    size_t half_row = size_t(1) << (height - 1);
    size_t last_row = size - ((size_t(1) << height) - 1);
    return half_row - 1 + std::min(last_row, half_row);
}

template<typename K1, typename K2, typename V, typename UIntSize, size_t NumGroups, size_t GroupSize, class KeyVector1, class KeyVector2, class IndexVector>
inline size_t MultiOrthogonalMaxSearchTree<K1, K2, V, UIntSize, NumGroups, GroupSize, KeyVector1, KeyVector2, IndexVector>::depth(size_t x) {
    return 63 - __builtin_clzll(x + 1);
}

template<typename K1, typename K2, typename V, typename UIntSize, size_t NumGroups, size_t GroupSize, class KeyVector1, class KeyVector2, class IndexVector>
inline size_t MultiOrthogonalMaxSearchTree<K1, K2, V, UIntSize, NumGroups, GroupSize, KeyVector1, KeyVector2, IndexVector>::left(size_t x) {
    return 2 * x + 1;
}

template<typename K1, typename K2, typename V, typename UIntSize, size_t NumGroups, size_t GroupSize, class KeyVector1, class KeyVector2, class IndexVector>
inline size_t MultiOrthogonalMaxSearchTree<K1, K2, V, UIntSize, NumGroups, GroupSize, KeyVector1, KeyVector2, IndexVector>::right(size_t x) {
    return 2 * x + 2;
}

template<typename K1, typename K2, typename V, typename UIntSize, size_t NumGroups, size_t GroupSize, class KeyVector1, class KeyVector2, class IndexVector>
inline size_t MultiOrthogonalMaxSearchTree<K1, K2, V, UIntSize, NumGroups, GroupSize, KeyVector1, KeyVector2, IndexVector>::parent(size_t x) {
    return (x - 1) / 2;
}

template<typename K1, typename K2, typename V, typename UIntSize, size_t NumGroups, size_t GroupSize, class KeyVector1, class KeyVector2, class IndexVector>
template<class Func>
void MultiOrthogonalMaxSearchTree<K1, K2, V, UIntSize, NumGroups, GroupSize, KeyVector1, KeyVector2, IndexVector>::in_order(size_t size, const Func& func) {
    if (size == 0) {
        return;
    }
    // start at the leftmost node
    size_t x = 0;
    while (left(x) < size) {
        x = left(x);
    }
    for (size_t r = 0; r < size; ++r) {
        func(r, x);
        if (right(x) < size) {
            // the successor is the leftmost node in the right subtree
            x = right(x);
            while (left(x) < size) {
                x = left(x);
            }
        }
        else {
            // the successor is the first ancestor that we reach from a left child
            while (x != 0 && x % 2 == 0) {
                x = parent(x);
            }
            x = parent(x);
        }
    }
}

template<typename K1, typename K2, typename V, typename UIntSize, size_t NumGroups, size_t GroupSize, class KeyVector1, class KeyVector2, class IndexVector>
inline V MultiOrthogonalMaxSearchTree<K1, K2, V, UIntSize, NumGroups, GroupSize, KeyVector1, KeyVector2, IndexVector>::value_of(size_t i, size_t group, size_t channel) const {
    return value[i * NumChannels + group * GroupSize + channel];
}

template<typename K1, typename K2, typename V, typename UIntSize, size_t NumGroups, size_t GroupSize, class KeyVector1, class KeyVector2, class IndexVector>
inline void MultiOrthogonalMaxSearchTree<K1, K2, V, UIntSize, NumGroups, GroupSize, KeyVector1, KeyVector2, IndexVector>::prefetch_value(size_t i, size_t group) const {
    const V* group_values = value.data() + i * NumChannels + group * GroupSize;
    prefetch_range(group_values, group_values + GroupSize);
}

template<typename K1, typename K2, typename V, typename UIntSize, size_t NumGroups, size_t GroupSize, class KeyVector1, class KeyVector2, class IndexVector>
inline void MultiOrthogonalMaxSearchTree<K1, K2, V, UIntSize, NumGroups, GroupSize, KeyVector1, KeyVector2, IndexVector>::prefetch_descendants(const Subtree& subtree, size_t x) const {
    size_t first = ((x + 1) << prefetch_depth) - 1;
    if (first < subtree.size) {
        size_t last = std::min<size_t>(subtree.size, first + (size_t(1) << prefetch_depth));
        cross_levels.prefetch(subtree.depth, subtree.lo + first, subtree.lo + last);
    }
}

template<typename K1, typename K2, typename V, typename UIntSize, size_t NumGroups, size_t GroupSize, class KeyVector1, class KeyVector2, class IndexVector>
size_t MultiOrthogonalMaxSearchTree<K1, K2, V, UIntSize, NumGroups, GroupSize, KeyVector1, KeyVector2, IndexVector>::huge_page_bytes() const {
    return cross_levels.huge_page_bytes();
}

template<typename K1, typename K2, typename V, typename UIntSize, size_t NumGroups, size_t GroupSize, class KeyVector1, class KeyVector2, class IndexVector>
inline void MultiOrthogonalMaxSearchTree<K1, K2, V, UIntSize, NumGroups, GroupSize, KeyVector1, KeyVector2, IndexVector>::consider(size_t i, size_t group, size_t channel,
                                                                                                                                   std::array<std::pair<V, UIntSize>, GroupSize>& maxes) const {
    V val = value_of(i, group, channel);
    auto& current = maxes[channel];
    if (current.second == size() || val > current.first) {
        current.first = val;
        current.second = i;
    }
}

template<typename K1, typename K2, typename V, typename UIntSize, size_t NumGroups, size_t GroupSize, class KeyVector1, class KeyVector2, class IndexVector>
inline void MultiOrthogonalMaxSearchTree<K1, K2, V, UIntSize, NumGroups, GroupSize, KeyVector1, KeyVector2, IndexVector>::consider(size_t i, size_t group,
                                                                                                                                   std::array<std::pair<V, UIntSize>, GroupSize>& maxes) const {
    for (size_t c = 0; c < GroupSize; ++c) {
        consider(i, group, c, maxes);
    }
}

template<typename K1, typename K2, typename V, typename UIntSize, size_t NumGroups, size_t GroupSize, class KeyVector1, class KeyVector2, class IndexVector>
inline void MultiOrthogonalMaxSearchTree<K1, K2, V, UIntSize, NumGroups, GroupSize, KeyVector1, KeyVector2, IndexVector>::subtree_max(const Subtree& subtree, const K2& hi2, size_t group,
                                                                                                                                      std::array<std::pair<V, UIntSize>, GroupSize>& maxes) const {

    if (subtree.size >= min_cross_size) {
        // descend the cross tree, which is sorted by key2
        size_t d = subtree.depth;
        size_t x = 0;
        while (x < subtree.size) {
            prefetch_descendants(subtree, x);
            size_t pos = subtree.lo + x;
            if (K2(cross_levels.key2(d, pos)) < hi2) {
                // this node and its entire left subtree are in the range, so find their entries and request
                // their values before considering them in the same order as otherwise
                size_t i = cross_levels.index(d, pos);
                prefetch_value(i, group);
                std::array<size_t, GroupSize> max_entries;
                size_t l = left(x);
                if (l < subtree.size) {
                    for (size_t c = 0; c < GroupSize; ++c) {
                        size_t max_pos = subtree.lo + cross_levels.max_pos(d, subtree.lo + l, group * GroupSize + c);
                        max_entries[c] = cross_levels.index(d, max_pos);
                        prefetch_value(max_entries[c], group);
                    }
                }
                consider(i, group, maxes);
                if (l < subtree.size) {
                    for (size_t c = 0; c < GroupSize; ++c) {
                        consider(max_entries[c], group, c, maxes);
                    }
                }
                x = right(x);
            }
            else {
                x = left(x);
            }
        }
    }
    else {
        // the subtree is small, so we just scan it, one row at a time
        for (size_t row_begin = subtree.node, row_size = 1; row_begin < size(); row_begin = left(row_begin), row_size *= 2) {
            for (size_t x = row_begin, row_end = std::min(size(), row_begin + row_size); x < row_end; ++x) {
                if (K2(key2[x]) < hi2) {
                    consider(x, group, maxes);
                }
            }
        }
    }
}

template<typename K1, typename K2, typename V, typename UIntSize, size_t NumGroups, size_t GroupSize, class KeyVector1, class KeyVector2, class IndexVector>
UIntSize MultiOrthogonalMaxSearchTree<K1, K2, V, UIntSize, NumGroups, GroupSize, KeyVector1, KeyVector2, IndexVector>::find(const K1& search_key1, const K2& search_key2) const {
    size_t x = 0;
    while (x < size()) {
        K1 k1 = key1[x];
        K2 k2 = key2[x];
        if (k1 == search_key1 && k2 == search_key2) {
            return x;
        }
        else if (search_key1 < k1 || (search_key1 == k1 && search_key2 < k2)) {
            x = left(x);
        }
        else {
            x = right(x);
        }
    }
    return size();
}

template<typename K1, typename K2, typename V, typename UIntSize, size_t NumGroups, size_t GroupSize, class KeyVector1, class KeyVector2, class IndexVector>
void MultiOrthogonalMaxSearchTree<K1, K2, V, UIntSize, NumGroups, GroupSize, KeyVector1, KeyVector2, IndexVector>::update(UIntSize i, const std::array<V, NumChannels>& values) {

    std::array<bool, NumChannels> increased;
    bool any_increased = false;
    for (size_t ch = 0; ch < NumChannels; ++ch) {
        auto& val = value[size_t(i) * NumChannels + ch];
        increased[ch] = (values[ch] > val);
        if (increased[ch]) {
            val = values[ch];
            any_increased = true;
        }
    }
    if (!any_increased) {
        return;
    }

    const auto& const_key2 = key2;
    K2 k2 = const_key2[i];

    // update the cross trees of all of the subtrees that contain this entry, from the root downward
    size_t depth_i = depth(i);
    Subtree subtree = root();
    for (size_t d = 0; d <= depth_i; ++d) {
        if (d != 0) {
            size_t ancestor = ((size_t(i) + 1) >> (depth_i - d)) - 1;
            subtree = (ancestor == left(subtree.node)) ? left_subtree(subtree) : right_subtree(subtree);
        }
        if (subtree.size < min_cross_size) {
            // all of the deeper subtrees are smaller, so they don't have cross trees either
            break;
        }

        // find the entry's position in the cross tree, which is sorted by key2 and then index
        size_t p = 0;
        while (true) {
            prefetch_descendants(subtree, p);
            size_t pos = subtree.lo + p;
            UIntSize j = cross_levels.index(d, pos);
            if (j == i) {
                break;
            }
            K2 kj = cross_levels.key2(d, pos);
            p = (k2 < kj || (k2 == kj && i < j)) ? left(p) : right(p);
        }

        for (size_t g = 0; g < NumGroups; ++g) {
            for (size_t c = 0; c < GroupSize; ++c) {
                if (!increased[g * GroupSize + c]) {
                    continue;
                }
                V val = value_of(i, g, c);
                // walk upward until reaching a subtree whose max is already at least as high
                size_t x = p;
                while (true) {
                    size_t max_pos = cross_levels.max_pos(d, subtree.lo + x, g * GroupSize + c);
                    if (max_pos != p) {
                        if (val > value_of(cross_levels.index(d, subtree.lo + max_pos), g, c)) {
                            cross_levels.set_max_pos(d, subtree.lo + x, g * GroupSize + c, p);
                        }
                        else {
                            // values only increase, so the higher subtrees' maxes are still correct
                            break;
                        }
                    }
                    if (x == 0) {
                        break;
                    }
                    x = parent(x);
                }
            }
        }
    }
}

template<typename K1, typename K2, typename V, typename UIntSize, size_t NumGroups, size_t GroupSize, class KeyVector1, class KeyVector2, class IndexVector>
std::array<std::pair<V, UIntSize>, GroupSize>
MultiOrthogonalMaxSearchTree<K1, K2, V, UIntSize, NumGroups, GroupSize, KeyVector1, KeyVector2, IndexVector>::range_max(const K1& lo1, const K1& hi1, const K2& hi2, size_t group) const {

    std::array<std::pair<V, UIntSize>, GroupSize> maxes;
    for (auto& max_val : maxes) {
        max_val.first = std::numeric_limits<V>::lowest();
        max_val.second = size();
    }

    // traverse downward until finding the first value in the key1 interval
    Subtree cursor = root();
    while (cursor.node < size()) {
        K1 k1 = key1[cursor.node];
        if (k1 >= hi1) {
            cursor = left_subtree(cursor);
        }
        else if (k1 < lo1) {
            cursor = right_subtree(cursor);
        }
        else {
            break;
        }
    }

    if (cursor.node >= size()) {
        // we never found a value that was in the interval
        return maxes;
    }

    if (K2(key2[cursor.node]) < hi2) {
        consider(cursor.node, group, maxes);
    }

    // split off the leftward traversal where rightward off-path subtrees are included entirely
    Subtree left_cursor = left_subtree(cursor);
    while (left_cursor.node < size()) {
        if (K1(key1[left_cursor.node]) >= lo1) {
            if (K2(key2[left_cursor.node]) < hi2) {
                consider(left_cursor.node, group, maxes);
            }
            Subtree r = right_subtree(left_cursor);
            if (r.node < size()) {
                subtree_max(r, hi2, group, maxes);
            }
            left_cursor = left_subtree(left_cursor);
        }
        else {
            left_cursor = right_subtree(left_cursor);
        }
    }

    // and next the rightward traversal where leftward off-path subtrees are included entirely
    Subtree right_cursor = right_subtree(cursor);
    while (right_cursor.node < size()) {
        if (K1(key1[right_cursor.node]) < hi1) {
            if (K2(key2[right_cursor.node]) < hi2) {
                consider(right_cursor.node, group, maxes);
            }
            Subtree l = left_subtree(right_cursor);
            if (l.node < size()) {
                subtree_max(l, hi2, group, maxes);
            }
            right_cursor = right_subtree(right_cursor);
        }
        else {
            right_cursor = left_subtree(right_cursor);
        }
    }

    return maxes;
}

template<typename K1, typename K2, typename V, typename UIntSize, size_t NumGroups, size_t GroupSize, class KeyVector1, class KeyVector2, class IndexVector>
K1 MultiOrthogonalMaxSearchTree<K1, K2, V, UIntSize, NumGroups, GroupSize, KeyVector1, KeyVector2, IndexVector>::key1_at(UIntSize i) const {
    return key1[i];
}

template<typename K1, typename K2, typename V, typename UIntSize, size_t NumGroups, size_t GroupSize, class KeyVector1, class KeyVector2, class IndexVector>
K2 MultiOrthogonalMaxSearchTree<K1, K2, V, UIntSize, NumGroups, GroupSize, KeyVector1, KeyVector2, IndexVector>::key2_at(UIntSize i) const {
    return key2[i];
}

template<typename K1, typename K2, typename V, typename UIntSize, size_t NumGroups, size_t GroupSize, class KeyVector1, class KeyVector2, class IndexVector>
V MultiOrthogonalMaxSearchTree<K1, K2, V, UIntSize, NumGroups, GroupSize, KeyVector1, KeyVector2, IndexVector>::value_at(UIntSize i, size_t group, size_t channel) const {
    return value_of(i, group, channel);
}

template<typename K1, typename K2, typename V, typename UIntSize, size_t NumGroups, size_t GroupSize, class KeyVector1, class KeyVector2, class IndexVector>
size_t MultiOrthogonalMaxSearchTree<K1, K2, V, UIntSize, NumGroups, GroupSize, KeyVector1, KeyVector2, IndexVector>::size() const {
    return key1.size();
}

template<typename K1, typename K2, typename V, typename UIntSize, size_t NumGroups, size_t GroupSize, class KeyVector1, class KeyVector2, class IndexVector>
bool MultiOrthogonalMaxSearchTree<K1, K2, V, UIntSize, NumGroups, GroupSize, KeyVector1, KeyVector2, IndexVector>::empty() const {
    return key1.empty();
}

template<typename K1, typename K2, typename V, typename UIntSize, size_t NumGroups, size_t GroupSize, class KeyVector1, class KeyVector2, class IndexVector>
size_t MultiOrthogonalMaxSearchTree<K1, K2, V, UIntSize, NumGroups, GroupSize, KeyVector1, KeyVector2, IndexVector>::memory_size() const {
    size_t mem_size = get_vec_memory_size(key1) + get_vec_memory_size(key2) + get_vec_memory_size(value);
    mem_size += cross_levels.memory_size();
    return mem_size;
}

}

#endif /* centrolign_multi_orthogonal_max_search_tree_hpp */
