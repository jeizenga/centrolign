#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include <random>
#include <limits>
#include <set>
#include <utility>

#include "centrolign/multi_orthogonal_max_search_tree.hpp"
#include "centrolign/packed_vector.hpp"
#include "centrolign/vector_support.hpp"

using namespace std;
using namespace centrolign;

using key1_t = pair<int32_t, uint32_t>;

const size_t NumGroups = 2;
const size_t GroupSize = 3;
const size_t NumChannels = NumGroups * GroupSize;

using UnpackedTree = MultiOrthogonalMaxSearchTree<key1_t, uint32_t, float, uint32_t, NumGroups, GroupSize>;
using PackedTree = MultiOrthogonalMaxSearchTree<key1_t, uint32_t, float, uint32_t, NumGroups, GroupSize,
                                                VectorPair<SignedPackedVector, PackedVector, int32_t, uint32_t>, PackedVector, PackedVector>;

const float mininf = numeric_limits<float>::lowest();

struct Entry {
    key1_t key1;
    uint32_t key2;
    array<float, NumChannels> values;
};

template<class Tree>
void check_query(const Tree& tree, const vector<Entry>& entries,
                 const key1_t& lo1, const key1_t& hi1, uint32_t hi2, size_t group) {

    auto result = tree.range_max(lo1, hi1, hi2, group);
    for (size_t c = 0; c < GroupSize; ++c) {
        // brute force
        bool found = false;
        float best = mininf;
        for (const auto& entry : entries) {
            if (entry.key1 >= lo1 && entry.key1 < hi1 && entry.key2 < hi2) {
                if (!found || entry.values[group * GroupSize + c] > best) {
                    best = entry.values[group * GroupSize + c];
                    found = true;
                }
            }
        }

        bool tree_found = (result[c].second != tree.size());
        if (found != tree_found ||
            (found && (result[c].first != best ||
                       tree.value_at(result[c].second, group, c) != best ||
                       tree.key1_at(result[c].second) < lo1 || !(tree.key1_at(result[c].second) < hi1) ||
                       !(tree.key2_at(result[c].second) < hi2)))) {
            cerr << "failed query [(" << lo1.first << "," << lo1.second << "), (" << hi1.first << "," << hi1.second << ")) x [0, " << hi2 << ") on group " << group << " channel " << c << '\n';
            cerr << "expected found? " << found << ", value " << best << '\n';
            cerr << "got found? " << tree_found << ", value " << result[c].first << '\n';
            cerr << "entries:\n";
            for (const auto& entry : entries) {
                cerr << "\t(" << entry.key1.first << "," << entry.key1.second << ") " << entry.key2;
                for (auto v : entry.values) {
                    cerr << ' ' << v;
                }
                cerr << '\n';
            }
            exit(1);
        }
    }
}

template<class Tree>
void do_test(size_t size, size_t key_range, size_t num_updates, size_t num_queries, mt19937& gen) {

    uniform_int_distribution<int32_t> shift_distr(-int32_t(key_range), key_range);
    uniform_int_distribution<uint32_t> id_distr(0, 3);
    uniform_int_distribution<uint32_t> offset_distr(0, key_range);
    uniform_real_distribution<float> value_distr(-100.0, 100.0);

    // make unique key pairs
    set<pair<key1_t, uint32_t>> keys;
    while (keys.size() < size) {
        keys.emplace(key1_t(shift_distr(gen), id_distr(gen)), offset_distr(gen));
    }
    vector<pair<key1_t, uint32_t>> data(keys.begin(), keys.end());
    shuffle(data.begin(), data.end(), gen);

    vector<Entry> entries;
    for (const auto& item : data) {
        entries.emplace_back();
        entries.back().key1 = item.first;
        entries.back().key2 = item.second;
        entries.back().values.fill(mininf);
    }

    Tree tree(data, mininf);
    if (tree.size() != size) {
        cerr << "wrong size " << tree.size() << " for " << size << " entries\n";
        exit(1);
    }

    for (const auto& entry : entries) {
        auto i = tree.find(entry.key1, entry.key2);
        if (i == tree.size() || tree.key1_at(i) != entry.key1 || tree.key2_at(i) != entry.key2) {
            cerr << "failed to find entry (" << entry.key1.first << "," << entry.key1.second << ") " << entry.key2 << '\n';
            exit(1);
        }
    }

    auto random_query = [&]() {
        size_t group = gen() % NumGroups;
        uint32_t hi2 = offset_distr(gen) + 1;
        key1_t k1(shift_distr(gen), id_distr(gen));
        key1_t k2(shift_distr(gen), id_distr(gen));
        switch (gen() % 3) {
            case 0:
                // left half-line, like the anchorer
                check_query(tree, entries, key1_t(numeric_limits<int32_t>::min(), 0), k1, hi2, group);
                break;
            case 1:
                // right half-line, like the anchorer
                check_query(tree, entries, k1, key1_t(numeric_limits<int32_t>::max(), numeric_limits<uint32_t>::max()), hi2, group);
                break;
            default:
                check_query(tree, entries, min(k1, k2), max(k1, k2), hi2, group);
                break;
        }
    };

    for (size_t u = 0; u < num_updates && size != 0; ++u) {
        auto& entry = entries[gen() % entries.size()];
        array<float, NumChannels> values;
        for (size_t ch = 0; ch < NumChannels; ++ch) {
            values[ch] = value_distr(gen);
            // values that don't increase are ignored
            entry.values[ch] = max(entry.values[ch], values[ch]);
        }
        tree.update(tree.find(entry.key1, entry.key2), values);
        for (size_t q = 0; q < num_queries; ++q) {
            random_query();
        }
    }
    // also check a fresh tree with no updates
    for (size_t q = 0; q < num_queries; ++q) {
        random_query();
    }
}

int main(int argc, char* argv[]) {

    mt19937 gen(1037);

    for (size_t size : {0, 1, 2, 3, 7, 8, 9, 15, 16, 17, 31, 64, 100, 257}) {
        for (size_t key_range : {2, 10, 1000}) {
            if (key_range * key_range * 4 < size) {
                continue;
            }
            for (size_t rep = 0; rep < 5; ++rep) {
                do_test<UnpackedTree>(size, key_range, 2 * size + 5, 4, gen);
                do_test<PackedTree>(size, key_range, 2 * size + 5, 4, gen);
            }
        }
    }

    // a few larger trees with deeper cross trees
    for (size_t size : {1000, 3000}) {
        do_test<UnpackedTree>(size, 500, 3 * size, 1, gen);
        do_test<PackedTree>(size, 500, 3 * size, 1, gen);
    }

    cerr << "passed all tests!" << endl;
    return 0;
}
