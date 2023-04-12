#include <cstdio>
#include <cstdlib>
#include <cstdlib>

#include <iostream>
#include <vector>
#include <algorithm>
#include <random>

#include "centrolign/max_search_tree.hpp"

using namespace std;
using namespace centrolign;

const pair<int, int> null(-1000000, -1000000);



vector<pair<int, pair<int, int>>> random_key_value_pairs(size_t size, default_random_engine& gen) {
    
    uniform_int_distribution<int> distr(0, size-1);
    
    // ensure unique keys to disambiguate range maxes
    vector<pair<int, pair<int, int>>> key_val_pairs;
    for (int i = 0; i < size; ++i) {
        int key = distr(gen);
        int val = distr(gen);
        key_val_pairs.emplace_back(key, make_pair(val, i));
    }
    return key_val_pairs;
}

vector<pair<int, pair<int, int>>> random_updates(size_t size, size_t number, default_random_engine& gen) {
    
    uniform_int_distribution<int> distr(0, size-1);
    
    vector<pair<int, pair<int, int>>> updates(number);
    for (auto& update : updates) {
        update.first = distr(gen);
        update.second.first = distr(gen);
    }
    return updates;
}

pair<size_t, size_t> equal_range(vector<pair<int, pair<int, int>>>& sorted,
                                 int key) {
    
    size_t begin = -1;
    for (size_t i = 0; i < sorted.size(); ++i) {
        if (sorted[i].first == key) {
            begin = i;
            break;
        }
    }
    if (begin != -1) {
        size_t end = begin + 1;
        while (end < sorted.size() && sorted[end].first == key) {
            ++end;
        }
        return make_pair(begin, end);
    }
    return pair<size_t, size_t>(sorted.size(), sorted.size());
}

void print_kv_pairs(vector<pair<int, pair<int, int>>>& key_val_pairs) {
    for (auto kv : key_val_pairs) {
        cerr << kv.first << ": " << kv.second.first << ", " << kv.second.second << '\n';
    }
}

pair<int, int> range_max(vector<pair<int, pair<int, int>>>& key_val_pairs,
                         int lo, int hi) {

    auto m = null;
    for (size_t i = 0; i < key_val_pairs.size(); ++i) {
        if (key_val_pairs[i].first >= lo && key_val_pairs[i].first < hi) {
            for (size_t j = i; j < key_val_pairs.size(); ++j) {
                if (key_val_pairs[j].first >= hi) {
                    break;
                }
                if (m == null || m < key_val_pairs[j].second) {
                    m = key_val_pairs[j].second;
                }
            }
            break;
        }
    }
    return m;
}

bool test_queries(MaxSearchTree<int, pair<int, int>>& tree,
                  vector<pair<int, pair<int, int>>>& key_val_pairs) {
    
    
    size_t i = 0;
    for (auto val : tree) {
        if (val != key_val_pairs[i]) {
            cerr << "failed iteration test on iteration " << i << "\n";
            return false;
        }
        ++i;
    }
    assert(i == key_val_pairs.size());
    
    for (int k = -1; k < (int) key_val_pairs.size() + 1; ++k) {
        auto vec_range = equal_range(key_val_pairs, k);
        auto tree_range = tree.equal_range(k);
        size_t i = 0;
        for (auto it = tree_range.first; it != tree_range.second; ++it) {
            if (*it != key_val_pairs[vec_range.first + i]) {
                cerr << "failed equal range test on key " << k << "\n";
                return false;
            }
            ++i;
        }
        if (i != vec_range.second - vec_range.first) {
            cerr << "incomplete equal range test on key " << k << " for data:\n";
            print_kv_pairs(key_val_pairs);
            exit(1);
        }
    }
    
    for (int k_lo = -1; k_lo < (int) key_val_pairs.size() + 1; ++k_lo) {
        for (int k_hi = k_lo - 1; k_hi < (int) key_val_pairs.size() + 1; ++k_hi) {
            auto vec_rm = range_max(key_val_pairs, k_lo, k_hi);
            auto tree_rm = tree.range_max(k_lo, k_hi);
            if (vec_rm == null) {
                if (tree_rm != tree.end()) {
                    cerr << "non null range max on range " << k_lo << " " << k_hi << "\n";
                    return false;
                }
            }
            else if ((*tree_rm).second != vec_rm) {
                cerr << "incorrect range max on range " << k_lo << " " << k_hi << " \n";
                return false;
            }
        }
    }
    return true;
}

void do_test(vector<pair<int, pair<int, int>>>& key_val_pairs,
             vector<pair<int, pair<int, int>>>& updates) {
    
    
    MaxSearchTree<int, pair<int, int>> tree(key_val_pairs);
    auto copy = key_val_pairs;
    if (!test_queries(tree, key_val_pairs)) {
        cerr << "without updates on data:\n";
        print_kv_pairs(copy);
        exit(1);
    }
    
    for (auto& update : updates) {
        auto it = tree.begin();
        for (int i = 0; i < update.first; ++i) {
            ++it;
        }
        tree.update(it, update.second);
        key_val_pairs[update.first].second = update.second;
        if (!test_queries(tree, key_val_pairs)) {
            cerr << "after update " << update.first << " -> " << update.second.first << ", " << update.second.second << " in data:\n";
            print_kv_pairs(copy);
            cerr << "with updates:\n";
            print_kv_pairs(updates);
            exit(1);
        }
    }
}

int main(int argc, char* argv[]) {
    
    {
        vector<pair<int, pair<int, int>>> kv_pairs{
            {0, {2, 0}},
            {0, {1, 1}},
            {2, {2, 2}}
        };
        vector<pair<int, pair<int, int>>> updates{
            {1, {1, 0}},
            {1, {1, 0}},
            {2, {1, 0}},
            {2, {0, 0}},
            {1, {0, 0}}
        };
        do_test(kv_pairs, updates);
    }
    {
        vector<pair<int, pair<int, int>>> kv_pairs{
            {0, {1, 2}},
            {0, {0, 3}},
            {1, {1, 2}},
            {2, {7, 0}},
            {6, {1, 3}}
        };
        vector<pair<int, pair<int, int>>> updates{
            {1, {4, 1}},
            {0, {2, 7}},
            {1, {6, 1}},
            {3, {5, 3}},
            {2, {2, 5}}
        };
        do_test(kv_pairs, updates);
    }

    
    random_device rd;
    default_random_engine gen(rd());
    
    uniform_int_distribution<size_t> size_distr(1, 16);
    uniform_int_distribution<size_t> num_updates_distr(1, 20);
    size_t num_tests = 100;
    for (size_t i = 0; i < num_tests; ++i) {
        size_t size = size_distr(gen);
        size_t num_updates = num_updates_distr(gen);
        auto key_val_pairs = random_key_value_pairs(size, gen);
        auto updates = random_updates(size, num_updates, gen);
        do_test(key_val_pairs, updates);
    }
    
    
    cerr << "passed all tests!" << endl;
}
