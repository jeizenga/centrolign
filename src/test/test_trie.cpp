#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <unordered_set>
#include <vector>
#include <list>
#include <string>
#include <random>
#include <iostream>

#include "centrolign/trie.hpp"

using namespace std;
using namespace centrolign;


int main(int argc, char* argv[]) {
     
    random_device rd;
    default_random_engine gen(rd());
    
    {
        // hackily switch this from a sequence-based test to integer one
        Trie trie;
        assert(trie.node_size() == 1);
        auto path1 = trie.insert_sequence("1", vector<uint64_t>{'A', 'A', 'C'});
        assert(trie.node_size() == 4);
        auto path2 = trie.insert_sequence("2", vector<uint64_t>{'A', 'T', 'T'});
        assert(trie.node_size() == 6);
        auto path3 = trie.insert_sequence("3", vector<uint64_t>{'A', 'T', 'C'});
        assert(trie.node_size() == 7);
        
        for (uint64_t node_id = 0; node_id < trie.node_size(); ++node_id) {
            for (auto c : string("ACGT")) {
                auto next_id = trie.follow(node_id, c);
                if (next_id != -1) {
                    assert(trie.label(next_id) == c);
                }
            }
        }
        
        for (uint64_t n = 0; n < trie.node_size(); ++n) {
            auto prev = trie.previous(n);
            assert(prev.size() == trie.previous_size(n));
            assert(prev.size() <= 1);
            if (!prev.empty()) {
                auto p = prev.front();
                auto next = trie.next(p);
                auto it = find(next.begin(), next.end(), n);
                assert(it != next.end());
            }
        }
        
        assert(trie.path_size() == 3);
        
        uint64_t n1, n2, n3, n4, n5, n6, n7;
        n1 = trie.get_root();
        assert(trie.next_size(n1) == 1);
        assert(trie.previous_size(n1) == 0);
        n2 = trie.follow(n1, 'A');
        assert(n2 != -1);
        assert(trie.next_size(n2) == 2);
        assert(trie.previous_size(n2) == 1);
        n3 = trie.follow(n2, 'A');
        assert(n3 != -1);
        assert(trie.next_size(n3) == 1);
        assert(trie.previous_size(n3) == 1);
        n4 = trie.follow(n3, 'C');
        assert(n4 != -1);
        assert(trie.next_size(n4) == 0);
        assert(trie.previous_size(n4) == 1);
        n5 = trie.follow(n2, 'T');
        assert(n5 != -1);
        assert(trie.next_size(n5) == 2);
        assert(trie.previous_size(n5) == 1);
        n6 = trie.follow(n5, 'T');
        assert(n6 != -1);
        assert(trie.next_size(n6) == 0);
        assert(trie.previous_size(n6) == 1);
        n7 = trie.follow(n5, 'C');
        assert(n7 != -1);
        assert(trie.next_size(n7) == 0);
        assert(trie.previous_size(n7) == 1);
        
        vector<uint64_t> p1{n2, n3, n4};
        vector<uint64_t> p2{n2, n5, n6};
        vector<uint64_t> p3{n2, n5, n7};
        
        assert(trie.path(path1) == p1);
        assert(trie.path(path2) == p2);
        assert(trie.path(path3) == p3);
    }
    
    cerr << "passed all tests!" << endl;
}
