#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <unordered_set>
#include <vector>
#include <list>
#include <string>
#include <random>

#include "centrolign/trie.hpp"

using namespace std;
using namespace centrolign;


int main(int argc, char* argv[]) {
     
    random_device rd;
    default_random_engine gen(rd());
    
    {
        Trie trie;
        assert(trie.node_size() == 1);
        auto path1 = trie.insert_sequence("1", "AAC");
        assert(trie.node_size() == 4);
        auto path2 = trie.insert_sequence("2", "ATT");
        assert(trie.node_size() == 6);
        auto path3 = trie.insert_sequence("3", "ATC");
        assert(trie.node_size() == 7);
        
        for (uint64_t node_id = 0; node_id < trie.node_size(); ++node_id) {
            for (auto c : string("ACGT")) {
                auto next_id = trie.follow(node_id, c);
                if (next_id != -1) {
                    assert(trie.label(next_id) == c);
                }
            }
        }
        
        assert(trie.path_size() == 3);
        
        uint64_t n1, n2, n3, n4, n5, n6, n7;
        n1 = trie.get_root();
        assert(trie.children_size(n1) == 1);
        n2 = trie.follow(n1, 'A');
        assert(n2 != -1);
        assert(trie.children_size(n2) == 2);
        n3 = trie.follow(n2, 'A');
        assert(n3 != -1);
        assert(trie.children_size(n3) == 1);
        n4 = trie.follow(n3, 'C');
        assert(n4 != -1);
        assert(trie.children_size(n4) == 0);
        n5 = trie.follow(n2, 'T');
        assert(n5 != -1);
        assert(trie.children_size(n5) == 2);
        n6 = trie.follow(n5, 'T');
        assert(n6 != -1);
        assert(trie.children_size(n6) == 0);
        n7 = trie.follow(n5, 'C');
        assert(n7 != -1);
        assert(trie.children_size(n7) == 0);
        
        vector<uint64_t> p1{n2, n3, n4};
        vector<uint64_t> p2{n2, n5, n6};
        vector<uint64_t> p3{n2, n5, n7};
        
        assert(trie.path(path1) == p1);
        assert(trie.path(path2) == p2);
        assert(trie.path(path3) == p3);
    }
    
    cerr << "passed all tests!" << endl;
}
