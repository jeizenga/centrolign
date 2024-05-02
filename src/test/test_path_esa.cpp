#include <cstdio>
#include <cstdlib>
#include <cstdlib>

#include <iostream>
#include <vector>
#include <algorithm>
#include <random>

#include "centrolign/path_esa.hpp"
#include "centrolign/utility.hpp"
#include "centrolign/test_util.hpp"

using namespace std;
using namespace centrolign;

class TestPathESA : public PathESA {
public:
    using PathESA::construct_suffix_array;
    using PathESA::construct_lcp_array;
};

vector<size_t> direct_suffix_array(const string& str) {
    
    auto indexes = range_vector(str.size());
    
    sort(indexes.begin(), indexes.end(), [&](size_t i, size_t j) {
        return str.substr(i, str.size() - i) < str.substr(j, str.size() - j);
    });
    
    return indexes;
}

vector<size_t> direct_lcp_array(const string& str, const vector<size_t>& suffix_array) {
    
    vector<size_t> lcp_array(suffix_array.size(), 0);
    
    for (size_t i = 1; i < lcp_array.size(); ++i) {
        auto p1 = suffix_array[i];
        auto p2 = suffix_array[i - 1];
        size_t l = 0;
        while (str[p1 + l] == str[p2 + l]) {
            ++l;
        }
        lcp_array[i] = l;
    }
    
    return lcp_array;
}

void test_suffix_array_construction(const string& str) {
    
    auto direct = direct_suffix_array(str);
    auto sais = TestPathESA::construct_suffix_array(str);
    
    if (sais != direct) {
        cerr << "SA construction failed on input: " << str << '\n';
        cerr << "direct:\n";
        for (auto i : direct) {
            cerr << '\t' << i << "\t" << str.substr(i, str.size()) << '\n';
        }
        cerr << "sais:\n";
        for (auto i : sais) {
            cerr << '\t' << i << "\t" << str.substr(i, str.size()) << '\n';
        }
        exit(1);
    }
    
    auto direct_lcp = direct_lcp_array(str, sais);
    auto inverse_sa = invert(sais);
    auto kasai = TestPathESA::construct_lcp_array(str, sais, inverse_sa);
    
    if (direct_lcp != kasai) {
        cerr << "LCP construction failed on input: " << str << '\n';
        cerr << "direct:\n";
        for (size_t i = 0; i < str.size(); ++i) {
            cerr << i << '\t' << direct_lcp[i] << "\t" << str.substr(sais[i], str.size()) << '\n';
        }
        cerr << "kasai:\n";
        for (size_t i = 0; i < str.size(); ++i) {
            cerr << i << '\t' << kasai[i] << "\t" << str.substr(sais[i], str.size()) << '\n';
        }
        exit(1);
    }
}


int main(int argc, char* argv[]) {
    
    random_device rd;
    default_random_engine gen(rd());
    
    {
        // example from SAIS paper
        string str = "MMIISSIISSIIPPII$";
        test_suffix_array_construction(str);
    }
    
    {
        int num_reps = 10;
        vector<int> sizes{1, 2, 3, 6, 10, 20, 40, 80};
        for (auto s : sizes) {
            for (int rep = 0; rep < num_reps; ++rep) {
                {
                    auto seq = random_sequence(s, gen);
                    seq.push_back('$');
                    test_suffix_array_construction(seq);
                }
                {
                    auto seq = random_low_entropy_sequence(s, gen);
                    seq.push_back('$');
                    test_suffix_array_construction(seq);
                }
            }
        }
    }
    
    
    cerr << "passed all tests!" << endl;
}
