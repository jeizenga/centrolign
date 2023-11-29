//
//  compare_truth_aln.cpp
//  
//
//  Created by Jordan Eizenga on 11/21/23.
//

#include <tuple>
#include <utility>
#include <fstream>
#include <iostream>
#include <string>
#include <cstdint>
#include <cctype>
#include <stdio.h>

#include "centrolign/utility.hpp"

using namespace std;
using namespace centrolign;

void print_help() {
    cerr << "\n";
    cerr << "Usage: compare_truth_aln seq1_identity.txt seq2_identity.txt truth_cigar.txt aln_cigar.txt\n";
}

vector<size_t> parse_identity(ifstream& in) {
    
    vector<size_t> identity;
    
    string line;
    while (in) {
        getline(in, line);
        if (!line.empty()) {
            identity.push_back(parse_int(line));
        }
    }
    
    return identity;
}

vector<pair<int, char>> parse_cigar(ifstream& in) {
    
    string cigar_str;
    string line;
    while (in) {
        getline(in, line);
        cigar_str += line;
    }
    
    vector<pair<int, char>> parsed;
    for (size_t i = 0; i < cigar_str.size(); ++i) {
        size_t begin = i;
        while (isdigit(cigar_str[i])) {
            ++i;
        }
        parsed.emplace_back(parse_int(cigar_str.substr(begin, i - begin)), cigar_str[i]);
    }
    return parsed;
}

// tuple of (matches, seq len 1, seq len 2)
tuple<size_t, size_t, size_t> compute_consistency(const vector<size_t>& identity1,
                                                  const vector<size_t>& identity2,
                                                  const vector<pair<int, char>>& cigar) {
    
    size_t matches = 0;
    size_t len1 = 0;
    size_t len2 = 0;
    
    size_t i = 0;
    size_t j = 0;
    for (const auto& cigar_op : cigar) {
        switch (cigar_op.second) {
            case 'M':
            case 'X':
            case '=':
                for (int k = 0; k < cigar_op.second; ++k) {
                    if (identity1[i + k] == identity2[j + k]) {
                        ++matches;
                    }
                }
                len1 += cigar_op.first;
                len2 += cigar_op.first;
                break;
            case 'I':
            case 'H':
            case 'S':
                len2 += cigar_op.first;
                break;
            case 'D':
            case 'N':
                len1 += cigar_op.first;
                break;
            default:
                cerr << "error: unrecognized cigar operation " << cigar_op.second << '\n';
                exit(1);
                break;
        }
    }
    
    return make_tuple(matches, len1, len2);
}

int main(int argc, char* argv[]) {

    if (argc != 5) {
        print_help();
        return 1;
    }
    
    string ident_filename1 = argv[1];
    string ident_filename2 = argv[2];
    string truth_filename = argv[3];
    string aln_filename = argv[4];
    
    ifstream ident_in1(ident_filename1);
    ifstream ident_in2(ident_filename2);
    ifstream truth_in(truth_filename);
    ifstream aln_in(aln_filename);
    
    if (!ident_in1) {
        cerr << "error: could not open " << ident_filename1 << '\n';
        return 1;
    }
    if (!ident_in2) {
        cerr << "error: could not open " << ident_filename2 << '\n';
        return 1;
    }
    if (!truth_in) {
        cerr << "error: could not open " << truth_filename << '\n';
        return 1;
    }
    if (!aln_in) {
        cerr << "error: could not open " << aln_filename << '\n';
        return 1;
    }
    
    auto identity1 = parse_identity(ident_in1);
    auto identity2 = parse_identity(ident_in2);
    
    auto truth_cigar = parse_cigar(truth_in);
    auto aln_cigar = parse_cigar(aln_in);
    
    size_t truth_matches, aln_matches, len1, len2;
    tie(truth_matches, len1, len2) = compute_consistency(identity1, identity2, truth_cigar);
    tie(aln_matches, len1, len2) = compute_consistency(identity1, identity2, aln_cigar);
    
    cout << "truth matches: " << truth_matches << '\n';
    cout << "truth match rate: " << double(2 * truth_matches) / double(len1 + len2) << '\n';
    cout << "aln matches: " << aln_matches << '\n';
    cout << "aln match rate: " << double(2 * aln_matches) / double(len1 + len2) << '\n';
    
    return 0;
}

