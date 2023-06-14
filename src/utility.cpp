#include "centrolign/utility.hpp"

#include <stdexcept>
#include <string>

namespace centrolign {

using namespace std;


vector<pair<string, string>> parse_fasta(istream& in) {
    
    vector<pair<string, string>> parsed;
    
    string line;
    size_t prev_line_len = -1;
    size_t prev_prev_line_len = -1;
    uint64_t line_num = 0;
    while (in) {
        getline(in, line);
        ++line_num;
        if (line.front() == '>') {
            size_t space_pos = line.find(' ');
            parsed.emplace_back();
            parsed.back().first = line.substr(1, space_pos - 1);
            if (parsed.back().first.empty()) {
                throw runtime_error("FASTA input is missing sequence name at line " + to_string(line_num));
            }
            
            prev_line_len = -1;
            prev_prev_line_len = -1;
        }
        else {
            if (parsed.empty()) {
                throw runtime_error("FASTA input does not have sequence name line");
            }
            else if (prev_prev_line_len != -1 && prev_line_len != prev_prev_line_len
                     && !line.empty()) {
                throw runtime_error("Encountered sequence lines of unequal lengths that were not followed by a sequence name at line " + to_string(line_num) + " of FASTA input");
            }
            else if (prev_line_len != -1 && line.size() > prev_line_len) {
                throw runtime_error("Encountered adjacent sequence lines of increasing lengths at line " + to_string(line_num) + " of FASTA input");
            }
            
            
            parsed.back().second.append(line);
            
            prev_prev_line_len = prev_line_len;
            prev_line_len = line.size();
        }
    }
    if (parsed.empty()) {
        throw runtime_error("FASTA file is empty");
    }
    return parsed;
}

int64_t parse_int(const std::string& str) {
    size_t idx;
    int64_t parsed = stoll(str, &idx);
    if (idx != str.size()) {
        throw runtime_error("Could not parse \"" + str + "\" as an integer");
    }
    return parsed;
}

double parse_double(const std::string& str) {
    size_t idx;
    double parsed = stod(str, &idx);
    if (idx != str.size()) {
        throw runtime_error("Could not parse \"" + str + "\" as a decimal number");
    }
    return parsed;
}

string path_to_string(const BaseGraph& graph, const vector<uint64_t>& path) {
    string seq;
    for (auto i : path) {
        seq.push_back(graph.label(i));
    }
    return seq;
}

istream* get_input(const string& stream_name, ifstream& openable) {
    if (stream_name == "-") {
        return &cin;
    }
    else {
        openable.open(stream_name);
        if (!openable) {
            throw runtime_error("Could not open file " + stream_name);
        }
        return &openable;
    }
}

vector<size_t> range_vector(size_t size) {
    vector<size_t> range(size, 0);
    for (size_t i = 1; i < range.size(); ++i) {
        range[i] = i;
    }
    return range;
}

// modified from sdsl-lite by simon gog
// using built-in method or
// 64-bit version of 32-bit proposal of
// http://www-graphics.stanford.edu/~seander/bithacks.html
const uint32_t lt_hi[] = {
    0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7
};
uint32_t hi_bit(uint64_t x)
{
#ifdef __GNUG__
    if (x == 0)
        return 0;
    return 63 - __builtin_clzll(x);
#else
    uint64_t t,tt; // temporaries
    if ((tt = x >> 32)) { // hi >= 32
        if ((t = tt >> 16)) { // hi >= 48
            return (tt = t >> 8) ? 56 + lt_hi[tt] : 48 + lt_hi[t];
        } else { // hi < 48
            return (t = tt >> 8) ? 40 + lt_hi[t] : 32 + lt_hi[tt];
        }
    } else { // hi < 32
        if ((t = x >> 16)) { // hi >= 16
            return (tt = t >> 8) ? 24 + lt_hi[tt] : 16 + lt_hi[t];
        } else { // hi < 16
            return (tt = x >> 8) ?  8 + lt_hi[tt] : lt_hi[x];
        }
    }
#endif
}

const uint8_t enc[] = {
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, // 16
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, // 32
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, // 48
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, // 64
    5, 0, 5, 1, 5, 5, 5, 2, 5, 5, 5, 5, 5, 5, 4, 5, // 80
    5, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, // 96
    5, 0, 5, 1, 5, 5, 5, 2, 5, 5, 5, 5, 5, 5, 4, 5, // 112
    5, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, // 128
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5
};

uint8_t encode_base(char nt) {
    return enc[(uint8_t)nt];
}

std::string encode_seq(const std::string& seq) {
    std::string encoded(seq.size(), 0);
    for (size_t i = 0; i < seq.size(); ++i) {
        encoded[i] = encode_base(seq[i]);
    }
    return encoded;
}

const char dec[] = {'A', 'C', 'G', 'T', 'N', '\0'};

char decode_base(uint8_t nt) {
    return dec[nt];
}

std::string decode_seq(const std::string& seq) {
    std::string decoded(seq.size(), '\0');
    for (size_t i = 0; i < seq.size(); ++i) {
        decoded[i] = decode_base(seq[i]);
    }
    return decoded;
}

}
