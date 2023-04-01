#include "centrolign/utility.hpp"

#include <stdexcept>

namespace centrolign {

using namespace std;


vector<pair<string, string>> parse_fasta(istream& in) {
    
    vector<pair<string, string>> parsed;
    
    string line;
    while (in) {
        getline(in, line);
        if (line.front() == '>') {
            size_t space_pos = line.find(' ');
            parsed.emplace_back();
            parsed.back().first = line.substr(1, space_pos - 1);
            if (parsed.back().first.empty()) {
                throw runtime_error("FASTA file is missing sequence name");
            }
        }
        else {
            if (parsed.empty()) {
                throw runtime_error("FASTA file is does not have sequence name line");
            }
            // TODO: this is harder to check in the middle of a file...
            parsed.back().second.append(line);
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
