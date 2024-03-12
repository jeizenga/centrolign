#include "centrolign/utility.hpp"

#include <unistd.h>
#include <sys/resource.h>
#include <stdexcept>
#include <string>
#include <sstream>
#include <iomanip>

#ifdef __APPLE__
#include <mach/mach.h>
#endif

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

bool parse_bool(const std::string& str) {
    bool parsed = false;
    if (str == "true" || str == "True" || str == "TRUE" || str == "T" || str == "t" || str == "1") {
        parsed = true;
    }
    else if (str == "false" || str == "False" || str == "FALSE" || str == "F" || str == "f" || str == "0") {
        parsed = false;
    }
    else {
        throw runtime_error("Could not parse \"" + str + "\" as a Boolean value");
    }
    return parsed;
}

vector<string> tokenize(const string& str, char delim) {
    vector<string> tokens;
    size_t prev = 0;
    size_t cursor = 0;
    while (prev < str.size()) {
        while (cursor < str.size() && str[cursor] != delim) {
            ++cursor;
        }
        tokens.emplace_back(str.substr(prev, cursor - prev));
        prev = cursor + 1;
        cursor = prev;
    }
    return tokens;
}

std::string join(const std::vector<std::string>& values, const std::string& sep) {
    
    size_t size = 0;
    for (const auto& value : values) {
        size += value.size();
    }
    if (!values.empty()) {
        size += (values.size() - 1) * sep.size();
    }
    
    std::string joined(size, '\0');
    size_t i = 0;
    for (size_t j = 0; j < values.size(); ++j) {
        if (j != 0) {
            for (size_t k = 0; k < sep.size(); ++k) {
                joined[i++] = sep[k];
            }
        }
        const auto& value = values[j];
        for (size_t k = 0; k < value.size(); ++k) {
            joined[i++] = value[k];
        }
    }
    
    return joined;
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

int64_t current_memory_usage() {
#if defined(__APPLE__)
    // from https://en.wikichip.org/wiki/resident_set_size
    struct mach_task_basic_info info;
    mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
    if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, (task_info_t)&info, &count) == KERN_SUCCESS) {
        return info.resident_size;
    }
    return -1; /* query failed */
#elif defined(__linux__)
    string procfp = "/proc/" + to_string(getpid()) + "/stat";
    string rss;
    // this "file" is updated about 1x per second, and it seems that it's sometimes unavailable
    for (int attempt = 0; attempt < 3 && rss.empty(); ++attempt) {
        if (attempt != 0) {
            usleep(50);
        }
        ifstream proc(procfp);
        if (!proc) {
            continue;
        }
        // RSS is the 24th item, according to man PROC(5)
        for (int i = 0; i < 23 && proc.good(); ++i) {
            proc.ignore(numeric_limits<streamsize>::max(), ' ');
        }
        if (proc) {
            getline(proc, rss, ' ');
        }
    }
    if (rss.empty()) {
        return -1;
    }
    else {
        return parse_int(rss) * getpagesize();
    }
#else
    cerr << "Warning: Failed to measure current memory usage. Active memory monitoring only supported on Linux and Apple systems\n";
    return -1;
#endif
}

int64_t max_memory_usage() {
    
    struct rusage usage;
    int code = getrusage(RUSAGE_SELF, &usage);
    if (code != 0) {
        // failed to measure memory usage
        return -1;
    }
    else {
        int64_t max_mem = usage.ru_maxrss;
        // seems that mac and linux do this differently
#ifdef __linux__
        max_mem *= 1024;
#endif
        return max_mem;
    }
}

std::string format_memory_usage(int64_t mem) {
    if (mem < 0) {
        return "N/A";
    }
    std::string unit = "";
    double memd = mem;
    if (memd >= 1024.0) {
        memd /= 1024.0;
        unit = "k";
        if (memd >= 1024.0) {
            memd /= 1024.0;
            unit = "M";
            if (memd >= 1024.0) {
                memd /= 1024.0;
                unit = "G";
            }
            if (memd >= 1024.0) {
                memd /= 1024.0;
                unit = "T";
            }
        }
    }
    
    std::stringstream strm;
    strm << fixed << setprecision(2) << memd << ' ' << unit << 'B';
    return strm.str();
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
