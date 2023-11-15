#include <vector>
#include <string>
#include <list>
#include <iostream>
#include <cstdint>
#include <cassert>
#include <regex>
#include <random>
#include <cmath>
#include <algorithm>
#include <limits>
#include <getopt.h>

#include "centrolign/utility.hpp"
#include "centrolign/alignment.hpp"

using namespace std;
using namespace centrolign;

const int64_t default_num_generations = 100;
constexpr double default_hor_indel_rate = 1.0 / 5000000.0;
constexpr double default_exp_hor_indel = 2.0;
constexpr double default_monomer_indel_rate = 1.0 / 25000000.0;
constexpr double default_exp_monomer_indel = 3.0;
constexpr double default_point_indel_rate = 1.0 / 50000.0;
constexpr double default_exp_point_indel = 1.5;
constexpr double default_subs_rate = 0.0;
// TODO: length distribution parameters

// from the CentromereArchitect paper supplementary material
const string alpha_consensus = "AATCTGCAAGTGGACATTTGGAGCGCTTTGAGGCCTATGGTGGAAAAGGAAATATCTTCACATAAAAACTAGACAGAAGCATTCTC"
                               "AGAAACTTCTTTGTGATGTGTGCATTCAACTCACAGAGTTGAACCTTTCTTTTGATAGAGCAGTTTTGAAACACTCTTTTTGTAG";


void print_help() {
    
    cerr << "\n";
    cerr << "Usage: sim_centromere sequence.fasta annotation_bed\n";
    cerr << "\n";
    cerr << "Options:\n";
    cerr << " --generations / -g INT            Number of generations [" << default_num_generations << "]\n";
    cerr << " --hor-indel-rate / -h FLOAT       Rate of full HOR indels [" << default_hor_indel_rate << "]\n";
    cerr << " --hor-indel-size / -H FLOAT       Expected size of full HOR indels in HOR units [" << default_exp_hor_indel << "]\n";
    cerr << " --monomer-indel-rate / -m FLOAT   Rate of full monomer indels [" << default_monomer_indel_rate << "]\n";
    cerr << " --monomer-indel-size / -M FLOAT   Expected size of full monomer indels in monomer units [" << default_exp_monomer_indel << "]\n";
    cerr << " --point-indel-rate / -p FLOAT     Rate of point indels [" << default_point_indel_rate << "]\n";
    cerr << " --point-indel-size / -P FLOAT     Expected size of point indels in base pairs [" << default_exp_point_indel << "]\n";
    cerr << " --substitution-rate / -s FLOAT    Rate of substitutions [" << default_subs_rate << "]\n";
    cerr << " --seed / -z INT                   Seed for PRNG [random]\n";
    cerr << " --help                            Print this message and exit\n";
}

struct EvolvedBase {
    
    EvolvedBase(char base, size_t origin, size_t idx_in_monomer, size_t monomer_idx) :
        base(base), origin(origin), idx_in_monomer(idx_in_monomer), monomer_idx(monomer_idx) { }
    
    EvolvedBase(const EvolvedBase& other) = default;
    EvolvedBase& operator=(const EvolvedBase& other) = default;
    
    char base = '\0';
    size_t origin = -1;
    size_t idx_in_monomer = -1;
    size_t monomer_idx = -1;
};

vector<tuple<string, size_t, size_t, string>> parse_bed(istream& in) {
    
    static const string header_start = "track name";
    
    vector<tuple<string, size_t, size_t, string>> intervals;
    while (in) {
        string line;
        getline(in, line);
        if (line.substr(0, header_start.size()) == header_start) {
            continue;
        }
        auto tokens = tokenize(line);
        assert(tokens.size() >= 4);
        intervals.emplace_back(tokens[0], parse_int(tokens[1]), parse_int(tokens[2]), tokens[3]);
    }
    
    return intervals;
}

size_t parse_alpha_type(const string& type) {
    
    static const regex type_regex("S\\d+C[XYM0-9/]+H\\d-?[dLAB].(\\d+)");
    smatch sub;
    bool found_match = regex_search(type, sub, type_regex);
    if (!found_match || sub.size() < 2) {
        return -1;
    }
    else {
        string monomer_num = sub[1];
        return parse_int(monomer_num);
    }
}

list<EvolvedBase> initialize_sequence(const string& fasta, const string& bed) {
    
    ifstream fasta_in, bed_in;
    auto seqs = parse_fasta(*get_input(fasta, fasta_in));
    assert(seqs.size() == 1);
    auto& name = seqs.front().first;
    auto& seq = seqs.front().second;
    auto intervals = parse_bed(*get_input(bed, bed_in));
    for (const auto& interval : intervals) {
        assert(get<0>(interval) == name);
        assert(get<2>(interval) <= seq.size());
    }
    sort(intervals.begin(), intervals.end());
    
    size_t last_monomer = 0;
    
    list<EvolvedBase> annotated_seq;
    size_t seq_idx = 0;
    size_t interval_idx = 0;
    while (seq_idx < seq.size()) {
        
        // add any sequence between monomers
        size_t next_interval_begin = interval_idx < intervals.size() ? get<1>(intervals[interval_idx]) : seq.size();
        while (seq_idx < next_interval_begin) {
            annotated_seq.emplace_back(seq[seq_idx], seq_idx, last_monomer, alpha_consensus.size());
            ++seq_idx;
        }
        
        if (interval_idx < intervals.size()) {
            // add the sequence corresponding to a monomer
            const auto& interval = intervals[interval_idx++];
            
            // FIXME: i need to figure out which strand the HORs are on and possibly take the
            // reverse complement of the alpha consensus
            
            string monomer = seq.substr(get<1>(interval), get<2>(interval) - get<1>(interval));
            size_t monomer_type = parse_alpha_type(get<3>(interval));
            auto alpha_aln = align_ond(monomer, alpha_consensus);
            
            // add the sequence annotated with the positions
            size_t cons_pos = 0;
            for (size_t i = 0; i < alpha_aln.size(); ++i) {
                if (alpha_aln[i].node_id2 != AlignedPair::gap) {
                    cons_pos = alpha_aln[i].node_id2;
                }
                if (alpha_aln[i].node_id1 == AlignedPair::gap) {
                    continue;
                }
                annotated_seq.emplace_back(monomer[i], seq_idx + i, monomer_type, cons_pos);
            }
            if (monomer_type != -1) {
                last_monomer = monomer_type;
            }
            seq_idx += monomer.size();
        }
    }
    
    return annotated_seq;
}

// return pair of (monomers increasing in index, HOR size in monomers)
pair<bool, size_t> determine_hor(const list<EvolvedBase>& sequence) {
    
    // figure out the size the HOR
    size_t max_monomer = 0;
    size_t min_monomer = 0;
    for (const auto& base : sequence) {
        if (base.monomer_idx != -1) {
            max_monomer = max(base.monomer_idx, max_monomer);
            min_monomer = min(base.monomer_idx, min_monomer);
        }
    }
    
    // count the number of times the index goes up or down by 1
    size_t prev_monomer = -1;
    uint64_t count_increasing = 0;
    uint64_t count_decreasing = 0;
    for (const auto& base : sequence) {
        if (base.monomer_idx != -1) {
            if (prev_monomer != -1) {
                if (prev_monomer == base.monomer_idx - 1 ||
                    (prev_monomer == max_monomer && base.monomer_idx == min_monomer)) {
                    ++count_increasing;
                }
                else if (prev_monomer == base.monomer_idx + 1 ||
                         (prev_monomer == min_monomer && base.monomer_idx == max_monomer)) {
                    ++count_decreasing;
                }
            }
            prev_monomer = base.monomer_idx;
        }
    }
    
    return make_pair(count_increasing > count_decreasing, max_monomer);
}

void duplicate_subseq(list<EvolvedBase>& sequence,
                      list<EvolvedBase>::iterator begin, list<EvolvedBase>::iterator end) {
    
    for (auto it = begin; it != end; ++it) {
        sequence.insert(begin, *it);
    }
}

template<class Generator>
void point_insert(list<EvolvedBase>& sequence, list<EvolvedBase>::iterator pos,
                  size_t size, Generator& gen) {
    
    static const string alphabet = "ACGT";
    
    auto it = pos;
    if (it == sequence.end()) {
        // we don't allow indels to overhang the end
        return;
    }
    for (size_t i = 0; i < size; ++i) {
        char base = alphabet[uniform_int_distribution<int>(0, 3)(gen)];
        sequence.emplace(pos, base, it->origin, it->idx_in_monomer, it->monomer_idx);
    }
}

// TODO: do i really need this?
void delete_subseq(list<EvolvedBase>& sequence,
                   list<EvolvedBase>::iterator begin, list<EvolvedBase>::iterator end) {
    
    sequence.erase(begin, end);
}

template<class Generator>
char substitute(char base, Generator& gen) {
    static const string alphabet = "ACGT";
    char sub = base;
    while (sub == base) {
        sub = alphabet[uniform_int_distribution<int>(0, 3)(gen)];
    }
    return sub;
}

template<class Generator>
bool sample_prob(double probability, Generator& gen) {
    return uniform_real_distribution<double>(0.0, 1.0)(gen) < probability;
}

template <class Generator>
uint64_t sample_geom(double mean, bool from_0, Generator& gen) {
    double mu;
    if (from_0) {
        assert(mean >= 0.0);
        mu = mean;
    }
    else {
        assert(mean >= 1.0);
        mu = mean - 1.0;
    }
    
    if (mu == 0.0) {
        return from_0 ? 0 : 1;
    }
    
    double lam = log((mu + 1.0) / mu);
    double expo = exponential_distribution<double>(lam)(gen);
    uint64_t geom = (uint64_t) expo;
    return from_0 ? geom : geom + 1;
}

list<EvolvedBase>::iterator advance(const list<EvolvedBase>& sequence, list<EvolvedBase>::iterator it, size_t distance) {
    auto advanced = it;
    for (size_t i = 0; i < distance; ++i) {
        ++advanced;
        if (advanced == sequence.end()) {
            break;
        }
    }
    return advanced;
}


template<class Generator>
list<EvolvedBase>::iterator advance_monomers(const list<EvolvedBase>& sequence, list<EvolvedBase>::iterator pos,
                                             size_t num_monomers, Generator& gen) {
    
    assert(num_monomers > 0);
    
    // walk until you pass this position N monomers forward
    size_t num_mons_passed = 0;
    size_t prev_idx = -1;
    auto it = pos;
    for (; it != sequence.end(); ++it) {
        if (prev_idx != -1 && prev_idx > it->idx_in_monomer) {
            ++num_mons_passed;
        }
        if ((num_mons_passed == num_monomers && it->idx_in_monomer >= pos->idx_in_monomer) ||
            num_mons_passed > num_monomers) {
            break;
        }
        prev_idx = it->idx_in_monomer;
    }
    
    // find the set of bases aligned to the equivalent position in the source monomer
    vector<list<EvolvedBase>::iterator> equal_positions;
    auto prev_it = it;
    while (prev_it != pos) {
        --prev_it;
        if (prev_it == pos) {
            break;
        }
        else if (prev_it->idx_in_monomer == pos->idx_in_monomer) {
            equal_positions.push_back(prev_it);
        }
        else if (prev_it->idx_in_monomer < pos->idx_in_monomer) {
            break;
        }
    }
    
    if (equal_positions.empty()) {
        // there are no exactly aligned bases, so we just take an adjacent one
        return it;
    }
    else {
        // choose randomly among the aligned bases
        return equal_positions[uniform_int_distribution<size_t>(0, equal_positions.size() - 1)(gen)];
    }
}


template<class Generator>
list<EvolvedBase>::iterator advance_hors(const list<EvolvedBase>& sequence, list<EvolvedBase>::iterator pos,
                                         size_t num_hors, size_t hor_size, bool increasing, Generator& gen) {
    
    assert(pos->monomer_idx != -1);
    assert(num_hors != 0);
    
    // loop to walk to the
    size_t num_hors_passed = 0;
    size_t prev_idx = -1;
    size_t prev_monomer = -1;
    auto it = pos;
    list<EvolvedBase>::iterator final_hor_begin;
    for (; it != sequence.end(); ++it) {
        if (prev_monomer != -1 && it->monomer_idx != -1 &&
            (prev_monomer != it->monomer_idx || (prev_monomer == it->monomer_idx && prev_idx > it->idx_in_monomer))) {
            // we've moved onto a new monomer
            
            // compute distance in the direction of increasing monomer number
            int64_t fwd_dist, rev_dist;
            if (prev_monomer < it->monomer_idx) {
                fwd_dist = it->monomer_idx - prev_monomer;
                rev_dist = num_hors - it->monomer_idx + prev_monomer;
            }
            else {
                fwd_dist = num_hors - prev_monomer + it->monomer_idx;
                rev_dist = prev_monomer - it->monomer_idx;
            }
            
            if (increasing) {
                // we expect the increasing monomer direction to be smaller
                if (fwd_dist <= rev_dist) {
                    // we don't think this was a backward stutter
                    if ((prev_monomer < it->monomer_idx && pos->monomer_idx > prev_monomer && pos->monomer_idx <= it->monomer_idx) ||
                        (it->monomer_idx < prev_monomer && (pos->monomer_idx > prev_monomer || pos->monomer_idx <= it->monomer_idx))) {
                        // we either landed on the source monomer or overshot it
                        ++num_hors_passed;
                        if (num_hors_passed == num_hors) {
                            final_hor_begin = it;
                        }
                        else if (num_hors_passed > num_hors) {
                            break;
                        }
                    }
                }
            }
            else {
                // we expect the decreasing monomer direction to be smaller
                if (rev_dist <= fwd_dist) {
                    // we don't think this was a backward stutter
                    if ((prev_monomer > it->monomer_idx && pos->monomer_idx >= it->monomer_idx && pos->monomer_idx < prev_monomer) ||
                        (it->monomer_idx > prev_monomer && (pos->monomer_idx >= it->monomer_idx || pos->monomer_idx < prev_monomer))) {
                        // we either landed on the source monomer or overshot it
                        ++num_hors_passed;
                        if (num_hors_passed == num_hors) {
                            final_hor_begin = it;
                        }
                        else if (num_hors_passed > num_hors) {
                            break;
                        }
                    }
                }
            }
        }
        
        if (it->monomer_idx != -1) {
            prev_monomer = it->monomer_idx;
        }
    }
    
    // get a parse of the HOR into monomers
    auto final_hor_end = it;
    vector<list<EvolvedBase>::iterator> monomer_begins;
    prev_idx = -1;
    for (auto it = final_hor_begin; it != final_hor_end; ++it) {
        if (prev_idx != -1 && prev_idx > it->idx_in_monomer) {
            monomer_begins.push_back(it);
        }
        prev_idx = it->idx_in_monomer;
    }
    
    // find the monomer(s) of the correct family in this HOR
    vector<pair<list<EvolvedBase>::iterator, list<EvolvedBase>::iterator>> candidate_monomers;
    for (size_t i = 0; i < monomer_begins.size(); ++i) {
        if (monomer_begins[i]->monomer_idx == pos->monomer_idx) {
            candidate_monomers.emplace_back(monomer_begins[i],
                                            i + 1 < monomer_begins.size() ? monomer_begins[i + 1] : final_hor_end);
        }
    }
    
    if (candidate_monomers.empty()) {
        // there are no monomers of the correct family, find the closest ones there are
        size_t closest_idx = -1;
        size_t closest_dist = -1;
        for (size_t i = 0; i < monomer_begins.size(); ++i) {
            if (monomer_begins[i]->monomer_idx != -1) {
                size_t dist;
                // compute distance in the direction of increasing monomer number
                // TODO: repetitive with above
                size_t fwd_dist, rev_dist;
                if (monomer_begins[i]->monomer_idx < pos->monomer_idx) {
                    fwd_dist = pos->monomer_idx - monomer_begins[i]->monomer_idx;
                    rev_dist = num_hors - pos->monomer_idx + monomer_begins[i]->monomer_idx;
                }
                else {
                    fwd_dist = num_hors - monomer_begins[i]->monomer_idx + pos->monomer_idx;
                    rev_dist = monomer_begins[i]->monomer_idx - pos->monomer_idx;
                }
                dist = min(fwd_dist, rev_dist);
                if (closest_idx == -1 || dist < closest_dist) {
                    closest_idx = i;
                    closest_dist = dist;
                }
            }
        }
        
        // to handle strange cases where large parts of a HOR are deleted
        if (candidate_monomers.empty()) {
            closest_idx = 0;
        }
        
        // expand the window to any neighboring monomers with no HOR type classification
        size_t i = closest_idx;
        size_t j = closest_idx + 1;
        while (i != 0 && monomer_begins[i - 1]->monomer_idx == -1) {
            --i;
        }
        while (j != monomer_begins.size() && monomer_begins[j]->monomer_idx == -1) {
            ++j;
        }
        
        for (; i < j; ++i) {
            candidate_monomers.emplace_back(monomer_begins[i],
                                            i + 1 < monomer_begins.size() ? monomer_begins[i + 1] : final_hor_end);
        }
    }
    
    assert(!candidate_monomers.empty());
    
    // choose a random monomer
    auto monomer_range = candidate_monomers[uniform_int_distribution<size_t>(0, candidate_monomers.size() - 1)(gen)];
    
    // find the closest to in-register bases in in the chosen monomer
    vector<list<EvolvedBase>::iterator> candidate_bases;
    int64_t closest_base_dist = numeric_limits<int64_t>::max();
    for (auto it = monomer_range.first; it != monomer_range.second; ++it) {
        int64_t dist = abs(int64_t(it->idx_in_monomer - pos->idx_in_monomer));
        if (dist <= closest_base_dist) {
            // this is the closest we've seen so far (or tied)
            if (dist < closest_base_dist) {
                // this is closer than the previous, so we get rid of them
                candidate_bases.clear();
                closest_base_dist = dist;
            }
            candidate_bases.push_back(it);
        }
    }
    
    assert(!candidate_bases.empty());
    
    // choose a random base
    return candidate_bases[uniform_int_distribution<size_t>(0, candidate_bases.size() - 1)(gen)];
}

int main(int argc, char* argv[]) {
    
    int64_t num_generations = default_num_generations;
    
    double hor_indel_rate = default_hor_indel_rate;
    double monomer_indel_rate = default_monomer_indel_rate;
    double point_indel_rate = default_point_indel_rate;
    double subs_rate = default_subs_rate;
    
    double exp_hor_indel = default_exp_hor_indel;
    double exp_monomer_indel = default_exp_monomer_indel;
    double exp_point_indel = default_exp_point_indel;
    
    uint64_t seed = -1;
    bool set_seed = false;
    
    static const int help_opt = 1000;
    while (true)
    {
        static struct option options[] = {
            {"generations", required_argument, NULL, 'g'},
            {"hor-indel-rate", required_argument, NULL, 'h'},
            {"hor-indel-size", required_argument, NULL, 'H'},
            {"monomer-indel-rate", required_argument, NULL, 'm'},
            {"monomer-indel-size", required_argument, NULL, 'M'},
            {"point-indel-rate", required_argument, NULL, 'p'},
            {"point-indel-size", required_argument, NULL, 'P'},
            {"substitution-rate", required_argument, NULL, 's'},
            {"seed", required_argument, NULL, 'z'},
            {"help", no_argument, NULL, help_opt},
            {NULL, 0, NULL, 0}
        };
        int o = getopt_long(argc, argv, "g:h:H:m:M:p:P:s:z:", options, NULL);
        
        if (o == -1) {
            // end of uptions
            break;
        }
        switch (o)
        {
            case 'g':
                num_generations = parse_int(optarg);
                break;
            case 'h':
                hor_indel_rate = parse_double(optarg);
                break;
            case 'H':
                exp_hor_indel = parse_double(optarg);
                break;
            case 'm':
                monomer_indel_rate = parse_double(optarg);
                break;
            case 'M':
                exp_hor_indel = parse_double(optarg);
                break;
            case 'p':
                point_indel_rate = parse_double(optarg);
                break;
            case 'P':
                exp_point_indel = parse_double(optarg);
                break;
            case 's':
                subs_rate = parse_double(optarg);
                break;
            case 'z':
                seed = parse_int(optarg);
                set_seed = true;
                break;
            case help_opt:
                print_help();
                return 0;
            default:
                print_help();
                return 1;
        }
    }
    
    if (argc - optind != 2) {
        cerr << "error: expected 2 positional arguments but got " << (argc - optind) << "\n";
        print_help();
        return 1;
    }
    
    string fasta = argv[optind++];
    string bed = argv[optind];
        
    // make two copies of the sequence
    auto seq1 = initialize_sequence(fasta, bed);
    auto seq2 = seq1;
    
    if (!set_seed) {
        random_device rd;
        seed = rd();
    }
    
    mt19937 gen(seed);
    
    bool monomers_increasing;
    size_t hor_size;
    tie(monomers_increasing, hor_size) = determine_hor(seq1);
    
    for (size_t generation = 0; generation < num_generations; ++generation) {
        for (auto seq_ptr : {&seq1, &seq2}) {
            auto& seq = *seq_ptr;
            
            // add HOR indels
            for (auto it = seq.begin(); it != seq.end();) {
                // only sample them at monomers that can be in HOR register
                if (it->monomer_idx != -1 && sample_prob(hor_indel_rate, gen)) {
                    size_t size = sample_geom(exp_hor_indel, false, gen);
                    auto indel_end = advance_hors(seq, it, size, hor_size, monomers_increasing, gen);
                    if (indel_end == seq.end()) {
                        // we don't allow indels to hang over the edge
                        ++it;
                    }
                    else if (sample_prob(0.5, gen)) {
                        // a HOR insertion
                        duplicate_subseq(seq, it, indel_end);
                        ++it;
                    }
                    else {
                        // a HOR deletion
                        delete_subseq(seq, it, indel_end);
                        it = indel_end;
                    }
                }
                else {
                    ++it;
                }
            }
            
            // add monomer indels
            for (auto it = seq.begin(); it != seq.end();) {
                if (sample_prob(monomer_indel_rate, gen)) {
                    size_t size = sample_geom(exp_monomer_indel, false, gen);
                    auto indel_end = advance_monomers(seq, it, size, gen);
                    if (indel_end == seq.end()) {
                        // we don't allow indels to hang over the edge
                        ++it;
                    }
                    else if (sample_prob(0.5, gen)) {
                        // a monomer insertion
                        duplicate_subseq(seq, it, indel_end);
                        ++it;
                    }
                    else {
                        // a monomer deletion
                        delete_subseq(seq, it, indel_end);
                        it = indel_end;
                    }
                }
                else {
                    ++it;
                }
            }

            // add point indels
            for (auto it = seq.begin(); it != seq.end();) {
                if (sample_prob(point_indel_rate, gen)) {
                    size_t size = sample_geom(exp_point_indel, false, gen);
                    if (sample_prob(0.5, gen)) {
                        // a point insertion
                        point_insert(seq, it, size, gen);
                        ++it;
                    }
                    else {
                        // a point deletion
                        auto indel_end = advance(seq, it, size);
                        if (indel_end != seq.end()) {
                            delete_subseq(seq, it, indel_end);
                            it = indel_end;
                        }
                        else {
                            ++it;
                        }
                    }
                }
                else {
                    ++it;
                }
            }
            
            // add substitutions
            for (auto it = seq.begin(); it != seq.end(); ++it) {
                if (sample_prob(subs_rate, gen)) {
                    it->base = substitute(it->base, gen);
                }
            }
        }
    }
}
