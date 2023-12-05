#include <vector>
#include <string>
#include <list>
#include <iostream>
#include <sstream>
#include <cstdint>
#include <cassert>
#include <stdexcept>
#include <regex>
#include <random>
#include <cmath>
#include <algorithm>
#include <limits>
#include <getopt.h>

#include "centrolign/utility.hpp"
#include "centrolign/alignment.hpp"
#include "centrolign/tree.hpp"

using namespace std;
using namespace centrolign;

const int64_t default_num_generations = 100;
constexpr double default_small_hor_indel_rate = 1.0 / 1000000.0;
constexpr double default_exp_small_hor_indel = 1.25;
constexpr double default_large_hor_indel_rate = 1.0 / 8000000.0;
constexpr double default_exp_large_hor_indel = 8.0;
constexpr double default_large_hor_indel_tail_heaviness = 10.0;
constexpr double default_monomer_indel_rate = 1.0 / 25000000.0;
constexpr double default_exp_monomer_indel = 3.0;
constexpr double default_point_indel_rate = 1.0 / 1000000.0;
constexpr double default_exp_point_indel = 1.5;
constexpr double default_subs_rate = 1.0 / 100000.0;
// TODO: length distribution parameters

static const bool find_opt_alignment = true;

// from the CentromereArchitect paper supplementary material
const string alpha_consensus = "AATCTGCAAGTGGACATTTGGAGCGCTTTGAGGCCTATGGTGGAAAAGGAAATATCTTCACATAAAAACTAGACAGAAGCATTCTCAGAAACTTCTTTGTGATGTGTGCATTCAACTCACAGAGTTGAACCTTTCTTTTGATAGAGCAGTTTTGAAACACTCTTTTTGTAG";


void print_help() {
    cerr << "\n";
    cerr << "Usage: sim_centromere [options] sequence.fasta annotation_bed\n";
    cerr << "\n";
    cerr << "Options:\n";
    cerr << " --output / -o PREFIX               Prefix for output (required)\n";
    cerr << " --generations / -g INT             Number of generations [" << default_num_generations << "]\n";
    cerr << " --tree / -T FILE                   Tree in Newick format (overrides --generations)\n";
    cerr << " --hor-indel-small-rate / -h FLOAT  Rate of small full HOR indels per base per generation [" << default_small_hor_indel_rate << "]\n";
    cerr << " --hor-indel-small-size / -H FLOAT  Expected size of small full HOR indels in HOR units [" << default_exp_small_hor_indel << "]\n";
    cerr << " --hor-indel-large-rate / -r FLOAT  Rate of heavy-tailed full HOR indels per base per generation [" << default_large_hor_indel_rate << "]\n";
    cerr << " --hor-indel-large-size / -R FLOAT  Expected size of heavy-tailed full HOR indels in HOR units [" << default_exp_large_hor_indel << "]\n";
    cerr << " --hor-indel-heaviness / -t FLOAT   Tail heaviness of HOR indel distribution [" << default_large_hor_indel_tail_heaviness << "]\n";
    cerr << " --monomer-indel-rate / -m FLOAT    Rate of full monomer indels per base per generation [" << default_monomer_indel_rate << "]\n";
    cerr << " --monomer-indel-size / -M FLOAT    Expected size of full monomer indels in monomer units [" << default_exp_monomer_indel << "]\n";
    cerr << " --point-indel-rate / -p FLOAT      Rate of point indels per base per generation [" << default_point_indel_rate << "]\n";
    cerr << " --point-indel-size / -P FLOAT      Expected size of point indels in base pairs [" << default_exp_point_indel << "]\n";
    cerr << " --substitution-rate / -s FLOAT     Rate of substitutions per base per generation [" << default_subs_rate << "]\n";
    cerr << " --seed / -z INT                    Seed for PRNG [random]\n";
    cerr << " --help                             Print this message and exit\n";
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
        if (line.empty() || line.substr(0, header_start.size()) == header_start) {
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
    
    static const bool debug = false;
    
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
    
    size_t last_monomer = -1;
    
    list<EvolvedBase> annotated_seq;
    size_t seq_idx = 0;
    size_t interval_idx = 0;
    while (seq_idx < seq.size()) {
        
        // add any sequence between monomers
        size_t next_interval_begin = interval_idx < intervals.size() ? get<1>(intervals[interval_idx]) : seq.size();
        while (seq_idx < next_interval_begin) {
            if (debug) {
                cerr << "add base between monomers, seq idx " << seq_idx << ", last mon " << last_monomer << '\n';
            }
            annotated_seq.emplace_back(seq[seq_idx], seq_idx, alpha_consensus.size(), last_monomer);
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
            if (debug) {
                cerr << "adding monomer of type " << monomer_type << '\n';
            }
            
            // add the sequence annotated with the positions
            size_t cons_pos = 0;
            for (size_t i = 0; i < alpha_aln.size(); ++i) {
                if (alpha_aln[i].node_id2 != AlignedPair::gap) {
                    cons_pos = alpha_aln[i].node_id2;
                }
                if (alpha_aln[i].node_id1 == AlignedPair::gap) {
                    continue;
                }
                if (debug) {
                    cerr << "add monomer base " << (seq_idx + alpha_aln[i].node_id1) << ", mon type " << monomer_type << ", cons pos " << cons_pos  << '\n';
                }
                annotated_seq.emplace_back(monomer[alpha_aln[i].node_id1], seq_idx + alpha_aln[i].node_id1,
                                           cons_pos, monomer_type);
            }
            if (monomer_type != -1) {
                last_monomer = monomer_type;
            }
            seq_idx += monomer.size();
        }
    }
    
    return annotated_seq;
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

/*
 * Hurwitz zeta function adapted from scipy source code:
 * https://github.com/scipy/scipy/blob/f3cda9db9ad67cc897953fc2d1786edc651c1d4e/scipy/special/cephes/zeta.c
 */
double zeta(double x, double q)
{
    
    static double A[] = {
        12.0,
        -720.0,
        30240.0,
        -1209600.0,
        47900160.0,
        -1.8924375803183791606e9,    /*1.307674368e12/691 */
        7.47242496e10,
        -2.950130727918164224e12,    /*1.067062284288e16/3617 */
        1.1646782814350067249e14,    /*5.109094217170944e18/43867 */
        -4.5979787224074726105e15,    /*8.028576626982912e20/174611 */
        1.8152105401943546773e17,    /*1.5511210043330985984e23/854513 */
        -7.1661652561756670113e18    /*1.6938241367317436694528e27/236364091 */
    };
    
    /* 30 Nov 86 -- error in third coefficient fixed */
    
    int i;
    double a, b, k, s, t, w;
    
    if (x == 1.0)
        return numeric_limits<double>::infinity();
    
    if (x < 1.0) {
        return std::nan("0");
    }
    
    if (q <= 0.0) {
        if (q == floor(q)) {
            return numeric_limits<double>::infinity();
        }
        if (x != floor(x))
            return std::nan("0");    /* because q^-x not defined */
    }
    
    /* Asymptotic expansion
     * https://dlmf.nist.gov/25.11#E43
     */
    if (q > 1e8) {
        return (1/(x - 1) + 1/(2*q)) * pow(q, 1 - x);
    }
    
    /* Euler-Maclaurin summation formula */
    
    /* Permit negative q but continue sum until n+q > +9 .
     * This case should be handled by a reflection formula.
     * If q<0 and x is an integer, there is a relation to
     * the polyGamma function.
     */
    s = pow(q, -x);
    a = q;
    i = 0;
    b = 0.0;
    bool done = false;
    while (!done && ((i < 9) || (a <= 9.0))) {
        i += 1;
        a += 1.0;
        b = pow(a, -x);
        s += b;
        if (fabs(b / s) < numeric_limits<double>::epsilon())
            done = true;
    }
    
    if (!done) {
        w = a;
        s += b * w / (x - 1.0);
        s -= 0.5 * b;
        a = 1.0;
        k = 0.0;
        for (i = 0; i < 12; i++) {
            a *= x + k;
            b /= w;
            t = a * b / A[i];
            s = s + t;
            t = fabs(t / s);
            if (t < numeric_limits<double>::epsilon())
                break;
            k += 1.0;
            a *= x + k;
            b /= w;
            k += 1.0;
        }
    }
    return (s);
}

// root finding by bisection
template<class Func>
double bisect(const Func& func, double lower_bound, double upper_bound, double tol) {
    
    double hi = upper_bound;
    double lo = lower_bound;
    double fhi = func(upper_bound);
    double flo = func(lower_bound);
    
    assert((fhi > 0.0) != (flo > 0.0));
    if (fhi == 0.0) {
        return hi;
    }
    if (flo == 0.0) {
        return lo;
    }
    
    while (abs(hi - lo) >= tol) {
        double mid = (hi + lo) / 2.0;
        double fmid = func(mid);
        if ((fmid > 0.0) == (flo > 0.0)) {
            lo = mid;
            flo = fmid;
        }
        else {
            hi = mid;
            fhi = fmid;
        }
    }
    
    return (hi + lo) / 2.0;
}

double discrete_pareto_expected_value(double beta, double sigma) {
    return pow(sigma, beta) * zeta(beta, sigma);
}

uint64_t discrete_pareto_quantile(double q, double beta, double sigma) {
    double q_term = pow(1.0 - q, 1.0 / beta);
    return ceil(sigma * (1.0 - q_term) / q_term);
}

template<class Generator>
uint64_t sample_discrete_pareto(double beta, double sigma, Generator& gen) {
    return discrete_pareto_quantile(uniform_real_distribution<double>(0.0, 1.0)(gen), beta, sigma);
}

double choose_discrete_pareto_sigma(double expected_val, double beta) {
    assert(expected_val > 1.0);
    assert(beta > 1.0);
    
    auto f = [&](double s) {
        if (s == 0.0) {
            return 1.0 - expected_val;
        }
        else {
            return discrete_pareto_expected_value(beta, s) - expected_val;
        }
    };
    
    double hi_sigma = 1.0;
    double exp_val;
    while ((exp_val = discrete_pareto_expected_value(beta, hi_sigma)) < expected_val || isnan(exp_val)) {
        hi_sigma *= 2.0;
    }
    
    return bisect(f, 0.0, hi_sigma, 1e-6);
}






class Evolver {
public:
    Evolver() = default;
    ~Evolver() = default;
    
    double small_hor_indel_rate = default_small_hor_indel_rate;
    double large_hor_indel_rate = default_large_hor_indel_rate;
    double monomer_indel_rate = default_monomer_indel_rate;
    double point_indel_rate = default_point_indel_rate;
    double subs_rate = default_subs_rate;
    
    double exp_small_hor_indel = default_exp_small_hor_indel;
    double exp_monomer_indel = default_exp_monomer_indel;
    double exp_point_indel = default_exp_point_indel;
    
    double large_hor_indel_beta = 1.5;
    double large_hor_indel_sigma = 5.0;
    
    // determine sampling parameters for the sequence
    void determine_hor(const list<EvolvedBase>& sequence);
    
    template<class Generator>
    list<EvolvedBase> evolve(const list<EvolvedBase>& parent, uint64_t num_generations,
                             Generator& gen) const;
    
    
private:
    
    bool monomers_increasing = true;
    size_t hor_size = -1;
    
    template<class Generator>
    list<EvolvedBase>::iterator advance_hors(list<EvolvedBase>& sequence, list<EvolvedBase>::iterator pos,
                                             size_t num_hors, Generator& gen) const;
    
    template<class Generator>
    list<EvolvedBase>::iterator advance_monomers(const list<EvolvedBase>& sequence, list<EvolvedBase>::iterator pos,
                                                 size_t num_monomers, Generator& gen) const;
    
    list<EvolvedBase>::iterator advance(const list<EvolvedBase>& sequence, list<EvolvedBase>::iterator it, size_t distance) const;

    void duplicate_subseq(list<EvolvedBase>& sequence,
                          list<EvolvedBase>::iterator begin, list<EvolvedBase>::iterator end) const;
    
    template<class Generator>
    void point_insert(list<EvolvedBase>& sequence, list<EvolvedBase>::iterator pos,
                      size_t size, Generator& gen) const;
    
    void delete_subseq(list<EvolvedBase>& sequence,
                       list<EvolvedBase>::iterator begin, list<EvolvedBase>::iterator end) const;
};

void Evolver::determine_hor(const list<EvolvedBase>& sequence) {
    
    static const bool debug = false;
    
    // figure out the size the HOR
    size_t max_monomer = 0;
    size_t min_monomer = -1;
    for (const auto& base : sequence) {
        if (base.monomer_idx != -1) {
            max_monomer = max(base.monomer_idx, max_monomer);
            min_monomer = min(base.monomer_idx, min_monomer);
        }
    }
    
    if (debug) {
        cerr << "min and max monomers: " << min_monomer << ", " << max_monomer << '\n';
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
    
    if (debug) {
        cerr << "num increasing " << count_increasing << ", num decreasing " << count_decreasing << '\n';
    }
    
    monomers_increasing = (count_increasing > count_decreasing);
    hor_size = max_monomer - min_monomer + 1;
}


template<class Generator>
list<EvolvedBase> Evolver::evolve(const list<EvolvedBase>& parent, uint64_t num_generations,
                                  Generator& gen) const {
    
    if (hor_size == -1) {
        throw std::runtime_error("must determine HOR size before evolving");
    }
    
    auto seq = parent;
    for (size_t generation = 1; generation <= num_generations; ++generation) {
        if (generation % 10 == 0) {
            cerr << "generation " << generation << " of " << num_generations << '\n';
        }
        
        // add small HOR indels
        for (auto it = seq.begin(); it != seq.end();) {
            // only sample them at monomers that can be in HOR register
            if (it->monomer_idx != -1 && sample_prob(small_hor_indel_rate, gen)) {
                size_t size = sample_geom(exp_small_hor_indel, false, gen);
//                ++num_hor_indels;
//                size_hor_indels += size;
                auto indel_end = advance_hors(seq, it, size, gen);
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
        
        // add large HOR indels
        for (auto it = seq.begin(); it != seq.end();) {
            // only sample them at monomers that can be in HOR register
            if (it->monomer_idx != -1 && sample_prob(large_hor_indel_rate, gen)) {
                size_t size = sample_discrete_pareto(large_hor_indel_beta, large_hor_indel_sigma, gen);
//                ++num_hor_indels;
//                size_hor_indels += size;
                auto indel_end = advance_hors(seq, it, size, gen);
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
//                ++num_mon_indels;
//                size_mon_indels += size;
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
//                ++num_point_indels;
//                size_point_indels += size;
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
//                ++num_substitutions;
                it->base = substitute(it->base, gen);
            }
        }
        
    }
    
    return seq;
}

template<class Generator>
list<EvolvedBase>::iterator Evolver::advance_hors(list<EvolvedBase>& sequence, list<EvolvedBase>::iterator pos,
                                                  size_t num_hors, Generator& gen) const {
    
    static const bool debug = false;
    
    assert(pos->monomer_idx != -1);
    assert(num_hors != 0);
    
    if (debug) {
        cerr << "advancing " << num_hors << " in a HOR of size " << hor_size << " that is increasing? " << monomers_increasing << "\n";
        cerr << "starting position monomer " << pos->monomer_idx << ", index " << pos->idx_in_monomer << '\n';
    }
    
    // loop to walk to find the range of the HOR that the mutation ends on
    size_t num_hors_passed = 0;
    size_t prev_idx = pos->idx_in_monomer;
    size_t prev_monomer = pos->monomer_idx;
    size_t prev_advancing_monomer = pos->monomer_idx;
    auto it = pos;
    list<EvolvedBase>::iterator final_hor_begin = sequence.end(), final_hor_end = sequence.end();
    for (; it != sequence.end(); ++it) {
        //        if (debug) {
        //            cerr << "searching, seq idx " << it->origin << ", monomer " << it->monomer_idx << ", mon idx " << it->idx_in_monomer << ", prev monomer " << prev_monomer << ", prev mon idx " << prev_idx << '\n';
        //        }
        if (prev_monomer != it->monomer_idx || (prev_monomer == it->monomer_idx && prev_idx > it->idx_in_monomer)) {
            // we've moved onto a new monomer
            if (debug) {
                cerr << "entering new monomer\n";
                cerr << "seq idx " << it->origin << ", monomer " << it->monomer_idx << ", mon idx " << it->idx_in_monomer << ", prev monomer " << prev_monomer << ", prev mon idx " << prev_idx << ", prev advancing monomer " << prev_advancing_monomer << '\n';
            }
            // compute distance in the direction of increasing monomer number
            int64_t fwd_dist, rev_dist;
            if (prev_advancing_monomer < it->monomer_idx) {
                fwd_dist = it->monomer_idx - prev_advancing_monomer;
                rev_dist = hor_size - it->monomer_idx + prev_advancing_monomer;
            }
            else {
                fwd_dist = hor_size - prev_advancing_monomer + it->monomer_idx;
                rev_dist = prev_advancing_monomer - it->monomer_idx;
            }
            
            if (monomers_increasing) {
                // we expect the increasing monomer direction to be smaller
                if (fwd_dist <= rev_dist && fwd_dist > 0) {
                    // we don't think this was a backward stutter
                    if (debug) {
                        cerr << "appears to be a forward direction advancement of dist " << fwd_dist << "\n";
                    }
                    if ((prev_advancing_monomer < it->monomer_idx && pos->monomer_idx > prev_advancing_monomer && pos->monomer_idx <= it->monomer_idx) ||
                        (it->monomer_idx < prev_advancing_monomer && (pos->monomer_idx > prev_advancing_monomer || pos->monomer_idx <= it->monomer_idx))) {
                        // we either landed on the source monomer or overshot it
                        if (debug) {
                            cerr << "entering a new HOR\n";
                        }
                        ++num_hors_passed;
                        if (num_hors_passed == num_hors) {
                            if (debug) {
                                cerr << "this is the start of the final HOR\n";
                            }
                            final_hor_begin = it;
                        }
                        else if (num_hors_passed > num_hors) {
                            if (debug) {
                                cerr << "this is the end of the final HOR\n";
                            }
                            final_hor_end = it;
                            break;
                        }
                    }
                    
                    prev_advancing_monomer = it->monomer_idx;
                }
            }
            else {
                // we expect the decreasing monomer direction to be smaller
                if (rev_dist <= fwd_dist && rev_dist > 0) {
                    // we don't think this was a backward stutter
                    if (debug) {
                        cerr << "appears to be a backward direction advancement of dist " << rev_dist << "\n";
                    }
                    if ((prev_advancing_monomer > it->monomer_idx && pos->monomer_idx >= it->monomer_idx && pos->monomer_idx < prev_advancing_monomer) ||
                        (it->monomer_idx > prev_advancing_monomer && (pos->monomer_idx >= it->monomer_idx || pos->monomer_idx < prev_advancing_monomer))) {
                        // we either landed on the source monomer or overshot it
                        if (debug) {
                            cerr << "entering a new HOR\n";
                        }
                        ++num_hors_passed;
                        if (num_hors_passed == num_hors) {
                            if (debug) {
                                cerr << "this is the start of the final HOR\n";
                            }
                            final_hor_begin = it;
                        }
                        else if (num_hors_passed > num_hors) {
                            if (debug) {
                                cerr << "this is the end of the final HOR\n";
                            }
                            final_hor_end = it;
                            break;
                        }
                    }
                    
                    prev_advancing_monomer = it->monomer_idx;
                }
            }
        }
        
        prev_idx = it->idx_in_monomer;
        if (it->monomer_idx != -1) {
            prev_monomer = it->monomer_idx;
            if (prev_advancing_monomer == -1) {
                prev_advancing_monomer = it->monomer_idx;
            }
        }
    }
    
    if (final_hor_begin == sequence.end()) {
        // we we're able to walk far enough to get to the HOR we were aiming for
        if (debug) {
            cerr << "hit end of sequence before finding the target HOR\n";
        }
        return sequence.end();
    }
    
    
    if (debug) {
        cerr << "final hor is from monomer " << final_hor_begin->monomer_idx << ", idx " << final_hor_begin->idx_in_monomer << ", origin " << final_hor_begin->origin << " -> ";
        if (final_hor_end == sequence.end()) {
            cerr << "END";
        }
        else {
            cerr << "monomer " << final_hor_end->monomer_idx << ", idx " << final_hor_end->idx_in_monomer << ", origin " << final_hor_end->origin;
        }
        cerr << '\n';
    }
    
    // get a parse of the HOR into monomers
    vector<list<EvolvedBase>::iterator> monomer_begins;
    prev_idx = -1;
    for (auto it = final_hor_begin; it != final_hor_end; ++it) {
        if (prev_idx == -1 || prev_idx > it->idx_in_monomer) {
            monomer_begins.push_back(it);
        }
        prev_idx = it->idx_in_monomer;
    }
    
    if (debug) {
        cerr << "parse the hor into " << monomer_begins.size() << " monomers\n";
    }
    
    // find the monomer(s) of the correct family in this HOR
    vector<pair<list<EvolvedBase>::iterator, list<EvolvedBase>::iterator>> candidate_monomers;
    for (size_t i = 0; i < monomer_begins.size(); ++i) {
        if (monomer_begins[i]->monomer_idx == pos->monomer_idx) {
            candidate_monomers.emplace_back(monomer_begins[i],
                                            i + 1 < monomer_begins.size() ? monomer_begins[i + 1] : final_hor_end);
        }
    }
    if (debug) {
        cerr << "got " << candidate_monomers.size() << " candidates that match monomer idx " << pos->monomer_idx << "\n";
    }
    
    if (candidate_monomers.empty()) {
        // there are no monomers of the correct family, find the closest ones there are
        
        if (final_hor_end == sequence.end()) {
            // try to figure out if the reason there is no exact match is because this HOR is truncated
            
            if (monomers_increasing) {
                if (candidate_monomers.front().first->monomer_idx < candidate_monomers.back().first->monomer_idx) {
                    // the walked interval does not wrap
                    if (pos->monomer_idx > candidate_monomers.back().first->monomer_idx
                        || pos->monomer_idx < candidate_monomers.front().first->monomer_idx) {
                        return sequence.end();
                    }
                }
                else {
                    // the walked interval wraps
                    if (pos->monomer_idx > candidate_monomers.front().first->monomer_idx &&
                        pos->monomer_idx < candidate_monomers.back().first->monomer_idx) {
                        return sequence.end();
                    }
                }
            }
            else {
                if (candidate_monomers.front().first->monomer_idx > candidate_monomers.back().first->monomer_idx) {
                    // the walked interval does not wrap
                    if (pos->monomer_idx > candidate_monomers.front().first->monomer_idx ||
                        pos->monomer_idx < candidate_monomers.back().first->monomer_idx) {
                        return sequence.end();
                    }
                }
                else {
                    // the walked interval wraps
                    if (pos->monomer_idx < candidate_monomers.back().first->monomer_idx &&
                        pos->monomer_idx > candidate_monomers.front().first->monomer_idx) {
                        return sequence.end();
                    }
                }
            }
        }
        
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
        
        if (debug) {
            cerr << "after finding closest monomers, got " << candidate_monomers.size() << " candidates\n";
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
    
    if (debug) {
        cerr << "chose candidate monomers:\n";
        for (auto r : candidate_monomers) {
            auto l = r.second;
            --l;
            cerr << '\t' << r.first->origin << "," << r.first->monomer_idx << "," << r.first->idx_in_monomer << " -> " << l->origin << "," << l->monomer_idx << "," << l->idx_in_monomer << '\n';
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
    
    if (debug) {
        cerr << "candidate bases to end advancement:\n";
        for (auto it : candidate_bases) {
            cerr << '\t' << it->origin << ", " << it->monomer_idx << ", " << it->idx_in_monomer << '\n';
        }
    }
    
    // choose a random base
    return candidate_bases[uniform_int_distribution<size_t>(0, candidate_bases.size() - 1)(gen)];
}

template<class Generator>
list<EvolvedBase>::iterator Evolver::advance_monomers(const list<EvolvedBase>& sequence, list<EvolvedBase>::iterator pos,
                                                      size_t num_monomers, Generator& gen) const {
    
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

list<EvolvedBase>::iterator Evolver::advance(const list<EvolvedBase>& sequence, list<EvolvedBase>::iterator it, size_t distance) const {
    auto advanced = it;
    for (size_t i = 0; i < distance; ++i) {
        ++advanced;
        if (advanced == sequence.end()) {
            break;
        }
    }
    return advanced;
}

void Evolver::duplicate_subseq(list<EvolvedBase>& sequence,
                               list<EvolvedBase>::iterator begin, list<EvolvedBase>::iterator end) const {
    
    for (auto it = begin; it != end; ++it) {
        sequence.insert(begin, *it);
    }
}

template<class Generator>
void Evolver::point_insert(list<EvolvedBase>& sequence, list<EvolvedBase>::iterator pos,
                           size_t size, Generator& gen) const {
    
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
void Evolver::delete_subseq(list<EvolvedBase>& sequence,
                            list<EvolvedBase>::iterator begin, list<EvolvedBase>::iterator end) const {
    
    sequence.erase(begin, end);
}

string dummy_newick(uint64_t num_generations) {
    stringstream strm;
    strm << "(seq1:" << num_generations << ",seq2:" << num_generations << ");";
    return strm.str();
}

void write_fasta(const list<EvolvedBase>& seq, const string& name, ostream& out) {
    out << '>' << name << '\n';
    size_t i = 0;
    for (auto evolved_base : seq) {
        out << evolved_base.base;
        ++i;
        if (i % 80 == 0) {
            out << '\n';
        }
    }
    // add a terminating newline (if we haven't already)
    if (i % 80 != 0) {
        out << '\n';
    }
}

void write_identity(const list<EvolvedBase>& seq, ostream& out) {
    for (auto evolved_base : seq) {
        out << evolved_base.origin << '\n';
    }
}

int main(int argc, char* argv[]) {
    
    // world's worst unit test setup with a few values from scipy.special.zeta
    double x, q, expected, diff;
    x = 2.0; q = 5.0; expected = 0.22132295573711533;
    diff = fabs(zeta(x, q) - expected);
    assert(diff < 0.00001);
    x = 3.0; q = 10.0; expected = 0.005524917485401034;
    diff = fabs(zeta(x, q) - expected);
    assert(diff < 0.00001);
    x = 3.0; q = 1.0; expected = 1.202056903159594;
    diff = fabs(zeta(x, q) - expected);
    assert(diff < 0.00001);
    x = 20.0; q = 1.0; expected = 1.0000009539620338;
    diff = fabs(zeta(x, q) - expected);
    assert(diff < 0.00001);
    
    
    string prefix = "";
    
    int64_t num_generations = default_num_generations;
    
    Evolver evolver;
    
    double exp_large_hor_indel = default_exp_large_hor_indel;
    double hor_indel_tail_heaviness = default_large_hor_indel_tail_heaviness;
    
    string tree_file;
    
    uint64_t seed = -1;
    bool set_seed = false;
    
    static const int help_opt = 1000;
    while (true)
    {
        static struct option options[] = {
            {"output", required_argument, NULL, 'o'},
            {"generations", required_argument, NULL, 'g'},
            {"tree", required_argument, NULL, 'T'},
            {"hor-indel-small-rate", required_argument, NULL, 'h'},
            {"hor-indel-small-size", required_argument, NULL, 'H'},
            {"hor-indel-large-rate", required_argument, NULL, 'r'},
            {"hor-indel-large-size", required_argument, NULL, 'R'},
            {"hor-indel-heaviness", required_argument, NULL, 't'},
            {"monomer-indel-rate", required_argument, NULL, 'm'},
            {"monomer-indel-size", required_argument, NULL, 'M'},
            {"point-indel-rate", required_argument, NULL, 'p'},
            {"point-indel-size", required_argument, NULL, 'P'},
            {"substitution-rate", required_argument, NULL, 's'},
            {"seed", required_argument, NULL, 'z'},
            {"help", no_argument, NULL, help_opt},
            {NULL, 0, NULL, 0}
        };
        int o = getopt_long(argc, argv, "o:T:g:h:H:r:R:t:m:M:p:P:s:z:", options, NULL);
        
        if (o == -1) {
            // end of uptions
            break;
        }
        switch (o)
        {
            case 'o':
                prefix = optarg;
                break;
            case 'g':
                num_generations = parse_int(optarg);
                break;
            case 'T':
                tree_file = optarg;
                break;
            case 'h':
                evolver.small_hor_indel_rate = parse_double(optarg);
                break;
            case 'H':
                evolver.exp_small_hor_indel = parse_double(optarg);
            case 'r':
                evolver.large_hor_indel_rate = parse_double(optarg);
                break;
            case 'R':
                exp_large_hor_indel = parse_double(optarg);
                break;
            case 't':
                hor_indel_tail_heaviness = parse_double(optarg);
                break;
            case 'm':
                evolver.monomer_indel_rate = parse_double(optarg);
                break;
            case 'M':
                evolver.exp_monomer_indel = parse_double(optarg);
                break;
            case 'p':
                evolver.point_indel_rate = parse_double(optarg);
                break;
            case 'P':
                evolver.exp_point_indel = parse_double(optarg);
                break;
            case 's':
                evolver.subs_rate = parse_double(optarg);
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
    
    if (prefix.empty()) {
        cerr << "error: output prefix is required\n";
        return 1;
    }
    
    cerr << "parsing inputs\n";
    
    string fasta = argv[optind++];
    string bed = argv[optind];
    
    if (!set_seed) {
        random_device rd;
        seed = rd();
    }
    cerr << "seed is " << seed << '\n';
    
    assert(hor_indel_tail_heaviness > 0.0);
    evolver.large_hor_indel_beta = 1.0 + 1.0 / hor_indel_tail_heaviness;
    evolver.large_hor_indel_sigma = choose_discrete_pareto_sigma(exp_large_hor_indel, evolver.large_hor_indel_beta);
    
    mt19937 gen(seed);
    
    string newick;
    if (tree_file.empty()) {
        newick = dummy_newick(num_generations);
    }
    else {
        ifstream tree_in;
        stringstream sstrm;
        sstrm << get_input(tree_file, tree_in)->rdbuf();
        newick = sstrm.str();
    }
    Tree tree(newick);
    
    for (uint64_t node_id = 0; node_id < tree.node_size(); ++node_id) {
        if (tree.is_leaf(node_id) && tree.label(node_id).empty()) {
            throw std::runtime_error("leaf node in tree does not have a label");
        }
    }
    
    //    size_t num_hor_indels = 0;
    //    size_t size_hor_indels = 0;
    //    size_t num_mon_indels = 0;
    //    size_t size_mon_indels = 0;
    //    size_t num_point_indels = 0;
    //    size_t size_point_indels = 0;
    //    size_t num_substitutions = 0;
    
    
    vector<list<EvolvedBase>> sequences(tree.node_size());
    
    for (uint64_t node_id : tree.preorder()) {
        if (node_id == tree.get_root()) {
            cerr << "initializing root sequence (id " << node_id << ")\n";
            sequences[node_id] = initialize_sequence(fasta, bed);
        }
        else {
            const auto& parent_seq = sequences[tree.get_parent(node_id)];
            
            double num_gens = tree.distance(node_id);
            assert(num_gens == double(uint64_t(num_gens)));
            
            cerr << "mutating " << num_gens << " generations from id " << tree.get_parent(node_id) << " to id " << node_id;
            if (tree.is_leaf(node_id)) {
                cerr << " (" << tree.label(node_id) << ")";
            }
            cerr << '\n';
            
            sequences[node_id] = evolver.evolve(parent_seq, num_gens, gen);
        }
    }
    
    for (uint64_t node_id = 0; node_id < tree.node_size(); ++node_id) {
        if (tree.is_leaf(node_id)) {
            string fasta_filename = prefix + "_" + tree.label(node_id) + ".fasta";
            ofstream fasta_out(fasta_filename);
            if (!fasta_out) {
                cerr << "error: could not write to " << fasta_filename << '\n';
                return 1;
            }
            string id_filename = prefix + "_" + tree.label(node_id) + "_identity.txt";
            ofstream id_out(id_filename);
            if (!id_out) {
                cerr << "error: could not write to " << id_filename << '\n';
                return 1;
            }
            
            write_fasta(sequences[node_id], tree.label(node_id), fasta_out);
            write_identity(sequences[node_id], id_out);
            
            vector<size_t> origin1;
            for (const auto& base : sequences[node_id]) {
                origin1.push_back(base.origin);
            }
            
            // do pairwise alignments
            for (uint64_t other_id = node_id + 1; other_id < tree.node_size(); ++other_id) {
                if (tree.is_leaf(other_id)) {
                    vector<size_t> origin2;
                    for (const auto& base : sequences[other_id]) {
                        origin2.push_back(base.origin);
                    }
                    auto alignment = move(align_hs(origin1, origin2));
                    // write the CIGAR string
                    string cigar_filename = prefix + tree.label(node_id) + "_" + tree.label(other_id) + "_cigar.txt";
                    ofstream cigar_out(cigar_filename);
                    if (!cigar_out) {
                        cerr << "error: failed to write to " << cigar_filename << '\n';
                    }
                    cigar_out << cigar(alignment) << '\n';
                }
            }
        }
    }
    
//    // compile the summary info
//    stringstream info_strm;
//    info_strm << "summary:\n";
//    info_strm << "\tseed: " << seed << "\n";
//    info_strm << "\tsubstitutions: " << num_substitutions << '\n';
//    info_strm << "\tpoint indels: " << num_point_indels << ", " << size_point_indels << " bases\n";
//    info_strm << "\tmonomer indels: " << num_mon_indels << ", " << size_mon_indels << " monomers\n";
//    info_strm << "\tHOR indels: " << num_hor_indels << ", " << size_hor_indels << " HORs\n";
//    if (find_opt_alignment) {
//        info_strm << "\tmax recoverable aligned bases: " << num_matches << ", prop " << prop_matches << '\n';
//    }
//
//    // write it to stderr and to a file
//    string info = info_strm.str();
//    cerr << info;
//    string info_filename = prefix + "_info.txt";
//    ofstream info_out(info_filename);
//    if (!info_out) {
//        cerr << "error: failed to write to " << info_filename << '\n';
//    }
//    info_out << info;
    
}
