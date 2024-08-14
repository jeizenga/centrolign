#ifndef centrolign_bonder_hpp
#define centrolign_bonder_hpp

#include <vector>
#include <array>
#include <utility>
#include <stdexcept>
#include <limits>

#include "centrolign/anchorer.hpp"
#include "centrolign/partitioner.hpp"
#include "centrolign/utility.hpp"
#include "centrolign/step_index.hpp"
#include "centrolign/superbubble_distance_oracle.hpp"

namespace centrolign {


/*
 * An interval of two paths that we've identified as mergeable in a cycle motif
 */
struct bond_t {
    
    bond_t() noexcept = default;
    bond_t(const bond_t& other) noexcept = default;
    bond_t(bond_t&& other) noexcept = default;
    ~bond_t() = default;
    bond_t& operator=(const bond_t& other) noexcept = default;
    bond_t& operator=(bond_t&& other) noexcept = default;
    
    std::string path1;
    std::string path2;
    size_t offset1;
    size_t offset2;
    size_t length;
    double score;
};

/*
 * A series of bonds that represents a future merge, creating a cycle
 */
using bond_interval_t = std::vector<bond_t>;

/*
 * Class to identify cycles for the alignment from anchors
 */
class Bonder : public Extractor, public PartitionClient {
public:
    Bonder() = default;
    ~Bonder() = default;
    
    template<class BGraph, class XMerge>
    std::vector<bond_interval_t> identify_bonds(const BGraph& graph1, const BGraph& graph2,
                                                const SentinelTableau& tableau1, const SentinelTableau& tableau2,
                                                const XMerge& xmerge1, const XMerge& xmerge2,
                                                const std::vector<anchor_t>& opt_chain,
                                                const std::vector<anchor_t>& secondary_chain) const;
    
    void deduplicate_self_bonds(std::vector<bond_interval_t>& bonds) const;
    
    enum BondAlgorithm {Null, LongestNearOpt, LongestWindowedNearOpt, LongestNearOptDevConstrained};
    
    BondAlgorithm bond_algorithm = LongestNearOptDevConstrained;
    
    double min_opt_proportion = 0.2;
    
    bool include_gap_scores = true;
    
    double min_length = 100000.0;
    
    double window_length = 75000.0;
    
    double deviation_drift_factor = 150.0;
    
    double separation_drift_factor = 50.0;
    
    // TODO: extend this to the unwindowed algorithm?
    bool break_intervening_windows = true;
    
    double deduplication_slosh_proportion = 0.1;
    
protected:
    
    // passed in as records of (length, opt segment score, secondary segment score)
    std::vector<std::pair<size_t, size_t>> longest_partition(const std::vector<std::tuple<double, double, double>>& shared_subanchors,
                                                             const std::vector<std::tuple<double, double, double>>& intervening_segments) const;
    
    std::vector<std::pair<size_t, size_t>> longest_windowed_partition(const std::vector<std::tuple<double, double, double>>& shared_subanchors,
                                                                      const std::vector<std::tuple<double, double, double>>& intervening_segments) const;
    
    
    std::vector<std::pair<size_t, size_t>>
    longest_deviation_constrained_partition(const std::vector<std::tuple<double, double, double>>& shared_subanchors,
                                            const std::vector<std::tuple<double, double, double>>& intervening_segments,
                                            const std::vector<std::pair<int64_t, int64_t>>& deviation,
                                            const std::vector<std::tuple<uint64_t, uint64_t, uint64_t, uint64_t>>* shared_node_ids = nullptr,
                                            const SuperbubbleDistanceOracle* bond_oracle = nullptr) const;
    
    
    // a counter, solely for instrumentation/development
    static int output_num;
};



/*
 * Template implementations
 */

template<class BGraph, class XMerge>
std::vector<bond_interval_t> Bonder::identify_bonds(const BGraph& graph1, const BGraph& graph2,
                                                    const SentinelTableau& tableau1, const SentinelTableau& tableau2,
                                                    const XMerge& xmerge1, const XMerge& xmerge2,
                                                    const std::vector<anchor_t>& opt_chain,
                                                    const std::vector<anchor_t>& secondary_chain) const {
    
    
    logging::log(logging::Debug, "Identifying bonds on graphs containing " + graph1.path_name(0) + " and " + graph2.path_name(0));
    
    static const bool debug = false;
    
    std::vector<bond_interval_t> bonds;
    
    if (bond_algorithm != Null) {
        for (bool on_graph1 : {true, false}) {
            
            if (debug) {
                std::cerr << "looking for near opt intervals on graph " << (on_graph1 ? 1 : 2) << '\n';
            }
            
            // project overlapping anchors on one graph to find bonds on the other
            
            const auto& proj_graph = on_graph1 ? graph1 : graph2;
            const auto& bond_graph = on_graph1 ? graph2 : graph1;
            auto proj_walk = [&](const anchor_t& anchor) -> const std::vector<uint64_t>& {
                return on_graph1 ? anchor.walk1 : anchor.walk2;
            };
            auto bond_walk = [&](const anchor_t& anchor) -> const std::vector<uint64_t>& {
                return on_graph1 ? anchor.walk2 : anchor.walk1;
            };
            
            std::vector<std::pair<size_t, size_t>> node_locations(proj_graph.node_size(),
                                                                  std::pair<size_t, size_t>(-1, -1));
            
            // record where anchors from the optimal chain occur
            for (size_t i = 0; i < opt_chain.size(); ++i) {
                const auto& anchor = opt_chain[i];
                for (size_t j = 0; j < anchor.walk1.size(); ++j) {
                    node_locations[proj_walk(anchor)[j]] = std::make_pair(i, j);
                }
            }
            
            // records of (sec anchor, idx on sec anchor, opt anchor, idx on opt anchor, length)
            std::vector<std::tuple<size_t, size_t, size_t, size_t, size_t>> shared_subanchors;
            
            for (size_t i = 0; i < secondary_chain.size(); ++i) {
                size_t prev_k = -1, prev_l = -1;
                const auto& anchor = secondary_chain[i];
                for (size_t j = 0; j < anchor.walk1.size(); ++j) {
                    size_t k, l;
                    std::tie(k, l) = node_locations[proj_walk(anchor)[j]];
                    if (k != -1) {
                        if (prev_k == k && prev_l == l - 1) {
                            ++std::get<4>(shared_subanchors.back());
                        }
                        else {
                            shared_subanchors.emplace_back(i, j, k, l, 1);
                        }
                    }
                    prev_k = k;
                    prev_l = l;
                }
            }
            
            {
                // clear this memory
                auto dummy = std::move(node_locations);
            }
            
            if (!shared_subanchors.empty()) {
                if (debug) {
                    std::cerr << "shared subanchors:\n";
                    
                    for (size_t idx = 0; idx < shared_subanchors.size(); ++idx) {
                        size_t i, j, k, l, length;
                        std::tie(i, j, k, l, length) = shared_subanchors[idx];
                        std::cerr << idx << '\t' << i << '\t' << j << '\t' << secondary_chain[i].walk1[j] << '\t' << secondary_chain[i].walk2[j] << '\t' << k << '\t' << l << '\t' << opt_chain[k].walk1[l] << '\t' << opt_chain[k].walk2[l] << '\t' << length << '\n';
                    }
                }
                
                // measure the distance between anchors on the projecting graph
                std::vector<double> dist_between(opt_chain.size() - 1);
                {
                    auto subgraphs_between = extract_graphs_between(opt_chain, graph1, graph2, tableau1, tableau2,
                                                                    xmerge1, xmerge2);
                    for (size_t i = 1; i + 1 < subgraphs_between.size(); ++i) {
                        const auto& measurement_subgraph = on_graph1 ? subgraphs_between[i].first : subgraphs_between[i].second;
                        if (measurement_subgraph.subgraph.node_size() != 0) {
                            dist_between[i - 1] = source_sink_minmax(measurement_subgraph).first;
                        }
                        else {
                            dist_between[i - 1] = 0.0;
                        }
                    }
                }
                
                // tuples of (length, score opt, score secondary)
                std::vector<std::tuple<double, double, double>> shared_segments(shared_subanchors.size());
                std::vector<std::tuple<double, double, double>> intervening_segments(shared_subanchors.size() - 1);
                // pairs of (opt deviation, secondary deviation)
                std::vector<std::pair<int64_t, int64_t>> deviation;
                // for shared segments, tuples of (opt bond start id, opt bond final id, sec bond start id, sec bond final id)
                std::vector<std::tuple<uint64_t, uint64_t, uint64_t, uint64_t>> shared_node_ids;
                if (bond_algorithm == LongestNearOptDevConstrained) {
                    deviation.resize(intervening_segments.size(), std::pair<int64_t, int64_t>(0, 0));
                    
                    shared_node_ids.resize(shared_subanchors.size());
                    for (size_t i = 0; i < shared_node_ids.size(); ++i) {
                        const auto& shared = shared_subanchors[i];
                        shared_node_ids[i] = std::make_tuple(bond_walk(opt_chain[std::get<2>(shared)])[std::get<3>(shared)],
                                                             bond_walk(opt_chain[std::get<2>(shared)])[std::get<3>(shared) + std::get<4>(shared) - 1],
                                                             bond_walk(secondary_chain[std::get<0>(shared)])[std::get<1>(shared)],
                                                             bond_walk(secondary_chain[std::get<0>(shared)])[std::get<1>(shared) + std::get<4>(shared) - 1]);
                    }
                }
                
                for (size_t idx = 0; idx < shared_subanchors.size(); ++idx) {
                    
                    // shared segment is part of an anchor
                    
                    size_t i, j, k, l, length;
                    std::tie(i, j, k, l, length) = shared_subanchors[idx];
                    auto& segment = shared_segments[idx];
                    std::get<0>(segment) = length;
                    std::get<1>(segment) = (length * opt_chain[k].score) / opt_chain[k].walk1.size();
                    std::get<2>(segment) = (length * secondary_chain[i].score) / secondary_chain[i].walk1.size();
                    
                    if (idx != 0) {
                        // in-between segment is all partial anchors and gaps between anchors until next shared segment
                        
                        auto& between = intervening_segments[idx - 1];
                        
                        size_t pi, pj, pk, pl, plength;
                        std::tie(pi, pj, pk, pl, plength) = shared_subanchors[idx - 1];
                        
                        // add up the length and score on the optimal chain
                        if (pk == k) {
                            std::get<0>(between) = l - pl - plength;
                            std::get<1>(between) = (std::get<0>(between) * opt_chain[k].score) / opt_chain[k].walk1.size();
                        }
                        else {
                            for (size_t x = pk, offset = pl + plength; x <= k; ++x) {
                                size_t sublen = (x == k ? l : opt_chain[x].walk1.size() - offset);
                                std::get<0>(between) += sublen;
                                std::get<1>(between) += (sublen * opt_chain[x].score) / opt_chain[x].walk1.size();
                                if (x != k) {
                                    std::get<0>(between) += dist_between[x];
                                    if (include_gap_scores) {
                                        std::get<1>(between) += opt_chain[x].gap_score_after;
                                    }
                                    if (!deviation.empty()) {
                                        deviation[idx - 1].first += opt_chain[x].gap_after;
                                    }
                                }
                                offset = 0;
                            }
                        }
                        
                        // add up the length and score on the secondary chain
                        if (pi == i) {
                            std::get<2>(between) = ((j - pj - plength) * secondary_chain[i].score) / secondary_chain[i].walk1.size();
                        }
                        else {
                            for (size_t x = pi, offset = pj + plength; x <= i; ++x) {
                                size_t sublen = (x == i ? j : secondary_chain[x].walk1.size() - offset);
                                std::get<2>(between) += (sublen * secondary_chain[x].score) / secondary_chain[x].walk1.size();
                                if (x != i) {
                                    if (include_gap_scores) {
                                        std::get<2>(between) += secondary_chain[x].gap_score_after;
                                    }
                                    if (!deviation.empty()) {
                                        deviation[idx - 1].second += secondary_chain[x].gap_after;
                                    }
                                }
                                offset = 0;
                            }
                        }
                    }
                }
                
                
                static const bool instrument_segments = false;
                if (instrument_segments) {
                    for (size_t idx = 0; idx < shared_segments.size(); ++idx) {
                        if (idx != 0) {
                            size_t pi, pj, pk, pl, plen;
                            std::tie(pi, pj, pk, pl, plen) = shared_subanchors[idx - 1];
                            const auto& prev_opt_anchor = opt_chain[pk];
                            const auto& prev_sec_anchor = secondary_chain[pi];
                            double ls, o, s;
                            std::tie(ls, o, s) = intervening_segments[idx - 1];
                            std::cerr << '~' << '\t' << (2 * idx - 1) << '\t' << on_graph1 << '\t' << 'i' << '\t' << prev_opt_anchor.walk1[pl + plen - 1] << '\t' << prev_opt_anchor.walk2[pl + plen - 1] << '\t' << prev_sec_anchor.walk1[pj + plen - 1] << '\t' << prev_sec_anchor.walk2[pj + plen - 1] << '\t' << ls << '\t' << o << '\t' << s;
                            if (!deviation.empty()) {
                                std::cerr << '\t' << deviation[idx - 1].first << '\t' << deviation[idx - 1].second;
                            }
                            std::cerr << '\n';
                        }
                        size_t i, j, k, l, len;
                        std::tie(i, j, k, l, len) = shared_subanchors[idx];
                        const auto& opt_anchor = opt_chain[k];
                        const auto& sec_anchor = secondary_chain[i];
                        double ls, o, s;
                        std::tie(ls, o, s) = shared_segments[idx];
                        std::cerr << '~' << '\t' << (2 * idx) << '\t' << on_graph1 << '\t' << 's' << '\t' << opt_anchor.walk1[l] << '\t' << opt_anchor.walk2[l] << '\t' << sec_anchor.walk1[j] << '\t' << sec_anchor.walk2[j] << '\t' << ls << '\t' << o << '\t' << s;
                        if (!deviation.empty()) {
                            std::cerr << '\t' << 0 << '\t' << 0;
                        }
                        std::cerr << '\n';
                    }
                }
                
                // find the best partition into bonds
                std::vector<std::pair<size_t, size_t>> partition;
                if (bond_algorithm == LongestNearOpt) {
                    partition = std::move(longest_partition(shared_segments, intervening_segments));
                }
                else if (bond_algorithm == LongestWindowedNearOpt) {
                    partition = std::move(longest_windowed_partition(shared_segments, intervening_segments));
                }
                else if (bond_algorithm == LongestNearOptDevConstrained) {
                    
                    SuperbubbleDistanceOracle bond_oracle(proj_graph);
                    
                    partition = std::move(longest_deviation_constrained_partition(shared_segments, intervening_segments, deviation, &shared_node_ids, &bond_oracle));
                }
                else {
                    throw std::runtime_error("Unrecognized bond algorithm: " + std::to_string((int) bond_algorithm));
                }
                
                if (debug) {
                    std::cerr << "found " << partition.size() << " near opt intervals among " << shared_segments.size() << " shared segments\n";
                }
                
                
                if (!partition.empty()) {
                    
                    // we found bonds
                    
                    StepIndex step_index(bond_graph);
                    
                    for (const auto& interval : partition) {
                        
                        if (debug) {
                            std::cerr << "processing interval " << interval.first << ", " << interval.second << '\n';
                            double length = 0.0, opt_sc = 0.0, sec_sc = 0.0;
                            for (size_t i = interval.first; i < interval.second; ++i) {
                                const auto& shared = shared_segments[i];
                                length += std::get<0>(shared);
                                opt_sc += std::get<1>(shared);
                                sec_sc += std::get<2>(shared);
                                if (i != interval.first) {
                                    const auto& between = intervening_segments[i - 1];
                                    length += std::get<0>(between);
                                    opt_sc += std::get<1>(between);
                                    sec_sc += std::get<2>(between);
                                }
                            }
                            std::cerr << "length " << length << ", opt score " << opt_sc << ", secondary score " << sec_sc << '\n';
                            if (!deviation.empty()) {
                                int64_t opt_dev = 0, sec_dev = 0;
                                int64_t max_diff = 0, min_diff = 0;
                                for (size_t i = interval.first + 1; i < interval.second; ++i) {
                                    opt_dev += deviation[i - 1].first;
                                    sec_dev += deviation[i - 1].second;
                                    max_diff = std::max(max_diff, opt_dev - sec_dev);
                                    min_diff = std::min(min_diff, opt_dev - sec_dev);
                                }
                                std::cerr << "opt deviation " << opt_dev << ", secondary deviation " << sec_dev << ", min diff " << min_diff << ", max diff " << max_diff << '\n';
                            }
                            std::cerr << "gr" << '\t' << "opt" << '\t' << "po" << '\t' << "sec" << '\t' << "ps" << '\t' << "len" << '\t' << "oid" << '\t' << "sid" << '\n';
                        }
                        
                        // each interval of shared segments corresponds to a bond interval we need to emit
                        
                        bonds.emplace_back();
                        auto& bond_interval = bonds.back();
                        
                        for (size_t idx = interval.first; idx < interval.second; ++idx) {
                            // add a shared sub anchor
                            
                            size_t i, j, k, l, len;
                            std::tie(i, j, k, l, len) = shared_subanchors[idx];
                            
                            if (debug) {
                                std::cerr << on_graph1 << '\t' << k << '\t' << l << '\t' << i << '\t' << j << '\t' << len << '\t' << bond_walk(opt_chain[k])[l] << '\t' << bond_walk(secondary_chain[i])[j] << '\n';
                            }
                            
                            const auto& walk_opt = bond_walk(opt_chain[k]);
                            const auto& walk_sec = bond_walk(secondary_chain[i]);
                            
                            uint64_t curr_path_id1 = -1, curr_path_id2 = -1;
                            for (size_t x = 0; x < len; ++x) {
                                
                                uint64_t path_id1, path_id2;
                                size_t offset1, offset2;
                                
                                std::tie(path_id1, offset1) = step_index.path_steps(walk_opt[l + x]).front();
                                std::tie(path_id2, offset2) = step_index.path_steps(walk_sec[j + x]).front();
                                
                                if (bond_interval.empty() || path_id1 != curr_path_id1 || path_id2 != curr_path_id2
                                    || bond_interval.back().offset1 + bond_interval.back().length != offset1
                                    || bond_interval.back().offset2 + bond_interval.back().length != offset2) {
                                    
                                    // we need to start a new bond
                                    
                                    if (!bond_interval.empty()) {
                                        // annotate the previous one's score
                                        bond_interval.back().score = (bond_interval.back().length * secondary_chain[i].score) / walk_sec.size();
                                    }
                                    
                                    bond_interval.emplace_back();
                                    auto& bond = bond_interval.back();
                                    bond.path1 = bond_graph.path_name(path_id1);
                                    bond.path2 = bond_graph.path_name(path_id2);
                                    bond.offset1 = offset1;
                                    bond.offset2 = offset2;
                                    bond.length = 1;
                                }
                                else {
                                    // we can extend the previous bond
                                    ++bond_interval.back().length;
                                }
                                
                                if (!bond_interval.empty()) {
                                    // annotate the final bond's score
                                    bond_interval.back().score = (bond_interval.back().length * secondary_chain[i].score) / walk_sec.size();
                                }
                                
                                curr_path_id1 = path_id1;
                                curr_path_id2 = path_id2;
                            }
                        }
                    }
                }
            }
        }
    }
    
    static const bool instrument_bonds = true;
    if (instrument_bonds) {
        static const bool short_format = true;
        std::cerr << "instrumenting bonds (total " << bonds.size() << " discovered)\n";
        for (size_t i = 0; i < bonds.size(); ++i) {
            std::cerr << "bond interval " << i << '\n';
            const auto& bond_interval = bonds[i];
            if (short_format) {
                std::cerr << '<' << '\t' << bond_interval.front().path1  << '\t' << bond_interval.front().offset1 << '\t' << bond_interval.back().offset1 << '\t' << bond_interval.front().offset2 << '\t' << bond_interval.back().offset2 << '\n';
            }
            else {
                for (size_t j = 0; j < bond_interval.size(); ++j) {
                    const auto& bond = bond_interval[j];
                    std::cerr << '<' << '\t' << j << '\t' << bond.path1 << '\t' << bond.path2 << '\t' << bond.offset1 << '\t' << bond.offset2 << '\t' << bond.length << '\t' << bond.score << '\n';
                }
            }
        }
    }
    static const bool output_bonds = false;
    if (output_bonds) {
        std::unordered_map<uint64_t, uint64_t> paired_node1, paired_node2;
        for (const auto& anchor : opt_chain) {
            for (size_t i = 0; i < anchor.walk1.size(); ++i) {
                paired_node1[anchor.walk1[i]] = anchor.walk2[i];
                paired_node2[anchor.walk2[i]] = anchor.walk1[i];
            }
        }
        for (size_t i = 0; i < bonds.size(); ++i) {
            const auto& bond_interval = bonds[i];
            uint64_t first1 = -1, first2 = -1;
            for (bool do_opt : {true, false}) {
                for (size_t j = 0; j < bond_interval.size(); ++j) {
                    const auto& bond = bond_interval[j];
                    bool on1 = graph1.has_path(bond.path1);
                    const auto& graph = on1 ? graph1 : graph2;
                    const auto& paired_node = on1 ? paired_node1 : paired_node2;
                    const auto& path_name = do_opt ? bond.path1 : bond.path2;
                    const auto& offset = do_opt ? bond.offset1 : bond.offset2;
                    auto opt_node = graph.path(graph.path_id(bond.path1))[bond.offset1];
                    auto paired_id = paired_node.at(opt_node);
                    auto here_id = graph.path(graph.path_id(path_name))[offset];
                    if (!on1) {
                        std::swap(paired_id, here_id);
                    }
                    if (first1 == -1) {
                        first1 = here_id;
                        first2 = paired_id;
                    }
                    std::cout << here_id << '\t' << paired_id << '\t' << bond.length << '\t' << -(output_num + 1) << '\n';
                }
            }
            // close the bond to make an hourglass shape
            std::cout << first1 << '\t' << first2 << '\t' << 0 << '\t' << -(output_num + 1) << '\n';
            ++output_num;
        }
    }
    
    return bonds;
}



}

#endif /* centrolign_bond_identifier_hpp */
