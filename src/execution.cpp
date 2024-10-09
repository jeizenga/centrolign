#include "centrolign/execution.hpp"

#include <stdexcept>

#include "centrolign/logging.hpp"
#include "centrolign/gfa.hpp"

namespace centrolign {

using namespace std;

Execution::Execution(std::vector<std::pair<std::string, std::string>>&& names_and_sequences,
                     Tree&& tree_in) {
    
    // TODO: get rid of tree as a member
    
    auto sequences = move(names_and_sequences);
    tree = std::move(tree_in);
    
    unordered_map<string, size_t> name_to_idx;
    for (size_t i = 0; i < sequences.size(); ++i) {
        auto& name = sequences[i].first;
        if (name_to_idx.count(name)) {
            throw runtime_error("FASTA contains duplicate name " + name);
        }
        name_to_idx[name] = i;
    }
    
    // check the match between the fasta and the tree
    {
        std::vector<uint64_t> sequence_leaf_ids;
        for (size_t i = 0; i < sequences.size(); ++i) {
            
            const string& name = sequences[i].first;
            if (!tree.has_label(name)) {
                throw runtime_error("Guide tree does not include sequence " + name);
            }
            auto node_id = tree.get_id(name);
            if (!tree.is_leaf(node_id)) {
                throw runtime_error("Sequence " + name + " is not a leaf in the guide tree");
            }
            sequence_leaf_ids.push_back(node_id);
        }
        
        // get rid of samples we don't have the sequence for
        tree.prune(sequence_leaf_ids);
    }
    
    // get rid of non-branching paths
    tree.compact();
    
    if (logging::level >= logging::Debug) {
        Tree polytomized = tree;
        polytomized.polytomize();
        logging::log(logging::Debug, "Simplified and subsetted tree:\n" + polytomized.to_newick());
    }
    
    // convert into a binary tree
    tree.binarize();
    
    logging::log(logging::Debug, "Fully processed tree:\n" + tree.to_newick());
    
    log_memory_usage(logging::Debug);
    
    logging::log(logging::Basic, "Initializing leaf subproblems.");
    
    subproblems.resize(tree.node_size());
    for (uint64_t node_id = 0; node_id < tree.node_size(); ++node_id) {
        if (tree.is_leaf(node_id)) {
            const auto& name = tree.label(node_id);
            const auto& sequence = sequences[name_to_idx[name]].second;
            
            auto& subproblem = subproblems[node_id];
            
            subproblem.graph = make_base_graph(name, sequence);
            subproblem.tableau = add_sentinels(subproblem.graph, 5, 6);
            subproblem.name = name;
            subproblem.complete = true;
        }
    }
    
    // set up the execution order
    for (auto tree_id : tree.postorder()) {
        if (!tree.is_leaf(tree_id)) {
            execution_order.push_back(tree_id);
        }
    }
    next_subproblem = 0;
    
    log_memory_usage(logging::Debug);
}

bool Execution::finished() const {
    return next_subproblem > execution_order.size();
}

std::tuple<Subproblem*, Subproblem*, Subproblem*> Execution::next() {
    
    uint64_t node_id = execution_order[next_subproblem++];
    
    auto& next_problem = subproblems[node_id];
    const auto& children = tree.get_children(node_id);
    
    assert(children.size() == 2);
    
    auto& subproblem1 = subproblems[children.front()];
    auto& subproblem2 = subproblems[children.back()];
    
    return std::make_tuple(&next_problem, &subproblem1, &subproblem2);
}

std::tuple<const Subproblem*, const Subproblem*, const Subproblem*> Execution::current() const {
    
    std::tuple<const Subproblem*, const Subproblem*, const Subproblem*> problems(nullptr, nullptr, nullptr);
    
    if (next_subproblem != 0 || next_subproblem <= execution_order.size()) {
        uint64_t node_id = execution_order[next_subproblem - 1];
        std::get<0>(problems) = &subproblems[node_id];
        const auto& children = tree.get_children(node_id);
        std::get<1>(problems) = &subproblems[children.front()];
        std::get<2>(problems) = &subproblems[children.back()];
    }
    
    return problems;
}


std::vector<Subproblem*> Execution::leaf_subproblems() {
    
    std::vector<Subproblem*> leaves;
    leaves.reserve(execution_order.size() + 1);
    
    for (uint64_t tree_id = 0; tree_id < tree.node_size(); ++tree_id) {
        if (tree.is_leaf(tree_id)) {
            leaves.push_back(&subproblems[tree_id]);
        }
    }
    
    return leaves;
}

uint64_t Execution::get_tree_id(const Subproblem& subproblem) const {
    // compute the index in the subproblem vector
    // TODO: very hacky
    return (&subproblem - subproblems.data());
}


std::vector<std::string> Execution::leaf_descendents(const Subproblem& subproblem) const {
    
    auto tree_id = get_tree_id(subproblem);
    
    std::vector<std::string> descendents;
    if (tree.is_leaf(tree_id)) {
        descendents.push_back(tree.label(tree_id));
    }
    else {
        std::vector<uint64_t> stack = tree.get_children(tree_id);
        while (!stack.empty()) {
            auto here = stack.back();
            stack.pop_back();
            if (tree.is_leaf(here)) {
                descendents.push_back(tree.label(here));
            }
            else {
                for (auto next : tree.get_children(here)) {
                    stack.push_back(next);
                }
            }
        }
    }
    return descendents;
}


uint64_t Execution::subproblem_hash(const Subproblem& subproblem) const {
    auto seq_names = leaf_descendents(subproblem);
    sort(seq_names.begin(), seq_names.end());
    
    // create a hash digest of the sample names
    size_t hsh = 660422875706093811ull;
    for (auto& seq_name : seq_names) {
        hash_combine(hsh, size_t(2110260111091729000ull)); // spacer
        for (char c : seq_name) {
            hash_combine(hsh, c);
        }
    }
    return hsh;
}

Subproblem& Execution::final_subproblem() {
    return subproblems[tree.get_root()];
}

const Subproblem& Execution::final_subproblem() const {
    return subproblems[tree.get_root()];
}

const Subproblem& Execution::leaf_subproblem(const std::string& name) const {
    return subproblems[tree.get_id(name)];
}


const Tree& Execution::get_tree() const {
    return tree;
}


void Execution::restart(std::function<std::string(const Subproblem&)>& file_location,
                        bool preserve_leaves, bool preserve_internal_nodes) {
    
    int num_restarted = 0;
    int num_pruned = 0;
    for (auto node_id : tree.preorder()) {
        
        if (subproblems[node_id].complete) {
            if (!tree.is_leaf(node_id)) {
                ++num_pruned;
            }
            continue;
        }
        
        auto file_name = file_location(subproblems[node_id]);
        
        ifstream gfa_in(file_name);
        
        if (gfa_in) {
            // we have the results of this alignment problem saved
            
            ++num_restarted;
            logging::log(logging::Debug, "Loading previously completed subproblem " + file_name);
            
            auto& subproblem = subproblems[node_id];
            
            subproblem.graph = read_gfa(gfa_in);
            subproblem.tableau = add_sentinels(subproblem.graph, 5, 6);
            subproblem.complete = true;
            // FIXME: we dont' save the subproblems alignments, but for now that's not a problem
            //subproblem.alignment = ?
            log_memory_usage(logging::Debug);
            
            // mark descendents as complete
            auto stack = tree.get_children(node_id);
            while (!stack.empty()) {
                auto top = stack.back();
                stack.pop_back();
                subproblems[top].complete = true;
                if (!tree.is_leaf(top)) {
                    ++num_pruned;
                }
                if (!(preserve_leaves && tree.is_leaf(top)) &&
                    !(preserve_internal_nodes && !tree.is_leaf(top))) {
                    // clear out descendents
                    BaseGraph dummy = std::move(subproblems[top].graph);
                }
                for (auto child_id : tree.get_children(top)) {
                    stack.push_back(child_id);
                }
            }
            
        }
    }
    
    logging::log(logging::Basic, "Loaded results for " + std::to_string(num_restarted) + " subproblem(s) from previously completed run and pruned " + std::to_string(num_pruned) + " of their children as unnecessary.");
}

size_t Execution::memory_size() const {
    size_t size = 0;
    for (const auto& subproblem : subproblems) {
        size += subproblem.graph.memory_size();
        size += subproblem.alignment.capacity() * sizeof(decltype(subproblem.alignment)::value_type);
    }
    return size;
}

}
