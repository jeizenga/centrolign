#include "centrolign/execution.hpp"

#include <stdexcept>

#include "centrolign/logging.hpp"
#include "centrolign/gfa.hpp"
#include "centrolign/step_index.hpp"

namespace centrolign {

using namespace std;

Execution::Execution(bool suppress_logging) : suppress_logging(suppress_logging) {
    
}

void Execution::init(std::vector<std::pair<std::string, std::string>>&& names_and_sequences,
                     Tree&& tree_in) {
    // TODO: get rid of tree as a member
    
    auto sequences = std::move(names_and_sequences);
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
    
    logging::log(suppress_logging ? logging::Debug : logging::Basic, "Initializing leaf subproblems.");
    
    subproblems.resize(tree.node_size());
    subproblem_finished.resize(tree.node_size());
    for (uint64_t node_id = 0; node_id < tree.node_size(); ++node_id) {
        if (tree.is_leaf(node_id)) {
            const auto& name = tree.label(node_id);
            const auto& sequence = sequences[name_to_idx[name]].second;
            
            auto& subproblem = subproblems[node_id];
            
            subproblem.graph = make_base_graph(name, sequence);
            subproblem.tableau = add_sentinels(subproblem.graph, 5, 6);
            subproblem.name = name;
            subproblem_finished[node_id] = true;
        }
    }
    
    log_memory_usage(logging::Debug);
}


std::vector<Subproblem*> Execution::leaf_subproblems() {
    
    std::vector<Subproblem*> leaves;
    leaves.reserve(tree.node_size() / 2);
    
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

SubproblemScheduler* Execution::get_scheduler() {
    
    if (!scheduler.get()) {
        scheduler.reset(task_parallel ? nullptr : new SerialScheduler(tree, threads));
    }
    return scheduler.get();
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
                        bool preserve_leaves) {
    
    int num_restarted = 0;
    int num_pruned = 0;
    for (auto node_id : tree.preorder()) {
        
        if (subproblem_finished[node_id]) {
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
            subproblem_finished[node_id] = true;
            // FIXME: we dont' save the subproblems alignments, but for now that's not a problem
            //subproblem.alignment = ?
            log_memory_usage(logging::Debug);
            
            // mark descendents as complete
            auto stack = tree.get_children(node_id);
            while (!stack.empty()) {
                auto top = stack.back();
                stack.pop_back();
                subproblem_finished[top] = true;
                if (!tree.is_leaf(top)) {
                    ++num_pruned;
                }
                if (!(preserve_leaves && tree.is_leaf(top)) &&
                    !(preserve_subproblems && !tree.is_leaf(top))) {
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

void Execution::finish_subproblem(const Subproblem& subproblem) {

    uint64_t node_id = get_tree_id(subproblem);
    subproblem_finished[node_id] = true;
    
    if (!preserve_subproblems) {
        // clear out the children that we don't need anymore
        auto& child_subproblem1 = subproblems[tree.get_children(node_id).front()];
        auto& child_subproblem2 = subproblems[tree.get_children(node_id).back()];
        auto dummy_graph1 = std::move(child_subproblem1.graph);
        auto dummy_graph2 = std::move(child_subproblem2.graph);
        auto dummy_aln1 = std::move(child_subproblem1.alignment);
        auto dummy_aln2 = std::move(child_subproblem2.alignment);
    }
}


bool Execution::is_complete(const Subproblem& subproblem) {

    uint64_t node_id = get_tree_id(subproblem);
    return subproblem_finished[node_id];
}

bool Execution::finished() {
    
    return get_scheduler()->finished();
}

ProgressiveStep Execution::next() {
    
    // let the scheduler decide what comes next
    uint64_t node_id, task_threads;
    std::tie(node_id, task_threads) = get_scheduler()->next();
    
    // convert it into a progression MSA step
    ProgressiveStep next_step;
    next_step.parent = &subproblems[node_id];
    
    const auto& children = tree.get_children(node_id);
    
    if (children.size() != 2) {
        std::cerr << "error invalid tree: " << tree.to_newick() << '\n';
        throw std::runtime_error("Attempting execution with a tree that is not binary");
    }
    
    next_step.child1 = &subproblems[children.front()];
    next_step.child2 = &subproblems[children.back()];
    next_step.thread_budget = task_threads;
    
    if ((!suppress_logging && logging::level >= logging::Verbose) || logging::level == logging::Debug) {
        stringstream strm;
        strm << "Next subproblem contains sequences:\n";
        for (auto leaf_name : leaf_descendents(*next_step.parent)) {
            strm << '\t' << leaf_name << '\n';
        }
        logging::log(logging::Verbose, strm.str());
    }
    
    return next_step;
}

MinorExecution::MinorExecution() : Execution(true) {
    
}

MainExecution::MainExecution() : Execution(false) {
    
}

void MainExecution::finish_subproblem(const Subproblem& subproblem) {
    
    if (!subalignments_filepath.empty()) {
        
        uint64_t node_id = get_tree_id(subproblem);
        
        const auto& child1 = subproblems[tree.get_children(node_id).front()];
        const auto& child2 = subproblems[tree.get_children(node_id).back()];
        const auto& graph1 = child1.graph;
        const auto& graph2 = child2.graph;
        
        ofstream out(subalignments_filepath, ios_base::app);
        if (!out) {
            throw std::runtime_error("Failed to write to subalignment file " + subalignments_filepath);
        }
        
        out << "# sequence set 1\n";
        for (const auto& seq_name : leaf_descendents(child1)) {
            out << seq_name << '\n';
        }
        out << "# sequence set 2\n";
        for (const auto& seq_name : leaf_descendents(child2)) {
            out << seq_name << '\n';
        }
        
        StepIndex step_index1(graph1);
        StepIndex step_index2(graph2);
        out << "# alignment\n";
        for (const auto& aln_pair : subproblem.alignment) {
            if (aln_pair.node_id1 == AlignedPair::gap) {
                out << "-\t-\t-";
            }
            else {
                uint64_t path_id;
                size_t step;
                tie(path_id, step) = step_index1.path_steps(aln_pair.node_id1).front();
                out << graph1.path_name(path_id) << '\t' << step << '\t' << decode_base(graph1.label(graph1.path(path_id)[step]));
            }
            out << '\t';
            if (aln_pair.node_id2 == AlignedPair::gap) {
                out << "-\t-\t-";
            }
            else {
                uint64_t path_id;
                size_t step;
                tie(path_id, step) = step_index2.path_steps(aln_pair.node_id2).front();
                out << graph2.path_name(path_id) << '\t' << step << '\t' << decode_base(graph2.label(graph2.path(path_id)[step]));
            }
            out << '\n';
        }
        
    }
    
    Execution::finish_subproblem(subproblem);
}

SubproblemScheduler::SubproblemScheduler(uint64_t threads) : threads(threads) {
    
}

SerialScheduler::SerialScheduler(const Tree& tree, uint64_t threads) : SubproblemScheduler(threads) {
    
    // set up the execution order
    execution_order.reserve(tree.node_size() / 2);
    for (auto tree_id : tree.small_first_postorder()) {
        if (!tree.is_leaf(tree_id)) {
            execution_order.push_back(tree_id);
        }
    }
}

bool SerialScheduler::finished() const {
    return next_subproblem >= execution_order.size();
}

std::pair<uint64_t, uint64_t> SerialScheduler::next() {
    
    return std::pair<uint64_t, uint64_t>(execution_order[next_subproblem++], threads);
}

}
