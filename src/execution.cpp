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
    
    // collect leaf nodes
    std::vector<uint64_t> leaves;
    leaves.reserve(tree.node_size() / 2 + 1);
    subproblems.resize(tree.node_size());
    for (uint64_t node_id = 0; node_id < tree.node_size(); ++node_id) {
        if (tree.is_leaf(node_id)) {
            leaves.emplace_back(node_id);
        }
    }
    
    std::atomic<size_t> next_idx(0);
    
    // parallelizable init of subproblems in a loop
    auto init_subproblems = [&]() {
        while (true) {
            size_t idx = next_idx++;
            if (idx >= leaves.size()) {
                break;
            }
            uint64_t node_id = leaves[idx];
            
            const auto& name = tree.label(node_id);
            const auto& sequence = sequences[name_to_idx[name]].second;
            
            auto& subproblem = subproblems[node_id];
            
            subproblem.graph = make_base_graph(name, sequence);
            subproblem.tableau = add_sentinels(subproblem.graph, 5, 6);
            subproblem.name = name;
        }
    };
    
    // dispatch worker threads
    std::vector<std::thread> workers;
    for (size_t t = 0; t + 1 < threads; ++t) {
        workers.emplace_back(init_subproblems);
    }
    // do work in main thread
    init_subproblems();
    
    // barrier sync
    for (auto& worker : workers) {
        worker.join();
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
        // the scheduler has not been initialized

        // measure sequence size to calibrate a safe thread wait in the parallel scheduler
        size_t max_size = 0;
        if (task_parallel) {
            for (const auto& subproblem : subproblems) {
                max_size = std::max(max_size, subproblem.graph.node_size());
            }
        }
        
        // init the scheduler once according to the configuration
        if (task_parallel) {
            scheduler.reset(new ParallelScheduler(tree, threads, max_size));
        }
        else {
            scheduler.reset(new SerialScheduler(tree, threads));
        }
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

size_t Execution::memory_size() const {
    size_t size = 0;
    for (const auto& subproblem : subproblems) {
        size += subproblem.graph.memory_size();
        size += subproblem.alignment.capacity() * sizeof(decltype(subproblem.alignment)::value_type);
    }
    return size;
}

void Execution::execute(const std::function<void(const ProgressiveStep&)>& do_subproblem) {
    
    while (!get_scheduler()->finished()) {
        
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
                
        get_scheduler()->handle_task([next_step, node_id, do_subproblem, this]() {
            
            do_subproblem(next_step);
            
            if (!this->preserve_subproblems) {
                // clear out the children that we don't need anymore
                auto dummy_graph1 = std::move(next_step.child1->graph);
                auto dummy_graph2 = std::move(next_step.child2->graph);
                auto dummy_aln1 = std::move(next_step.child1->alignment);
                auto dummy_aln2 = std::move(next_step.child2->alignment);
            }
            
            return node_id;
        });
    }
}

MinorExecution::MinorExecution() : Execution(true) {
    
}

MainExecution::MainExecution() : Execution(false) {
    
}

void MainExecution::restart(bool preserve_leaves) {
    
    
    int num_restarted = 0;
    int num_pruned = 0;
    for (auto node_id : tree.preorder()) {
        
        if (get_scheduler()->is_complete(node_id)) {
            if (!tree.is_leaf(node_id)) {
                ++num_pruned;
            }
            continue;
        }
        
        auto file_name = subproblem_file_name(subproblems[node_id]);
        
        ifstream gfa_in(file_name);
        
        if (gfa_in) {
            // we have the results of this alignment problem saved
            
            ++num_restarted;
            logging::log(logging::Debug, "Loading previously completed subproblem " + file_name);
            
            auto& subproblem = subproblems[node_id];
            
            subproblem.graph = read_gfa(gfa_in);
            subproblem.tableau = add_sentinels(subproblem.graph, 5, 6);
            // FIXME: we dont' save the subproblems alignments, but for now that's not a problem
            //subproblem.alignment = ?
            log_memory_usage(logging::Debug);
            
            // mark descendents as complete in post-order
            std::vector<std::pair<uint64_t, size_t>> stack;
            stack.emplace_back(node_id, 0);
            while (!stack.empty()) {
                auto& top = stack.back();
                if (top.second == tree.get_children(top.first).size()) {
                    get_scheduler()->mark_complete(top.first);
                    if (!tree.is_leaf(top.first)) {
                        ++num_pruned;
                    }
                    if (top.first != node_id &&
                        !(preserve_leaves && tree.is_leaf(top.first)) &&
                        !(preserve_subproblems && !tree.is_leaf(top.first))) {
                        // clear out descendents
                        BaseGraph dummy = std::move(subproblems[top.first].graph);
                    }
                    stack.pop_back();
                }
                else {
                    stack.emplace_back(tree.get_children(top.first)[top.second++], 0);
                }
            }
        }
    }
    
    logging::log(logging::Basic, "Loaded results for " + std::to_string(num_restarted) + " subproblem(s) from previously completed run and pruned " + std::to_string(num_pruned) + " of their children as unnecessary.");
}


uint64_t MainExecution::subproblem_hash(const Subproblem& subproblem) const {
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

std::string MainExecution::subproblem_file_name(const Subproblem& subproblem) const {
    return subproblems_prefix + "_" + to_hex(subproblem_hash(subproblem)) + ".gfa";
}

std::string MainExecution::subproblem_info_file_name() const {
    return subproblems_prefix + "_info.txt";
}

void MainExecution::emit_subproblem(const Subproblem& subproblem) {
    
    auto gfa_file_name = subproblem_file_name(subproblem);
    auto info_file_name = subproblem_info_file_name();
    
    ofstream gfa_out(gfa_file_name);
    if (!gfa_out) {
        throw std::runtime_error("Failed to write to subproblem file " + gfa_file_name);
    }
    write_gfa(subproblem.graph, subproblem.tableau, gfa_out);
    
    auto sequences = leaf_descendents(subproblem);
    sort(sequences.begin(), sequences.end());
    {
        std::lock_guard<std::mutex> lock(subproblem_info_mutex);
        // check if the file already exists
        bool write_header = !(ifstream(info_file_name));
        ofstream info_out(info_file_name, ios_base::app);
        if (!info_out) {
            throw std::runtime_error("Failed to write to subproblem info file " + info_file_name);
        }
        if (write_header) {
            info_out << "filename\tsequences\n";
        }
        info_out << gfa_file_name << '\t' << join(sequences, ",") << '\n';
    }
    
}


void MainExecution::execute(const std::function<void(const ProgressiveStep&)>& do_subproblem) {

    auto do_subproblem_postscript = [do_subproblem, this](const ProgressiveStep& next_step){
        
        do_subproblem(next_step);

        if (!this->subproblems_prefix.empty()) {
            this->emit_subproblem(*next_step.parent);
        }
        
        if (!this->subalignments_filepath.empty()) {
            
            const auto& graph1 = next_step.child1->graph;
            const auto& graph2 = next_step.child2->graph;
            
            StepIndex step_index1(graph1);
            StepIndex step_index2(graph2);
            
            std::lock_guard<std::mutex> lock(this->subalignments_mutex);
            
            ofstream out(subalignments_filepath, ios_base::app);
            if (!out) {
                throw std::runtime_error("Failed to write to subalignment file " + subalignments_filepath);
            }
            
            out << "# sequence set 1\n";
            for (const auto& seq_name : leaf_descendents(*next_step.child1)) {
                out << seq_name << '\n';
            }
            out << "# sequence set 2\n";
            for (const auto& seq_name : leaf_descendents(*next_step.child2)) {
                out << seq_name << '\n';
            }
            out << "# alignment\n";
            for (const auto& aln_pair : next_step.parent->alignment) {
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
    };
    
    Execution::execute(do_subproblem_postscript);
}

SubproblemScheduler::SubproblemScheduler(const Tree& tree, uint64_t threads) : threads(threads), subproblem_finished(tree.node_size(), false) {
    
}

bool SubproblemScheduler::is_complete(uint64_t node_id) const {
    return subproblem_finished[node_id];
}

SerialScheduler::SerialScheduler(const Tree& tree, uint64_t threads) : SubproblemScheduler(tree, threads) {
    
    // set up the execution order
    execution_order.reserve(tree.node_size());
    for (auto tree_id : tree.small_first_postorder()) {
        execution_order.push_back(tree_id);
        if (tree.is_leaf(tree_id)) {
            mark_complete(tree_id);
        }
    }
}

bool SerialScheduler::finished() {
    return next_subproblem >= execution_order.size();
}

std::pair<uint64_t, uint64_t> SerialScheduler::next() {
    return std::pair<uint64_t, uint64_t>(execution_order[next_subproblem], threads);
}

void SerialScheduler::handle_task(const std::function<uint64_t(void)>& task) {
    mark_complete(task());
}

void SerialScheduler::mark_complete(uint64_t node_id) {
    subproblem_finished[node_id] = true;
    while (next_subproblem < execution_order.size() &&
           subproblem_finished[execution_order[next_subproblem]]) {
        ++next_subproblem;
    }
}

ParallelScheduler::ParallelScheduler(const Tree& tree, uint64_t threads, size_t problem_size) : SubproblemScheduler(tree, threads), sleep_ms(std::max<uint64_t>(1, problem_size * 0.0002)), threads_available(threads), tasks_executing(0), memory_executing(0) {
    
    scheduling_info.resize(tree.node_size(), SchedulingInfo());
    
    // count the number of sequences in each subtree (and gather some info while we're at it)
    std::vector<uint64_t> subtree_leaves(tree.node_size(), 0);
    for (uint64_t node_id : tree.postorder()) {
        
        scheduling_info[node_id].parent = tree.get_parent(node_id);
        scheduling_info[node_id].children_remaining = tree.get_children(node_id).size();
        scheduling_info[node_id].queue_iter = queue.end();
        
        if (tree.is_leaf(node_id)) {
            subtree_leaves[node_id] = 1;
        }
        else {
            for (uint64_t child_id : tree.get_children(node_id)) {
                subtree_leaves[node_id] += subtree_leaves[child_id];
            }
        }
    }
    
    // approximate memory (in arbitrary units) for each subproblem
    // for large problems, memory is dominated by the k_1 * k_2 sparse DP query structures
    // there is also a factor of N log N in matches, but we hope this is more constant
    // TODO: determine if (k_1 * G_1 + k_2 * G_2) ever dominates from merge and post-switch structs
    for (uint64_t node_id = 0; node_id < tree.node_size(); ++node_id) {
        if (tree.is_leaf(node_id)) {
            continue;
        }
        else {
            const auto& children = tree.get_children(node_id);
            if (children.size() != 2) {
                throw std::runtime_error("Attempted to determine subproblem schedule for a non-binary tree");
            }
            auto& node_info = scheduling_info[node_id];
            node_info.memory = subtree_leaves[children.front()] * subtree_leaves[children.back()];
            node_info.max_threads = std::max(subtree_leaves[children.front()], subtree_leaves[children.back()]);
            memory_limit = std::max(memory_limit, node_info.memory);
        }
    }
    memory_limit *= max_relative_memory;
    
    // mark all leaves as completed subproblems
    for (uint64_t node_id = 0; node_id < tree.node_size(); ++node_id) {
        if (tree.is_leaf(node_id)) {
//            std::cerr << "marking complete leaf node " << node_id << '\n';
            mark_complete(node_id);
        }
    }
}

bool ParallelScheduler::finished() {
    
    auto check_queue_empty = [&]() -> bool {
        std::lock_guard<std::mutex> lock(queue_mutex);
        return queue.empty();
    };
//    std::cerr << ("enter finished check with queue size " + std::to_string(queue.size()) + ", tasks executing " + std::to_string(tasks_executing.load()) + "\n");
    
    // note: tasks are added to the executing count before leaving the queue mutex and cleared from
    // the executing count after adding follow-on tasks to the queue, so there is never a time that
    // the queue is accessible and empty while there are still tasks remaining in the execution
    while (check_queue_empty() && tasks_executing.load() != 0) {
        // there are no remaining tasks in the queue, but there might be more added if we wait
        std::this_thread::sleep_for(std::chrono::milliseconds(sleep_ms));
    }
    
    return check_queue_empty();
}

std::pair<uint64_t, uint64_t> ParallelScheduler::next() {
    
    // wait until there is a thread available to execute
    while (threads_available.load() == 0) {
        std::this_thread::sleep_for(std::chrono::milliseconds(sleep_ms));
    }
    
    // get an upper bound on the number of executing tasks in the future
    uint64_t frontier_size = 0;
    {
        std::lock_guard<std::mutex> lock(queue_mutex);
        frontier_size = queue.size() + tasks_executing.load();
    }
    assert(frontier_size != 0);
    target_threads_per_task = std::max(target_threads_per_task, threads / frontier_size);
    
    uint64_t node_id, weight;
    bool got_next = false;
    while (!got_next) {
        // we are starting to face memory constraints, so we switch to large jobs with higher threads
        if (target_threads_per_task < threads) {
            // early phase of the scheduling when there are many small problems that use not-so-much memory
            // so we can devote threads to task parallelism
            queue_mutex.lock();
            auto next_iter = queue.begin();
            if (next_iter->first + memory_executing.load() < memory_limit) {
                std::tie(weight, node_id) = *next_iter;
                queue.erase(next_iter);
                ++tasks_executing;
                queue_mutex.unlock();
                got_next = true;
            }
            else {
                // we have hit a memory limit, increase the target threads per task
                queue_mutex.unlock();
                target_threads_per_task = std::min(target_threads_per_task + 1, threads);
                
                // wait until a task can fit in the memory limit (or we are as unconstrained as possible)
                // to avoid repeatedly hitting this condition and increasing when we repeat the outer loop
                auto can_fit_in_memory_limit = [this]() {
                    std::lock_guard<std::mutex> lock(this->queue_mutex);
                    auto mem_executing_now = this->memory_executing.load();
                    return this->queue.begin()->first < (this->memory_limit - mem_executing_now) || mem_executing_now == 0;
                };
                while (!can_fit_in_memory_limit()) {
                    std::this_thread::sleep_for(std::chrono::milliseconds(sleep_ms));
                }
            }
        }
        else if (tasks_executing.load() == 0) {
            // we will never have a better opportunity to execute the largest task, regardless of
            // its memory contribution
            queue_mutex.lock();
            auto next_iter = queue.end();
            --next_iter;
            std::tie(weight, node_id) = *next_iter;
            queue.erase(next_iter);
            ++tasks_executing;
            queue_mutex.unlock();
            got_next = true;
        }
        else {
            // find the largest task we can schedule within the budget
            int64_t memory_budget = memory_limit - memory_executing.load();
            queue_mutex.lock();
            auto next_iter = queue.upper_bound(memory_budget);
            if (next_iter != queue.begin()) {
                --next_iter;
                std::tie(weight, node_id) = *next_iter;
                queue.erase(next_iter);
                ++tasks_executing;
                queue_mutex.unlock();
                got_next = true;
            }
            else {
                // none are small enough, have to wait
                queue_mutex.unlock();
                std::this_thread::sleep_for(std::chrono::milliseconds(sleep_ms));
            }
        }
    }
    
    auto& task_info = scheduling_info[node_id];
    
    // check whether the chosen task is the only one that can run
    bool bottlenecked = (tasks_executing.load() == 1);
    queue_mutex.lock();
    task_info.queue_iter = queue.end(); // also grab this in the same lock
    bottlenecked = (bottlenecked && queue.empty());
    queue_mutex.unlock();
    
    if (bottlenecked) {
        // give all threads, irrespective of full ability to utilize them
        task_info.threads_assigned = threads_available.load();
    }
    else {
        // give the current task-per-thread, up to the tasks ability to use them and thread availability
        task_info.threads_assigned = std::min(target_threads_per_task, std::min(task_info.max_threads, threads_available.load()));
    }
    threads_available -= task_info.threads_assigned;
    memory_executing += weight;
//    std::cerr << ("scheduling task " + std::to_string(node_id) + " with " + std::to_string(task_info.threads_assigned) + " threads\n");
    return std::make_pair(node_id, task_info.threads_assigned);
}


void ParallelScheduler::handle_task(const std::function<uint64_t(void)>& task) {
    
    // do the task and then do bookkeeping when finished
    auto thread_task = [task, this]() {
        uint64_t node_id = task();
        mark_complete(node_id);
        auto& task_info = scheduling_info[node_id];
        --tasks_executing;
        memory_executing -= task_info.memory;
        threads_available += task_info.threads_assigned;
    };
    
    // launch and release the task
//    std::cerr << "launch task\n";
    std::thread(thread_task).detach();
}



void ParallelScheduler::mark_complete(uint64_t node_id) {

    if (subproblem_finished[node_id]) {
        return;
    }
    subproblem_finished[node_id] = true;
    
    auto& task_info = scheduling_info[node_id];
    {
        std::lock_guard<std::mutex> lock(queue_mutex);
        if (task_info.queue_iter != queue.end()) {
            queue.erase(task_info.queue_iter);
            task_info.queue_iter = queue.end();
        }
    }
    
    auto parent_id = task_info.parent;
    if (parent_id != -1) {
        auto& parent_info = scheduling_info[parent_id];
        --parent_info.children_remaining;
        if (parent_info.children_remaining == 0 && !subproblem_finished[parent_id]) {
//            std::cerr << ("adding " + std::to_string(parent_id) + " to queue\n");
            std::lock_guard<std::mutex> lock(queue_mutex);
            scheduling_info[parent_id].queue_iter = queue.emplace(parent_info.memory, parent_id);
        }
    }
}

}
