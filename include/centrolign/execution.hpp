#ifndef centrolign_execution_hpp
#define centrolign_execution_hpp

#include <vector>
#include <string>
#include <tuple>
#include <map>
#include <set>
#include <functional>
#include <atomic>
#include <thread>
#include <mutex>
#include <limits>

#include "centrolign/graph.hpp"
#include "centrolign/compressible_graph.hpp"
#include "centrolign/modify_graph.hpp"
#include "centrolign/alignment.hpp"
#include "centrolign/tree.hpp"

namespace centrolign {

/*
 * An alignment resulting from a subproblem in the process of making the full MSA
 */
struct Subproblem {
    Subproblem() noexcept = default;
    Subproblem(const Subproblem& other) noexcept = default;
    Subproblem(Subproblem&& other) noexcept = default;
    ~Subproblem() = default;
    Subproblem& operator=(const Subproblem& other) noexcept = default;
    Subproblem& operator=(Subproblem&& other) noexcept = default;
    
    // the graph and its sentinel tableau, which may be compressed while not in use
    CompressibleGraph graph;
    Alignment alignment;
    std::string name;
};

/*
 * An instance of an alignment step in the progressive alignment
 */
struct ProgressiveStep {
    ProgressiveStep(Subproblem* parent, Subproblem* child1, Subproblem* child2, uint64_t thread_budget, bool restrain_memory = false) : parent(parent), child1(child1), child2(child2), thread_budget(thread_budget), restrain_memory(restrain_memory) {}
    ProgressiveStep() noexcept = default;
    ProgressiveStep(const ProgressiveStep& other) noexcept = default;
    ProgressiveStep(ProgressiveStep&& other) noexcept = default;
    ~ProgressiveStep() = default;
    
    ProgressiveStep& operator=(const ProgressiveStep& other) noexcept = default;
    ProgressiveStep& operator=(ProgressiveStep&& other) noexcept = default;
    
    Subproblem* parent = nullptr;
    Subproblem* child1 = nullptr;
    Subproblem* child2 = nullptr;
    uint64_t thread_budget = 0;
    // use memory-efficient data structures because of the memory use of concurrently executing steps
    bool restrain_memory = false;
};


/*
 * An abstract queue for subproblems
 */
class SubproblemScheduler {
protected:
    
    SubproblemScheduler(const Tree& tree, uint64_t threads);
    
    SubproblemScheduler() = default;
    virtual ~SubproblemScheduler() = default;
    
    const uint64_t threads = 1;
    
    std::vector<bool> subproblem_finished;
    
public:
    // is iteration complete
    virtual bool finished() = 0;
    
    // the next subproblem, its thread allocation, and whether it should restrain its memory use
    virtual std::tuple<uint64_t, uint64_t, bool> next() = 0;
    
    // launch a task and do necessary bookkeeping
    virtual void handle_task(const std::function<uint64_t(void)>& task) = 0;
    
    // indicate that a subproblem has been completed
    virtual void mark_complete(uint64_t node_id) = 0;
    
    // has the subproblem been marked complete
    bool is_complete(uint64_t node_id) const;
    
    // the maximum total of (number of sequences 1) * (number of sequences 2) across all concurrently
    // executing subproblems before subproblems will be told to restrain memory
    double memory_restraint_budget = std::numeric_limits<double>::max();
};


/*
 * Class that keeps track of the inputs, state, and execution order of the MSA
 */
class Execution {
protected:
    
    Execution(bool suppress_logging);
    
    Execution() = default;
    virtual ~Execution() = default;
    
public:
    
    void init(std::vector<std::pair<std::string, std::string>>&& names_and_sequences,
              Tree&& tree);
    
    // the next subproblem and its two children (parent first)
    virtual void execute(const std::function<void(const ProgressiveStep&)>& do_subproblem);
    
    // get the subproblem corresponding to one input sequence
    const Subproblem& leaf_subproblem(const std::string& name) const;
    
    // get the last subproblem, which contains the full alignment after execution
    Subproblem& final_subproblem();
    const Subproblem& final_subproblem() const;
    
    // get all of the subproblems that correspond to an input sequence
    std::vector<Subproblem*> leaf_subproblems();
    
    // get the names of the sequences involved in a given subproblem
    std::vector<std::string> leaf_descendents(const Subproblem& subproblem) const;
    
    // TODO: i don't love exposing this, but it's currently needed to make the subpath
    // tree during cyclic graph polishing
    const Tree& get_tree() const;
    
    // should we keep intermediate subproblems in memory?
    bool preserve_subproblems = true;
    
    // should tasks be carried out in parallel
    bool task_parallel = false;
    
    // total size of thread pool
    uint64_t threads = 1;
    
    // the maximum total of (number of sequences 1) * (number of sequences 2) across all concurrently
    // executing subproblems before subproblems will be told to restrain memory
    double memory_restraint_budget = std::numeric_limits<double>::max();
    
protected:
    
    uint64_t get_tree_id(const Subproblem& subproblem) const;
    
    SubproblemScheduler* get_scheduler();
    
    bool suppress_logging = false;
        
    // the guide tree
    Tree tree;
    
    // the individual alignment subproblems (including single-sequence leaves)
    std::vector<Subproblem> subproblems;
    
    std::shared_ptr<SubproblemScheduler> scheduler;
        
public:
    
    size_t memory_size() const;
};

/*
 * The main execution of the algorithm
 */
class MainExecution : public virtual Execution {
public:
    
    MainExecution();
    ~MainExecution() = default;
    
    void execute(const std::function<void(const ProgressiveStep&)>& do_subproblem);
    
    // restart from saved partial results
    void restart(bool preserve_leaves);
    
    // if non-empty, file to write suproblem alignments to
    std::string subalignments_filepath;
    
    // if non-empty, prefix to give GFA output for all suproblems
    std::string subproblems_prefix;
    
    // deterministic high entropy file name to avoid collisions
    std::string subproblem_file_name(const Subproblem& subproblem) const;
    
private:
    
    // emit a subproblem GFA
    void emit_subproblem(const Subproblem& subproblem);
    
    // name for file to map file names to sample sets
    std::string subproblem_info_file_name() const;
    
    // get a hash identifier for a subproblem
    uint64_t subproblem_hash(const Subproblem& subproblem) const;
    
    std::mutex subalignments_mutex;
    std::mutex subproblem_info_mutex;
    
};

/*
 * Sub-executions while harmonizing cyclic motifs
 */
class MinorExecution : public virtual Execution {
public:
    
    MinorExecution();
    ~MinorExecution() = default;
        
};

/*
 * Helper class for scheduling subproblemms in task-serial (tasks may still be
 * internally multithreaded)
 */
class SerialScheduler : public SubproblemScheduler {
public:
    SerialScheduler(const Tree& tree, uint64_t threads);
    
    SerialScheduler() = default;
    ~SerialScheduler() = default;
    
    // is iteration complete
    bool finished();
    
    // the next subproblem, its thread allocation, and whether it should restrain its memory use
    std::tuple<uint64_t, uint64_t, bool> next();

    // launch a task and do necessary bookkeeping
    void handle_task(const std::function<uint64_t(void)>& task);
    
    // indicate that a subproblem has been completed
    void mark_complete(uint64_t node_id);
        
private:
    
    // a postorder of the inner nodes of the tree
    std::vector<uint64_t> execution_order;
    
    // the next subproblem in the execution order that we need to do
    size_t next_subproblem = 0;
    
};


/*
 * Helper class for scheduling subproblems in task-parallel
 */
class ParallelScheduler : public SubproblemScheduler {
public:
    ParallelScheduler(const Tree& tree, uint64_t threads, size_t problem_size);
    
    ParallelScheduler() = default;
    ~ParallelScheduler() = default;
    
    // is iteration complete
    bool finished();
    
    // the next subproblem, its thread allocation, and whether it should restrain its memory use
    std::tuple<uint64_t, uint64_t, bool> next();
    
    // launch a task and do necessary bookkeeping, task should return subproblem ID on completion
    void handle_task(const std::function<uint64_t(void)>& task);
    
    // indicate that a subproblem has been completed
    void mark_complete(uint64_t node_id);
    
private:
    
    struct SchedulingInfo {
        // memory use (arbitrary units)
        int64_t memory = 0;
        // max number of fully utilizable threads
        uint64_t max_threads = 1;
        // once scheduled, the number of threads assigned to the task
        uint64_t threads_assigned = 0;
        // once scheduled, the units of memory charged against the memory restraint budget
        int64_t restraint_charge = 0;
        // number of child tasks that have not yet completed
        uint64_t children_remaining = -1;
        // parent node ID (for signaling completion)
        uint64_t parent = -1;
        // the minimum memory of this subproblem and all of its ancestors
        int64_t min_memory_to_root = 0;
        // whether the subproblem is in the frontier (either queued or executing)
        bool in_frontier = false;
        // once queued, the iterator 
        std::multimap<int64_t, uint64_t>::iterator queue_iter;
    };
    
    static constexpr double max_relative_memory = 0.8;
    // approximate ratio of the memory use of unrestrained to restrained subproblems
    static constexpr int64_t unrestrained_memory_ratio = 2;
    const uint64_t sleep_ms = 1;
    
    std::vector<SchedulingInfo> scheduling_info;
    
    // queue with priority determined by memory footprint (memory -> node ID)
    std::multimap<int64_t, uint64_t> queue;
    // the minimum memory to the root of each subproblem in the frontier (queued or executing),
    // every future subproblem is an ancestor of one of these
    std::multiset<int64_t> frontier_min_memory;

    // the number of threads we are currently giving to small tasks
    uint64_t target_threads_per_task = 1;
    // lock for interacting with queue
    std::mutex queue_mutex;
    // number of free threads that are not executing tasks
    std::atomic<uint64_t> threads_available;
    // number of tasks currently executing
    std::atomic<uint64_t> tasks_executing;
    // the maximum units of memory that we will try to have in use at one time
    int64_t memory_limit = 0;
    // the units of memory currently executing
    std::atomic<int64_t> memory_executing;
    // the units of memory currently executing, as measured for the memory restraint budget (in which
    // restrained tasks are charged their weight and unrestrained tasks are charged the ratio times their weight)
    std::atomic<int64_t> restraint_memory_executing;
    
};

}
#endif /* centrolign_execution_hpp */
