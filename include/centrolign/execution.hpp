#ifndef centrolign_execution_hpp
#define centrolign_execution_hpp

#include <vector>
#include <string>
#include <tuple>
#include <functional>

#include "centrolign/graph.hpp"
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
    
    BaseGraph graph;
    SentinelTableau tableau;
    Alignment alignment;
    std::string name;
};

/*
 * An instance of an alignment step in the progressive alignment
 */
struct ProgressiveStep {
    ProgressiveStep(Subproblem* parent, Subproblem* child1, Subproblem* child2, uint64_t thread_budget) : parent(parent), child1(child1), child2(child2), thread_budget(thread_budget) {}
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
};


class SubproblemScheduler {
protected:
    
    SubproblemScheduler(uint64_t threads);
    
    SubproblemScheduler() = default;
    SubproblemScheduler(const SubproblemScheduler& other) = default;
    SubproblemScheduler(SubproblemScheduler&& other) = default;
    virtual ~SubproblemScheduler() = default;
    
    const uint64_t threads = 1;
    
public:
    // is iteration complete
    virtual bool finished() const = 0;
    
    // the next subproblem and its thread allocation
    virtual std::pair<uint64_t, uint64_t> next() = 0;
    
};


/*
 * Class that keeps track of the inputs, state, and execution order of the MSA
 */
class Execution {
protected:
    
    Execution(bool suppress_logging);
    
    Execution() = default;
    Execution(const Execution& other) = default;
    Execution(Execution&& other) = default;
    virtual ~Execution() = default;
    
    Execution& operator=(const Execution& other) = default;
    Execution& operator=(Execution&& other) = default;
    
public:
    
    void init(std::vector<std::pair<std::string, std::string>>&& names_and_sequences,
              Tree&& tree);
    
    // true when the entire execution has completed
    bool finished();
    
    // the next subproblem and its two children (parent first)
    ProgressiveStep next();
    
    // mark a subproblem as completed
    virtual void finish_subproblem(const Subproblem& subproblem);
    
    // TODO: do this inside the execution (need to proactively skip complete nexts during iteration)
    // true if a subproblem has already been aligned
    bool is_complete(const Subproblem& subproblem);
    
    // restart from saved partial results
    void restart(std::function<std::string(const Subproblem&)>& file_location,
                 bool preserve_leaves);
    
    // get the subproblem corresponding to one input sequence
    const Subproblem& leaf_subproblem(const std::string& name) const;
    
    // get the last subproblem, which contains the full alignment after execution
    Subproblem& final_subproblem();
    const Subproblem& final_subproblem() const;
    
    // get all of the subproblems that correspond to an input sequence
    std::vector<Subproblem*> leaf_subproblems();
    
    // get a hash identifier for a subproblem
    uint64_t subproblem_hash(const Subproblem& subproblem) const;
    
    // get the names of the sequences involved in a given subproblem
    std::vector<std::string> leaf_descendents(const Subproblem& subproblem) const;
    
    // TODO: i don't love exposing this, but it's currently needed to make the subpath
    // tree during cyclic graph polishing
    const Tree& get_tree() const;
    
    bool preserve_subproblems = true;
    
    // should tasks be carried out in parallel
    bool task_parallel = false;
    
    // total size of thread pool
    uint64_t threads = 1;
    
protected:
    
    uint64_t get_tree_id(const Subproblem& subproblem) const;
    
    SubproblemScheduler* get_scheduler();
    
    bool suppress_logging = false;
        
    // the guide tree
    Tree tree;
    
    // the individual alignment subproblems (including single-sequence leaves)
    std::vector<Subproblem> subproblems;
    
    std::vector<bool> subproblem_finished;
    
    std::shared_ptr<SubproblemScheduler> scheduler;
        
public:
    
    size_t memory_size() const;
};

class SerialScheduler : public SubproblemScheduler {
public:
    SerialScheduler(const Tree& tree, uint64_t threads);
    
    SerialScheduler() = default;
    SerialScheduler(const SerialScheduler& other) = default;
    SerialScheduler(SerialScheduler&& other) = default;
    ~SerialScheduler() = default;
    
    // is iteration complete
    bool finished() const;
    
    // the next subproblem and its thread allocation
    std::pair<uint64_t, uint64_t> next();
        
private:
    
    // a postorder of the inner nodes of the tree
    std::vector<uint64_t> execution_order;
    
    
    // the next subproblem in the execution order that we need to do
    size_t next_subproblem = 0;
    
};

class MainExecution : public virtual Execution {
public:
    
    MainExecution();
    MainExecution(const MainExecution& other) = default;
    MainExecution(MainExecution&& other) = default;
    ~MainExecution() = default;
    
    MainExecution& operator=(const MainExecution& other) = default;
    MainExecution& operator=(MainExecution&& other) = default;
    
    void finish_subproblem(const Subproblem& subproblem);
    
    std::string subalignments_filepath;
    
};


class MinorExecution : public virtual Execution {
public:
    
    MinorExecution();
    MinorExecution(const MinorExecution& other) = default;
    MinorExecution(MinorExecution&& other) = default;
    ~MinorExecution() = default;
    
    MinorExecution& operator=(const MinorExecution& other) = default;
    MinorExecution& operator=(MinorExecution&& other) = default;
    
private:
    
};


//class TaskSerialExecution : public virtual Execution {
//public:
//    
//    TaskSerialExecution(std::vector<std::pair<std::string, std::string>>&& names_and_sequences,
//                        Tree&& tree, bool suppress_logging = false);
//    TaskSerialExecution() = default;
//    TaskSerialExecution(const TaskSerialExecution& other) = default;
//    TaskSerialExecution(TaskSerialExecution&& other) = default;
//    ~TaskSerialExecution() = default;
//    
//    TaskSerialExecution& operator=(const TaskSerialExecution& other) = default;
//    TaskSerialExecution& operator=(TaskSerialExecution&& other) = default;
//    
//    
//    bool finished() const;
//    
//    // the next subproblem and its two children (parent first)
//    ProgressiveStep next();
//    
//    // FIXME: remove this when subalignments move into execution
//    std::tuple<const Subproblem*, const Subproblem*, const Subproblem*> current() const;
//    
////    // mark a subproblem as completed
////    void finish_subproblem(const Subproblem& subproblem);
//    
//private:
//    
//    // a postorder of the inner nodes of the tree
//    std::vector<uint64_t> execution_order;
//    
//    // the next subproblem in the execution order that we need to do
//    size_t next_subproblem = 0;
//    
//};

//class MainSerialExecution : public TaskSerialExecution, public MainExecution {
//public:
//    MainSerialExecution(std::vector<std::pair<std::string, std::string>>&& names_and_sequences,
//                        Tree&& tree);
//    MainSerialExecution() = default;
//    MainSerialExecution(const MainSerialExecution& other) = default;
//    MainSerialExecution(MainSerialExecution&& other) = default;
//    ~MainSerialExecution() = default;
//    
//    MainSerialExecution& operator=(const MainSerialExecution& other) = default;
//    MainSerialExecution& operator=(MainSerialExecution&& other) = default;
//};


//class TaskParallelExecution : public Execution {
//public:
//    
//    TaskParallelExecution(std::vector<std::pair<std::string, std::string>>&& names_and_sequences,
//                          Tree&& tree, uint64_t threads, bool suppress_logging = false) : Execution(std::move(names_and_sequences), std::move(tree), threads, suppress_logging) {}
//    
//    bool finished() const;
//    
//    ProgressiveStep next();
//    
//    // mark a subproblem as completed
//    void finish_subproblem(const Subproblem& subproblem);
//
//}

}
#endif /* centrolign_execution_hpp */
