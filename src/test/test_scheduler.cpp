#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <random>
#include <iostream>
#include <atomic>
#include <mutex>
#include <set>
#include <thread>
#include <chrono>
#include <functional>

#include "centrolign/execution.hpp"
#include "centrolign/tree.hpp"

using namespace std;
using namespace centrolign;

int leaf_counter = 0;

string random_newick(int num_leaves, default_random_engine& gen) {
    if (num_leaves == 1) {
        return "L" + to_string(leaf_counter++);
    }
    int left = uniform_int_distribution<int>(1, num_leaves - 1)(gen);
    return "(" + random_newick(left, gen) + "," + random_newick(num_leaves - left, gen) + ")";
}

string caterpillar_newick(int num_leaves) {
    string newick = "L" + to_string(leaf_counter++);
    for (int i = 1; i < num_leaves; ++i) {
        newick = "(" + newick + ",L" + to_string(leaf_counter++) + ")";
    }
    return newick;
}

string balanced_newick(int num_leaves) {
    if (num_leaves == 1) {
        return "L" + to_string(leaf_counter++);
    }
    return "(" + balanced_newick(num_leaves / 2) + "," + balanced_newick(num_leaves - num_leaves / 2) + ")";
}

// execute the tree with fake tasks that check their thread allocation at the start of several stages,
// returns whether any task saw its allocation increase while it executed
bool test_parallel_scheduler(const string& newick, uint64_t threads, default_random_engine& gen) {

    Tree tree(newick + ";");
    ParallelScheduler scheduler(tree, threads, 0);
    scheduler.suppress_logging = true;

    vector<atomic<int>> done(tree.node_size());
    size_t num_internal = 0;
    for (uint64_t node_id = 0; node_id < tree.node_size(); ++node_id) {
        done[node_id] = tree.is_leaf(node_id) ? 1 : 0;
        num_internal += !tree.is_leaf(node_id);
    }

    mutex running_mutex;
    set<uint64_t> running;
    atomic<bool> failed(false);
    atomic<bool> saw_increase(false);

    auto check_allocations = [&]() {
        // the executing tasks should never be allocated more threads than exist in total
        lock_guard<mutex> lock(running_mutex);
        uint64_t total = 0;
        for (auto node_id : running) {
            total += scheduler.thread_grant(node_id)->load();
        }
        if (total > threads) {
            cerr << "error: " << total << " threads allocated to executing tasks with " << threads << " threads available\n";
            failed = true;
        }
    };

    size_t launched = 0;
    while (!scheduler.finished()) {
        uint64_t node_id, task_threads;
        bool restrain;
        tie(node_id, task_threads, restrain) = scheduler.next();
        ++launched;

        if (task_threads == 0 || task_threads > threads || done[node_id]) {
            cerr << "error: invalid assignment of " << task_threads << " threads to node " << node_id << '\n';
            failed = true;
        }
        for (auto child_id : tree.get_children(node_id)) {
            if (!done[child_id]) {
                cerr << "error: node " << node_id << " scheduled before its child " << child_id << '\n';
                failed = true;
            }
        }
        {
            lock_guard<mutex> lock(running_mutex);
            running.insert(node_id);
        }
        check_allocations();

        vector<int> stage_ms(3);
        for (auto& ms : stage_ms) {
            ms = uniform_int_distribution<int>(0, 3)(gen);
        }
        const atomic<uint64_t>* grant = scheduler.thread_grant(node_id);
        scheduler.handle_task([&, node_id, stage_ms, task_threads, grant]() {
            uint64_t prev_threads = task_threads;
            for (auto ms : stage_ms) {
                // a parallel section checks its allocation when it starts
                uint64_t stage_threads = grant->load();
                if (stage_threads < prev_threads) {
                    cerr << "error: thread allocation of node " << node_id << " decreased while executing\n";
                    failed = true;
                }
                if (stage_threads > prev_threads) {
                    saw_increase = true;
                }
                prev_threads = stage_threads;
                check_allocations();
                this_thread::sleep_for(chrono::milliseconds(ms));
            }
            {
                lock_guard<mutex> lock(running_mutex);
                running.erase(node_id);
            }
            done[node_id] = 1;
            return node_id;
        });
    }

    if (launched != num_internal) {
        cerr << "error: launched " << launched << " tasks for " << num_internal << " subproblems\n";
        failed = true;
    }
    for (uint64_t node_id = 0; node_id < tree.node_size(); ++node_id) {
        if (!done[node_id]) {
            cerr << "error: node " << node_id << " was not completed\n";
            failed = true;
        }
    }
    if (failed) {
        cerr << "failed on tree " << newick << " with " << threads << " threads\n";
        exit(1);
    }
    return saw_increase;
}

int main(int argc, char* argv[]) {

    random_device rd;
    default_random_engine gen(rd());

    bool saw_increase = false;
    for (int num_leaves : {2, 3, 4, 7, 16, 33, 64}) {
        for (int shape = 0; shape < 3; ++shape) {
            for (uint64_t threads : {1, 2, 3, 8, 13}) {
                leaf_counter = 0;
                string newick;
                if (shape == 0) {
                    newick = random_newick(num_leaves, gen);
                }
                else if (shape == 1) {
                    newick = caterpillar_newick(num_leaves);
                }
                else {
                    newick = balanced_newick(num_leaves);
                }
                saw_increase |= test_parallel_scheduler(newick, threads, gen);
            }
        }
    }

    if (!saw_increase) {
        cerr << "error: no task's thread allocation was increased while it executed\n";
        return 1;
    }

    cerr << "passed all tests!" << endl;
    return 0;
}
