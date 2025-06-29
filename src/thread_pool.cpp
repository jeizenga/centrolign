#include "centrolign/thread_pool.hpp"

namespace centrolign {

void SpinLock::lock() {
    while (flag.test_and_set(std::memory_order_acquire)) {
        std::this_thread::yield();
    }
}

void SpinLock::unlock() {
    flag.clear(std::memory_order_release);
}

ThreadPool::ThreadPool(size_t thread_count) : stop(false) {
    // make workers that do tasks forever
    for (size_t i = 0; i < thread_count; ++i) {
        workers.emplace_back([this]() {
            do_tasks(false);
        });
    }
}

ThreadPool::~ThreadPool() {
    // flag the workers to exit their infinite loop
    stop.store(true);
    // and wait for them to do it
    for (auto& worker : workers) {
        worker.join();
    }
}

void ThreadPool::do_tasks(bool exit_on_empty) {
    while (!stop.load()) {
        std::function<void()> task;
        
        // try to get a task from the queue
        {
            task_lock.lock();
            if (!tasks.empty()) {
                task = std::move(tasks.front());
                tasks.pop();
                task_lock.unlock();
            }
            else {
                task_lock.unlock();
                if (exit_on_empty) {
                    break;
                }
            }
        }
        
        if (task) {
            // do the task
            task();
        }
        else {
            // no task, spin briefly
            std::this_thread::yield();
        }
    }
}

void ThreadPool::submit(std::function<void()>&& task) {
    task_lock.lock();
    tasks.emplace(std::move(task));
    task_lock.unlock();
}

void ThreadPool::sync() {
    // do tasks in the leader thread until the queue empties
    do_tasks(true);
}


}
