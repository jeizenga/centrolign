#ifndef centrolign_thread_pool_hpp
#define centrolign_thread_pool_hpp

#include <vector>
#include <queue>
#include <thread>
#include <atomic>
#include <functional>
#include <future>
#include <stdexcept>
#include <chrono>

namespace centrolign {

/*
 * A busy-waiting lock
 */
class SpinLock {
public:
    SpinLock() = default;
    ~SpinLock() = default;
    
    // spin until the lock becomes available
    void lock();
    // release the lock
    void unlock();
    
private:
    std::atomic_flag flag = ATOMIC_FLAG_INIT;
};

/*
 * A busy-waiting thread-pool for task-level parallelism
 */
class ThreadPool {
public:
    
    ThreadPool(size_t thread_count);
    ThreadPool() = delete;
    ~ThreadPool();
    
    
    
    // submit a task to be donw
    void submit(std::function<void()>&& task);
    
    // wait until the queue is cleared
    void sync();
    
private:
    
    void do_tasks(bool exit_on_empty);
    
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;
    
    SpinLock task_lock;
    std::atomic<bool> stop;
    std::atomic<size_t> checked_in;
};

}
#endif /* centrolign_thread_pool_hpp */
