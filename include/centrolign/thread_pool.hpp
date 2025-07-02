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

class AbstractPool {
public:
    
    virtual ~AbstractPool() = default;
    
    // submit a task to be done
    virtual void submit(std::function<void()>&& task) = 0;
    
    // wait until the queue is cleared
    virtual void sync() = 0;
};


/*
 * A busy-waiting thread-pool for task-level parallelism
 */
class ThreadPool : public AbstractPool {
public:
    
    ThreadPool(size_t thread_count);
    ThreadPool() = delete;
    ~ThreadPool();
    
    // submit a task to be done
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

class SerialPool : public AbstractPool {
public:
    
    SerialPool() = default;
    ~SerialPool() = default;
    
    
    // submit a task to be done
    void submit(std::function<void()>&& task);
    
    // wait until the queue is cleared
    void sync();
};

}
#endif /* centrolign_thread_pool_hpp */
