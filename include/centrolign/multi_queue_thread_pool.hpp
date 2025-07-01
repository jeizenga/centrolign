#ifndef centrolign_multi_queue_thread_pool_hpp
#define centrolign_multi_queue_thread_pool_hpp

#include <vector>
#include <list>
#include <queue>
#include <thread>
#include <atomic>
#include <functional>

#include "centrolign/thread_pool.hpp"

namespace centrolign {

class MultiQueueThreadPool; // forward declaration

/*
 * Token that represents a distinct queue managed by a thread pool
 */
class QueueToken {
public:
    ~QueueToken() = default;
private:
    using queue_ptr_t = std::list<std::pair<std::queue<std::function<void()>>, std::atomic<size_t>>>::iterator;
    QueueToken() = delete;
    QueueToken(queue_ptr_t it) : it(it) {}
    
    queue_ptr_t it;
    friend class MultiQueueThreadPool;
};


/*
 * A busy-waiting thread-pool for task-level parallelism across multiple
 * simultaneous queues
 */
class MultiQueueThreadPool {
public:
    
    MultiQueueThreadPool(size_t thread_count);
    MultiQueueThreadPool() = delete;
    ~MultiQueueThreadPool();
    
    // create a new queue and get the token to submit to that queue
    QueueToken init_queue();
    
    // submit a task to be done to the indicated queue
    void submit(QueueToken token, std::function<void()>&& task);
    
    // prioritize this queue and work on it until it is cleared. invalidates token
    void sync(QueueToken token);
    
    // work on tasks ans wait to clear all queues. invalidates all tokens
    void sync_all();
    
private:
    
    // move this task queue to the back of the queue of queues
    void rotate_out(std::list<std::pair<std::queue<std::function<void()>>, std::atomic<size_t>>>::iterator q);
    
    // move this task queue to the front of the queue of queues
    void rotate_in(std::list<std::pair<std::queue<std::function<void()>>, std::atomic<size_t>>>::iterator q);
    
    void do_tasks(QueueToken* token);
    
    std::vector<std::thread> workers;
    std::list<std::pair<std::queue<std::function<void()>>, std::atomic<size_t>>> queues;
    SpinLock task_lock;
    size_t num_queues = 0;
    std::atomic<bool> stop;
};

}
#endif /* centrolign_multi_queue_thread_pool_hpp */
