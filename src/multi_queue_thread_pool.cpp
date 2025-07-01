#include "centrolign/multi_queue_thread_pool.hpp"

namespace centrolign {

MultiQueueThreadPool::MultiQueueThreadPool(size_t thread_count) : stop(false){
    // make workers that do tasks forever (on no particualr queue)
    for (size_t i = 0; i < thread_count; ++i) {
        workers.emplace_back([this]() {
            do_tasks(nullptr);
        });
    }
}

MultiQueueThreadPool::~MultiQueueThreadPool() {
    // finish all queues
    sync_all();
    // flag the workers to exit their infinite loop
    stop.store(true);
    // and wait for them to do it
    for (auto& worker : workers) {
        worker.join();
    }
}

QueueToken MultiQueueThreadPool::init_queue() {
    // make a new queue
    task_lock.lock();
    auto it = queues.emplace(queues.end());
    it->second.store(0);
    task_lock.unlock();
    ++num_queues;
    return QueueToken(it);
}

void MultiQueueThreadPool::rotate_out(std::list<std::pair<std::queue<std::function<void()>>, std::atomic<size_t>>>::iterator q) {
    queues.splice(queues.end(), queues, q);
}

void MultiQueueThreadPool::rotate_in(std::list<std::pair<std::queue<std::function<void()>>, std::atomic<size_t>>>::iterator q) {
    queues.splice(queues.begin(), queues, q);
}

void MultiQueueThreadPool::do_tasks(QueueToken* token) {
    // until we tear down the thread pool
    while (!stop.load()) {
        
        // try to get a task from a queue
        std::function<void()> task;
        std::atomic<size_t>* check_in = nullptr;
        {
            task_lock.lock();
            if (queues.empty()) {
                // there are no active queues
                task_lock.unlock();
            }
            else {
                // get either the first queue or a specific one
                auto it = token ? token->it : queues.begin();
                if (!it->first.empty()) {
                    // there are tasks in this queue, check one of them out
                    ++it->second;
                    task = std::move(it->first.front());
                    it->first.pop();
                    task_lock.unlock();
                    check_in = &it->second;
                }
                else {
                    // there are no tasks in this queue, deprioritize this queue
                    rotate_out(it);
                    task_lock.unlock();
                    if (token) {
                        // we were focusing on one queue, which is empty now, we
                        // can exit
                        return;
                    }
                }
            }
        }
        
        if (task) {
            // do the task
            task();
            --(*check_in);
        }
        else {
            // no task, spin briefly
            std::this_thread::yield();
        }
    }
}

void MultiQueueThreadPool::submit(QueueToken token, std::function<void()>&& task) {
    task_lock.lock();
    token.it->first.emplace(std::move(task));
    task_lock.unlock();
}

void MultiQueueThreadPool::sync(QueueToken token) {
    // prioritize this queue
    task_lock.lock();
    rotate_in(token.it);
    task_lock.unlock();
    // work on it until it's emptied
    do_tasks(&token);
    // wait for outstanding tasks in this queue to complete
    while (token.it->second.load() != 0) {
        std::this_thread::yield();
    }
    // clear this queue
    task_lock.lock();
    queues.erase(token.it);
    task_lock.unlock();
    --num_queues;
}

void MultiQueueThreadPool::sync_all() {
    while (num_queues != 0) {
        // get the next queue's token
        task_lock.lock();
        QueueToken token(queues.begin());
        task_lock.unlock();
        // clear it out
        sync(token);
    }
}


}
