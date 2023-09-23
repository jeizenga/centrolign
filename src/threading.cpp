#include "centrolign/threading.hpp"

#include <vector>
#include <atomic>
#include <thread>

namespace centrolign {

int Threading::num_threads = 1;

int Threading::get_num_threads() {
    return num_threads;
}

void Threading::set_num_threads(int new_num_threads) {
    if (new_num_threads <= 0) {
        throw std::runtime_error("Invalid value for number of threads: " + std::to_string(new_num_threads));
    }
    num_threads = new_num_threads;
}

void Threading::parallel_for(size_t size, const std::function<void(size_t)> lambda) {
    
    std::atomic<size_t> iter(0);
    
    std::vector<std::thread> workers;
    workers.reserve(num_threads);
    for (int i = 0; i < num_threads; ++i) {
        workers.emplace_back([&iter,&lambda,&size](){
            size_t idx;
            while ((idx = iter.fetch_add(1)) < size) {
                lambda(idx);
            }
        });
    }
    
    for (auto& worker : workers) {
        worker.join();
    }
}

}
