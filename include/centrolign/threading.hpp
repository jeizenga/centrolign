#ifndef centrolign_threading_hpp
#define centrolign_threading_hpp

#include <functional>
#include <cstdint>

namespace centrolign {

/*
 * Struct namespace for multithreading functions
 */
struct Threading {
public:
    
    static int get_num_threads();
    static void set_num_threads(int new_num_threads);
    static void parallel_for(size_t size, const std::function<void(size_t)> lambda);
    
private:
    
    static int num_threads;
};

}

#endif /* centrolign_threading_hpp */
