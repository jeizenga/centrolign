#ifndef centrolign_utility_hpp
#define centrolign_utility_hpp

#include <vector>
#include <cstdint>

namespace centrolign {

// return a vector of 0,1,...,size-1
std::vector<size_t> range_vector(size_t size);

// returns the highest 1-bit (or 0 in case of 0)
uint32_t hi_bit(uint64_t x);

}

#endif /* centrolign_utility_hpp */
