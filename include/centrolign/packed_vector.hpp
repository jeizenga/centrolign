#ifndef centrolign_packed_vector_hpp
#define centrolign_packed_vector_hpp

#include <stdexcept>
#include <string>
#include <cstdint>
#include <cstdlib>

namespace centrolign {


/*
 * A fixed size integer vector that automatically and dynamically bit-compresses the entries
 */
class PackedVector {
public:
    PackedVector() noexcept = default;
    PackedVector(size_t size) noexcept;
    PackedVector(size_t size, uint8_t width) noexcept;
    PackedVector(const PackedVector& other) noexcept;
    PackedVector(PackedVector&& other) noexcept;
    ~PackedVector() noexcept;
    
    PackedVector& operator=(const PackedVector& other) noexcept;
    PackedVector& operator=(PackedVector&& other) noexcept;
    
    // get a value
    inline uint64_t at(size_t i) const;
    // set a value
    inline void set(size_t i, uint64_t value);
    // the length of the vector
    inline size_t size() const;
    // true if the vector has no entries
    inline bool empty() const;
    
private:

    inline static void set_internal(uint64_t*& array, const uint8_t& width, const size_t& i, const uint64_t& value);
    inline static size_t internal_length(const size_t& s, const uint8_t& w);
    
    void repack(uint8_t new_width);
    
    uint64_t* array = nullptr;
    size_t _size = 0;
    uint8_t width = 1;
    
public:
    size_t memory_size() const;
    
};




/*
 * Template and inline implementations
 */

inline uint64_t PackedVector::at(size_t i) const {
    if (i >= size()) {
        throw std::out_of_range("Out of bounds index " + std::to_string(i) + " in PackedVector of size " + std::to_string(_size));
    }
    uint64_t mask = uint64_t(-1) >> (64 - width);
    size_t b = i * width;
    size_t k = b / 64;
    if (k == (b + width - 1) / 64) {
        // whole value is in one 64-bit int
        size_t s = 64 - width - (b % 64);
        return ((array[k] & (mask << s)) >> s);
    }
    else {
        // value is spread across two 64-bit ints
        size_t s1 = ((b % 64) + width - 64);
        size_t s2 = (k + 2) * 64 - b - width;
        return ((array[k] & (mask >> s1)) << s1) | ((array[k + 1] & (mask << s2)) >> s2);
    }
}

inline void PackedVector::set_internal(uint64_t*& array, const uint8_t& width, const size_t& i, const uint64_t& value) {
    // TODO: repetitive with ::at()
    uint64_t mask = uint64_t(-1) >> (64 - width);
    size_t b = i * width;
    size_t k = b / 64;
    if (k == (b + width - 1) / 64) {
        // whole value is in one 64-bit int
        size_t s = 64 - width - (b % 64);
        array[k] = (array[k] & ~(mask << s)) | ((value & mask) << s);
    }
    else {
        // value is spread across two 64-bit ints
        size_t s1 = ((b % 64) + width - 64);
        array[k] = (array[k] & ~(mask >> s1)) | ((mask & value) >> s1);
        size_t s2 = (k + 2) * 64 - b - width;
        array[k + 1] = (array[k + 1] & ~(mask << s2)) | ((mask & value) << s2);
    }
}

inline void PackedVector::set(size_t i, uint64_t value) {
    if (i >= size()) {
        throw std::out_of_range("Out of bounds index " + std::to_string(i) + " in PackedVector of size " + std::to_string(_size));
    }
    uint8_t w = width;
    while (w < 64 && (value & (uint64_t(-1) << w)) != 0) {
        ++w;
    }
    if (w != width) {
        repack(w);
    }
    set_internal(array, width, i, value);
}

inline size_t PackedVector::size() const {
    return _size;
}

inline bool PackedVector::empty() const {
    return _size == 0;
}

inline size_t PackedVector::internal_length(const size_t& s, const uint8_t& w) {
    return (s * w + 63) / 64;
}


}

#endif /* centrolign_packed_vector_hpp */
