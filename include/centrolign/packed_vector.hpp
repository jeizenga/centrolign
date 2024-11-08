#ifndef centrolign_packed_vector_hpp
#define centrolign_packed_vector_hpp

#include <stdexcept>
#include <string>
#include <cstdint>
#include <cstdlib>

#include "centrolign/vector_support.hpp"

namespace centrolign {

/*
 * A fixed size integer array that automatically and dynamically bit-compresses the entries
 */
class PackedArray {
public:
    PackedArray() noexcept = default;
    PackedArray(size_t size) noexcept;
    PackedArray(size_t size, uint8_t width) noexcept;
    PackedArray(const PackedArray& other) = delete;
    PackedArray(PackedArray&& other) noexcept;
    ~PackedArray() noexcept;
    
    PackedArray& operator=(const PackedArray& other) = delete;
    PackedArray& operator=(PackedArray&& other) noexcept;
    
    // get a value
    inline uint64_t at(size_t i) const;
    // set a value
    inline void set(size_t i, uint64_t value, size_t size);
    
    inline uint8_t width() const;
    
    using value_type = uint64_t;
    
private:
    
    inline static void set_internal(uint64_t*& array, const uint8_t& width, const size_t& i, const uint64_t& value);
    
    inline void repack(uint8_t new_width, size_t size);
    
    uint64_t* array = nullptr;
    uint8_t _width = 1;
    
public:
    size_t memory_size(size_t size) const;
    
};

/*
 * A fixed size integer vector that automatically and dynamically bit-compresses the entries
 */
class PackedVector {
public:
    PackedVector() noexcept = default;
    PackedVector(size_t size) noexcept;
    PackedVector(size_t size, uint8_t width) noexcept;
    PackedVector(const PackedVector& other) noexcept;
    PackedVector(PackedVector&& other) noexcept = default;
    ~PackedVector() noexcept = default;
    
    PackedVector& operator=(const PackedVector& other) noexcept;
    PackedVector& operator=(PackedVector&& other) noexcept = default;
    
    // get a value
    inline uint64_t at(size_t i) const;
    // set a value
    inline void set(size_t i, uint64_t value);
    // the length of the vector
    inline size_t size() const;
    // true if the vector has no entries
    inline bool empty() const;
    
    // get a value
    inline uint64_t operator[](size_t i) const;
    // get or set a value
    inline IntVectorSetter<PackedVector> operator[](size_t i);
    
    
    using value_type = uint64_t;
    
private:
    
    PackedArray array;
    size_t _size = 0;
    
public:
    size_t memory_size() const;
    
};

class SignedPackedVector {
public:
    SignedPackedVector() noexcept = default;
    SignedPackedVector(size_t size) noexcept;
    SignedPackedVector(size_t size, uint8_t width) noexcept;
    SignedPackedVector(const SignedPackedVector& other) noexcept = default;
    SignedPackedVector(SignedPackedVector&& other) noexcept = default;
    ~SignedPackedVector() noexcept = default;
    
    SignedPackedVector& operator=(const SignedPackedVector& other) noexcept = default;
    SignedPackedVector& operator=(SignedPackedVector&& other) noexcept = default;
    
    // get a value
    inline int64_t at(size_t i) const;
    // set a value
    inline void set(size_t i, int64_t value);
    // the length of the vector
    inline size_t size() const;
    // true if the vector has no entries
    inline bool empty() const;
    
    // get a value
    inline int64_t operator[](size_t i) const;
    // get or set a value
    inline IntVectorSetter<SignedPackedVector> operator[](size_t i);
    
    using value_type = int64_t;
private:
    
    inline static uint64_t encode(int64_t val);
    inline static int64_t decode(uint64_t val);
    
    PackedVector vec;
    
};


/*
 * Template and inline implementations
 */

inline uint64_t PackedVector::at(size_t i) const {
    if (i >= _size) {
        throw std::out_of_range("Out of bounds index " + std::to_string(i) + " in PackedVector of size " + std::to_string(_size));
    }
    return array.at(i);
}

inline void PackedVector::set(size_t i, uint64_t value) {
    if (i >= _size) {
        throw std::out_of_range("Out of bounds index " + std::to_string(i) + " in PackedVector of size " + std::to_string(_size));
    }
    array.set(i, value, _size);
}

inline size_t PackedVector::size() const {
    return _size;
}

inline bool PackedVector::empty() const {
    return _size == 0;
}

inline uint64_t PackedVector::operator[](size_t i) const {
    return array.at(i);
}

inline IntVectorSetter<PackedVector> PackedVector::operator[](size_t i) {
    return IntVectorSetter<PackedVector>(*this, i);
}

inline uint8_t PackedArray::width() const {
    return _width;
}

inline uint64_t PackedArray::at(size_t i) const {
    uint64_t mask = uint64_t(-1) >> (64 - _width);
    size_t b = i * _width;
    size_t k = b / 64;
    if (k == (b + _width - 1) / 64) {
        // whole value is in one 64-bit int
        size_t s = 64 - _width - (b % 64);
        return ((array[k] & (mask << s)) >> s);
    }
    else {
        // value is spread across two 64-bit ints
        size_t s1 = ((b % 64) + _width - 64);
        size_t s2 = (k + 2) * 64 - b - _width;
        return ((array[k] & (mask >> s1)) << s1) | ((array[k + 1] & (mask << s2)) >> s2);
    }
}

inline void PackedArray::set_internal(uint64_t*& array, const uint8_t& width, const size_t& i, const uint64_t& value) {
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

inline void PackedArray::set(size_t i, uint64_t value, size_t size) {
    uint8_t w = _width;
    while (w < 64 && (value & (uint64_t(-1) << w)) != 0) {
        ++w;
    }
    if (w != _width) {
        repack(w, size);
    }
    set_internal(array, _width, i, value);
}

inline void PackedArray::repack(uint8_t new_width, size_t size) {
    uint64_t* new_array = (uint64_t*) calloc(((size * new_width + 63) / 64), sizeof(uint64_t));
    for (size_t i = 0; i < size; ++i) {
        set_internal(new_array, new_width, i, at(i));
    }
    free(array);
    array = new_array;
    _width = new_width;
}

inline int64_t SignedPackedVector::at(size_t i) const {
    return decode(vec.at(i));
}

inline void SignedPackedVector::set(size_t i, int64_t value) {
    vec.set(i, encode(value));
}

inline size_t SignedPackedVector::size() const {
    return vec.size();
}

inline bool SignedPackedVector::empty() const {
    return vec.empty();
}

inline int64_t SignedPackedVector::operator[](size_t i) const {
    return at(i);
}

inline IntVectorSetter<SignedPackedVector> SignedPackedVector::operator[](size_t i) {
    return IntVectorSetter<SignedPackedVector>(*this, i);
}

inline uint64_t SignedPackedVector::encode(int64_t val) {
    if (val >= 0) {
        return val << 1;
    }
    else {
        return ((-val) << 1) | 1;
    }
}

inline int64_t SignedPackedVector::decode(uint64_t val) {
    if (val & 1) {
        return -(val >> 1);
    }
    else {
        return val >> 1;
    }
}

}

#endif /* centrolign_packed_vector_hpp */
