#include "centrolign/packed_vector.hpp"

namespace centrolign {


PackedVector::PackedVector(size_t s) noexcept : PackedVector(s, 1) {
    
}

PackedVector::PackedVector(size_t s, uint8_t w) noexcept : array(s, w), _size(s) {
    
}

PackedArray::PackedArray(size_t size) noexcept : PackedArray(size, 1) {
    
}

PackedArray::PackedArray(size_t size, uint8_t width) noexcept : _width(width) {
    array = (uint64_t*) calloc((size * width + 63) / 64, sizeof(uint64_t));
}

PackedArray::~PackedArray() noexcept {
    free(array);
}

PackedArray::PackedArray(PackedArray&& other) noexcept : _width(other._width), array(other.array) {
    other.array = nullptr;
    other._width = 1;
}

PackedArray& PackedArray::operator=(PackedArray&& other) noexcept {
    if (this != &other) {
        free(array);
        _width = other._width;
        array = other.array;
        other.array = nullptr;
        other._width = 1;
    }
    return *this;
}

size_t PackedArray::memory_size(size_t size) const {
    return sizeof(PackedArray) + ((size * _width + 63) / 64) * sizeof(uint64_t);
}

PackedVector::PackedVector(const PackedVector& other) noexcept : _size(other._size), array(other._size, other.array.width()) {
    for (size_t i = 0; i < _size; ++i) {
        array.set(i, other.array.at(i), _size);
    }
}

PackedVector& PackedVector::operator=(const PackedVector& other) noexcept {
    PackedArray new_array(other.size(), other.array.width());
    for (size_t i = 0; i < other.size(); ++i) {
        new_array.set(i, other.array.at(i), other.size());
    }
    array = std::move(new_array);
    _size = other._size;
    return *this;
}

size_t PackedVector::memory_size() const {
    return array.memory_size(_size) + sizeof(_size);
}

SignedPackedVector::SignedPackedVector(size_t size) noexcept : vec(size) {
    
}

SignedPackedVector::SignedPackedVector(size_t size, uint8_t width) noexcept : vec(size, width) {
    
}


}
