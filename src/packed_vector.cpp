#include "centrolign/packed_vector.hpp"

namespace centrolign {


PackedVector::PackedVector(size_t s) noexcept : PackedVector(s, 1) {
}

PackedVector::PackedVector(size_t s, uint8_t width) noexcept : width(width), _size(s) {
    array = (uint64_t*) calloc(internal_length(s, width), sizeof(uint64_t));
}

PackedVector::~PackedVector() noexcept {
    free(array);
}

PackedVector::PackedVector(const PackedVector& other) noexcept : _size(other._size), width(other.width) {
    size_t l = internal_length(_size, width);
    array = (uint64_t*) malloc(l * sizeof(uint64_t));
    for (size_t i = 0; i < l; ++i) {
        array[i] = other.array[i];
    }
}

PackedVector::PackedVector(PackedVector&& other) noexcept : _size(other._size), width(other.width) {
    array = other.array;
    other.array = nullptr;
    other._size = 0;
}

PackedVector& PackedVector::operator=(const PackedVector& other) noexcept {
    free(array);
    _size = other._size;
    width = other.width;
    size_t l = internal_length(_size, width);
    array = (uint64_t*) malloc(l * sizeof(uint64_t));
    for (size_t i = 0; i < l; ++i) {
        array[i] = other.array[i];
    }
    return *this;
}

PackedVector& PackedVector::operator=(PackedVector&& other) noexcept {
    free(array);
    _size = other._size;
    width = other.width;
    array = other.array;
    other.array = nullptr;
    other._size = 0;
    return *this;
}


void PackedVector::repack(uint8_t new_width) {
    uint64_t* new_array = (uint64_t*) calloc(internal_length(_size, new_width), sizeof(uint64_t));
    for (size_t i = 0; i < size(); ++i) {
        set_internal(new_array, new_width, i, at(i));
    }
    free(array);
    array = new_array;
    width = new_width;
}

size_t PackedVector::memory_size() const {
    return sizeof(PackedVector) + internal_length(_size, width) * sizeof(uint64_t);
}

}
