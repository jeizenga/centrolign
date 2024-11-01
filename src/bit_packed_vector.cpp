#include "centrolign/bit_packed_vector.hpp"

namespace centrolign {


PackedVector::PackedVector(size_t s) : PackedVector(s, 1) {
}

PackedVector::PackedVector(size_t s, uint8_t width) : width(width), _size(s) {
    array = (uint64_t*) calloc((s * width + 63) / 64, sizeof(uint64_t));
}

PackedVector::~PackedVector() {
    free(array);
}

void PackedVector::repack(uint8_t new_width) {
    uint64_t* new_array = (uint64_t*) calloc((size() * new_width + 63) / 64, sizeof(uint64_t));
    for (size_t i = 0; i < size(); ++i) {
        set_internal(new_array, new_width, i, at(i));
    }
    free(array);
    array = new_array;
    width = new_width;
}

}
