#ifndef centrolign_paged_vector_hpp
#define centrolign_paged_vector_hpp

#include <vector>
#include <cstdint>

#include "centrolign/packed_vector.hpp"

namespace centrolign {


/*
 * A fixed size integer vector that automatically and dynamically bit-compresses and
 * windowed-difference compresses the entries (this is effective when values tend to
 * be locally correlated with each other
 */
template<size_t PageSize, size_t TiltBias = 1>
class PagedVector {
public:
    PagedVector() = default;
    PagedVector(size_t size);
    PagedVector(const PagedVector& other) = default;
    PagedVector(PagedVector&& other) = default;
    ~PagedVector() = default;
    
    PagedVector& operator=(const PagedVector& other) = default;
    PagedVector& operator=(PagedVector&& other) = default;
    
    // get a value
    inline uint64_t at(size_t i) const;
    // set a value
    inline void set(size_t i, uint64_t value);
    // the length of the vector
    inline size_t size() const;
    // true if the vector has no entries
    inline bool empty() const;
    
private:
    
    inline static uint64_t to_diff(const uint64_t& anchor, const uint64_t& value);
    inline static uint64_t from_diff(const uint64_t& anchor, const uint64_t& diff);
    
    size_t _size = 0;
    PackedVector anchors;
    std::vector<PackedVector> pages;
    
};




/*
 * Template and inline implementations
 */

template<size_t PageSize, size_t TiltBias>
PagedVector<PageSize, TiltBias>::PagedVector(size_t size) : _size(size), anchors((size + PageSize - 1) / PageSize), pages(anchors.size(), PackedVector(PageSize)){
    
}

template<size_t PageSize, size_t TiltBias>
inline uint64_t PagedVector<PageSize, TiltBias>::at(size_t i) const {
    size_t k = i / PageSize;
    return from_diff(anchors.at(k), pages[k].at(i % PageSize));
}

template<size_t PageSize, size_t TiltBias>
inline void PagedVector<PageSize, TiltBias>::set(size_t i, uint64_t value) {
    size_t k = i / PageSize;
    if (anchors.at(k) == 0) {
        // note: only 0 values can already be found in this page, and they stay valid because of
        // of the difference encoding
        anchors.set(k, value);
    }
    pages[k].set(i % PageSize, to_diff(anchors.at(k), value));
}

template<size_t PageSize, size_t TiltBias>
inline size_t PagedVector<PageSize, TiltBias>::size() const  {
    return _size;
}

template<size_t PageSize, size_t TiltBias>
inline bool PagedVector<PageSize, TiltBias>::empty() const  {
    return _size == 0;
}

template<size_t PageSize, size_t TiltBias>
inline uint64_t PagedVector<PageSize, TiltBias>::to_diff(const uint64_t& anchor, const uint64_t& value) {
    if (value == 0) {
        return 0;
    }
    else if (value > anchor) {
        return ((value - anchor - 1) * (TiltBias + 1)) / TiltBias + 2;
    }
    else {
        return (anchor - value) * (TiltBias + 1) + 1;
    }
}

template<size_t PageSize, size_t TiltBias>
inline uint64_t PagedVector<PageSize, TiltBias>::from_diff(const uint64_t& anchor, const uint64_t& diff) {
    if (diff == 0){
        return 0;
    }
    else if ((diff - 1) % (TiltBias + 1) != 0) {
        return anchor + ((diff - 1) * TiltBias) / (TiltBias + 1) + 1;
    }
    else {
        return anchor - (diff - 1) / (TiltBias + 1);
    }
}

}

#endif /* centrolign_paged_vector_hpp */
