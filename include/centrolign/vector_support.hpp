#ifndef centrolign_vector_support_hpp
#define centrolign_vector_support_hpp

#include <cstdint>
#include <stdexcept>
#include <utility>

#include "centrolign/utility.hpp"

namespace centrolign {

/*
 * Functions as a reference for integer vector types for the purpose of operator[]
 */
template<class IntVectorType>
class IntVectorSetter {
public:
    
    IntVectorSetter(IntVectorType& vec, size_t i) : i(i), vec(&vec) {}
    IntVectorSetter() = default;
    ~IntVectorSetter() = default;
    
    void operator++() {
        vec->set(i, vec->at(i) + 1);
    }
    
    void operator--() {
        vec->set(i, vec->at(i) - 1);
    }
    
    template<typename T>
    inline void operator=(const T& value) {
        vec->set(i, typename IntVectorType::value_type(value));
    }
    
    template<typename T>
    operator T() const {
        return T(vec->at(i));
    }
    
private:
    
    IntVectorType* vec = nullptr;
    size_t i = 0;
};

/*
 * Functions as a reference for a pair of vector being used as a vector of pairs
 */
template<class VectorPairType>
class VectorPairSetter {
    
public:
    VectorPairSetter(VectorPairType& pair_vec, size_t i) : pair_vec(&pair_vec), i(i) { }
    VectorPairSetter() = default;
    ~VectorPairSetter() = default;
    
    template<typename Pair>
    inline void operator=(const Pair& value) {
        pair_vec->vec1[i] = value.first;
        pair_vec->vec2[i] = value.second;
    }
    
    template<typename T>
    operator T() const {
        return T(pair_vec->vec1.at(i), pair_vec->vec2.at(i));
    }
    
private:
    
    VectorPairType* pair_vec = nullptr;
    size_t i = 0;
};

/*
 * A pair of vector-like containers being used as a vector of pairs
 */
template<class Vec1, class Vec2, typename T1 = typename Vec1::value_type, typename T2 = typename Vec2::value_type>
class VectorPair {
public:
    VectorPair() = default;
    VectorPair(size_t size) : vec1(size), vec2(size) {}
    VectorPair(Vec1&& vec1, Vec2&& vec2) : vec1(vec1), vec2(vec2) {
        if (this->vec1.size() != this->vec2.size()) {
            throw std::runtime_error("VectorPair requires vectors of equal length");
        }
    }
    
    using value_type = std::pair<T1, T2>;
    using reference = VectorPairSetter<VectorPair<Vec1, Vec2, T1, T2>>;

    value_type at(size_t i) const {
        return value_type(vec1.at(i), vec2.at(i));
    }
    
    value_type operator[](size_t i) const {
        return at(i);
    }
    
    reference operator[](size_t i) {
        return reference(*this, i);
    }

    size_t size() const {
        return vec1.size();
    }

    bool empty() const {
        return vec1.empty();
    }
    
    size_t memory_size() const {
        return get_vec_memory_size(vec1) + get_vec_memory_size(vec2);
    }

private:
    
    friend reference;

    Vec1 vec1;
    Vec2 vec2;
};

}

#endif /* centrolign_vector_support_hpp */
