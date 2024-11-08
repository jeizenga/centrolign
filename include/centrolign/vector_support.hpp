#ifndef centrolign_vector_support_hpp
#define centrolign_vector_support_hpp

#include <cstdint>
#include <stdexcept>
#include <utility>

namespace centrolign {

/*
 * Functions as a reference for integer vector types for the purpose of operator[]
 */
template<class IntVectorType>
class IntVectorSetter {
public:
    
    IntVectorSetter(IntVectorType& vec, size_t i) : i(i), vec(&vec) {}
    IntVectorSetter() = default;
    ~IntVectorSetter() = default;\
    
    void operator++() {
        vec->set(i, vec->at(i) + 1);
    }
    
    void operator--() {
        vec->set(i, vec->at(i) - 1);
    }
    
    #define _assign(Type) inline void operator=(const Type& value) { vec->set(i, typename IntVectorType::value_type(value)); }
    _assign(uint8_t);
    _assign(uint16_t);
    _assign(uint32_t);
    _assign(uint64_t);
    _assign(int8_t);
    _assign(int16_t);
    _assign(int32_t);
    _assign(int64_t);
    _assign(float);
    _assign(double);
    _assign(long double);
    _assign(size_t);
    #undef _assign

    #define _convert(Type) inline operator Type() { return vec->at(i); }
    _convert(uint8_t);
    _convert(uint16_t);
    _convert(uint32_t);
    _convert(uint64_t);
    _convert(int8_t);
    _convert(int16_t);
    _convert(int32_t);
    _convert(int64_t);
    _convert(float);
    _convert(double);
    _convert(long double);
    _convert(size_t);
    #undef _convert
    
private:
    
    IntVectorType* vec = nullptr;
    size_t i = 0;
};

//template<class Vec1, class Vec2>
//class PairVector {
//public:
//    PairVector() = default;
//    PairVector(size_t size) : vec1(size), vec2(size) {}
//    PairVector(Vec1&& vec1, Vec2&& vec2) : vec1(vec1), vec2(vec2) {
//        if (this->vec1.size() != this->vec2.size()) {
//            throw std::runtime_error("PairVector requires vectors of equal length");
//        }
//    }
//
//    std::pair<typename Vec1::value_type, typename Vec2::value_type> at(size_t i) const {
//        return std::pair<typename Vec1::value_type, typename Vec2::value_type>(vec1.at(i), vec2.at());
//    }
//
//    std::pair<typename Vec1::value_type, typename Vec2::value_type> operator[](size_t i) const {
//        return std::pair<typename Vec1::value_type, typename Vec2::value_type>(vec1.at(i), vec2.at());
//    }
//
//    size_t size() const {
//        return vec1.size();
//    }
//
//    bool empty() const {
//        return vec1.empty();
//    }
//
//private:
//
//    Vec1 vec1;
//    Vec2 vec2;
//};

}

#endif /* centrolign_vector_support_hpp */
