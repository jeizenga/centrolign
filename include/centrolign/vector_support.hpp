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
