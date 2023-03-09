#ifndef centrolign_sparse_table_hpp
#define centrolign_sparse_table_hpp

#include <vector>

#include "utility"

namespace centrolign {

/*
 * sparse table data structure to compute O(1) RMQ with O(n log n) space
 */
template<typename T>
class SparseTable {
public:
    
    SparseTable(const std::vector<T>& arr);
    SparseTable() = default;
    ~SparseTable() = default;
    
    // assumes but does not check that begin < end
    size_t range_arg_min(size_t begin, size_t end) const;
    
private:
    
    // TODO: if I built in indirection to the original array i wouldn't
    // need to copy this in the main RMQ data structure
    const std::vector<T>* arr = nullptr;
    // gives RMQ(i, i+2^k) at index [k][i]
    std::vector<std::vector<size_t>> table;
    
}

/*
 * Template implementations
 */

template<typename T>
SparseTable<T>::SparseTable(const std::vector<T>& arr) : arr(&arr) {
    
    // the length 1 stride vector is trivial
    table.emplace_back(std::move(range_vector(arr.size())));
    
    size_t stride = 2;
    for (size_t stride = 2; stride <= arr.size(); stride *= 2) {
        // add a new row to the table
        table.emplace_back(arr.size() - stride + 1);
        auto& curr = table.back();
        auto& prev = table[table.size() - 2];
        
        size_t half_stride = stride / 2;
        for (size_t i = 0; i < curr.size(); ++i) {
            // the new entry is composed of 2 previous entries
            size_t left_arg_min = prev[i];
            size_t right_arg_min = prev[i + half_stride];
            curr[i] = arr[left_arg_min] <= arr[right_arg_min] ? left_arg_min : right_arg_min;
        }
    }
}

template<typename T>
size_t SparseTable<T>::range_arg_min(size_t begin, size_t end) const  {
    
    size_t width = end - begin;
    size_t log_width = hi_bit(width);
    
    size_t query_size = 1 << log_width;
    if (width == query_size) {
        // the range corresponds to exactly one block
        return table[log_width][begin];
    }
    else {
        // we compose 2 overlapping blocks
        size_t left_arg_min = table[log_width][begin];
        size_t right_arg_min = table[log_width][end - query_size];
        return (*arr)[left_arg_min] <= (*arr)[right_arg_min] ? left_arg_min : right_arg_min;
    }
    
}

}

#endif /* centrolign_sparse_table_hpp */
