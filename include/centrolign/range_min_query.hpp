#ifndef centrolign_range_min_query_hpp
#define centrolign_range_min_query_hpp

#include <vector>
#include <cmath>
#include <iostream>
#include <bitset>

#include "centrolign/utility.hpp"

namespace centrolign {

/*
 * takes O(n) space and computes RMQ in O(1) time
 */
template<typename T>
class RMQ {
public:
    
    // note: data structure becomes invalid if the array is destroyed
    // or modified
    RMQ(const std::vector<T>& arr);
    RMQ() = default;
    ~RMQ() = default;
    
    RMQ(RMQ<T>&& other) = default;
    RMQ<T>& operator=(RMQ<T>&& other) = default;
    
    // return the index of the minimum in the range [begin, end)
    size_t range_arg_min(size_t begin, size_t end) const;
    
protected:
    
    /*
     * table of all possible RMQ queries
     */
    class ExhaustiveRMQ {
    public:
        
        ExhaustiveRMQ() noexcept = default;
        ~ExhaustiveRMQ() noexcept = default;
        
        void initialize(typename std::vector<T>::const_iterator begin,
                        typename std::vector<T>::const_iterator end);
        
        // only valid if initialized
        size_t range_arg_min(size_t begin, size_t end) const;
        
        bool initialized() const;
    private:
        
        // indexed by [end][begin]
        std::vector<std::vector<size_t>> table;
    };
    
    /*
     * sparse table data structure to compute O(1) RMQ with O(n log n) space
     */
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
        
    };
    
    // original array
    const std::vector<T>* arr = nullptr;
    // length of the blocks we divide into
    size_t block_size = 0;
    
    // the min value of each block
    std::vector<T> block_min;
    // the location of the min value of each block
    std::vector<size_t> block_arg_min;
    // an RMQ over the block mins
    SparseTable block_sparse_table;
    
    // the Cartesian tree bit code of each block
    std::vector<uint16_t> block_tree_code;
    // exhaustively enumerated tables for each Cartesian tree
    // indexed by bit code
    std::vector<ExhaustiveRMQ> tree_tables;
    // the last one is a special case because it might not have the same
    // size as the others
    ExhaustiveRMQ final_table;
    
    /*
     * binary tree that completely specifies RMQ results
     */
    class CartesianTree {
    public:
        // the range must be non-empty to be valid
        CartesianTree(typename std::vector<T>::const_iterator begin,
                      typename std::vector<T>::const_iterator end);
        CartesianTree() = default;
        ~CartesianTree() = default;
        
        // return an integer encoding of the tree's topology
        uint16_t bit_code() const;
        
    private:
        struct Node {
            
            Node() = default;
            ~Node() = default;
            
            size_t left = -1;
            size_t right = -1;
            size_t parent = -1;
        };
        
        size_t root;
        std::vector<Node> nodes;
    };
};

/*
 * Template implementations
 */

static bool debug_rmq = false;

template<typename T>
RMQ<T>::RMQ(const std::vector<T>& arr) : block_size(ceil(log2(arr.size()) / 4.0)), arr(&arr) {
    
    
    // initialize block-wise records
    size_t num_blocks = (arr.size() + block_size - 1) / block_size;
    block_arg_min.resize(num_blocks);
    block_min.resize(num_blocks);
    block_tree_code.resize(num_blocks - 1); // the final one is outside the vector
    
    if (debug_rmq) {
        std::cerr << "constructing RMQ for " << arr.size() << " elements with " << num_blocks << " blocks of size " << block_size << '\n';
    }
    
    // tree bit codes are upperbounded by 2^2n
    tree_tables.resize(1 << (2 * block_size));
    
    // process each block
    for (size_t i = 0; i < num_blocks; ++i) {
        
        size_t arr_begin = i * block_size;
        size_t arr_end = std::min(arr_begin + block_size, arr.size());
        
        // TODO: we can skip this part on either the first or the last, since they aren't
        // queried but it might not actually save much time
        
        // record min of block
        size_t arg_min = arr_begin;
        for (size_t j = arr_begin + 1; j < arr_end; ++j) {
            if (arr[j] < arr[arg_min]) {
                arg_min = j;
            }
        }
        block_min[i] = arr[arg_min];
        block_arg_min[i] = arg_min;
        
        if (i + 1 < num_blocks) {
            // compute cartesian tree and initialize tree structures
            CartesianTree cartesian_tree(arr.begin() + arr_begin, arr.begin() + arr_end);
            uint16_t tree_code = cartesian_tree.bit_code();
            block_tree_code[i] = tree_code;
            if (!tree_tables[tree_code].initialized()) {
                if (debug_rmq) {
                    std::cerr << "initializing table for tree " << tree_code << " with interval " << arr_begin << ":" << arr_end << '\n';
                }
                tree_tables[tree_code].initialize(arr.begin() + arr_begin, arr.begin() + arr_end);
            }
        }
        else {
            // the last one might not be the same size, so the tree codes are inconsistent
            // and we have to handle it separately
            if (debug_rmq) {
                std::cerr << "initializing final table for interval " << arr_begin << ":" << arr_end << '\n';
            }
            final_table.initialize(arr.begin() + arr_begin, arr.begin() + arr_end);
        }
    }
    
    if (debug_rmq) {
        std::cerr << "block mins:" << '\n';
        for (size_t i = 0; i < block_min.size(); ++i) {
            if (i) {
                std::cerr << ' ';
            }
            std::cerr << block_min[i];
        }
        std::cerr << '\n';
        std::cerr << "block arg mins:" << '\n';
        for (size_t i = 0; i < block_arg_min.size(); ++i) {
            if (i) {
                std::cerr << ' ';
            }
            std::cerr << block_arg_min[i];
        }
        std::cerr << '\n';
        std::cerr << "tree codes:" << '\n';
        for (size_t i = 0; i < block_tree_code.size(); ++i) {
            if (i) {
                std::cerr << ' ';
            }
            std::cerr << std::bitset<8>(block_tree_code[i]);
        }
        std::cerr << '\n';
    }
    
    // build a sparse table over block-wise minimums
    block_sparse_table = SparseTable(block_min);
}

template<typename T>
size_t RMQ<T>::range_arg_min(size_t begin, size_t end) const {
    
    size_t first_block = begin / block_size;
    size_t final_block = (end - 1) / block_size;
    
    size_t arg_min;
    if (first_block == final_block) {
        // look up in one tree
        const auto& table = final_block == block_tree_code.size() ? final_table : tree_tables[block_tree_code[first_block]];
        size_t block_begin = first_block * block_size;
        arg_min = block_begin + table.range_arg_min(begin - block_begin, end - block_begin);
    }
    else {
        // look up in two trees
        const auto& table1 = tree_tables[block_tree_code[first_block]];
        size_t block_begin1 = first_block * block_size;
        arg_min = block_begin1 + table1.range_arg_min(begin - block_begin1, block_size);
        
        const auto& table2 = final_block == block_tree_code.size() ? final_table : tree_tables[block_tree_code[final_block]];
        size_t block_begin2 = final_block * block_size;
        size_t arg_min2 = block_begin2 + table2.range_arg_min(0, end - block_begin2);
        if (final_block != first_block + 1) {
            // also look up in the sparse table
            size_t min_block = block_sparse_table.range_arg_min(first_block + 1, final_block);
            size_t arg_min3 = block_arg_min[min_block];
            if ((*arr)[arg_min3] < (*arr)[arg_min]) {
                arg_min = arg_min3;
            }
        }
        // check min here to ensure that we always choose the leftmost
        if ((*arr)[arg_min2] < (*arr)[arg_min]) {
            arg_min = arg_min2;
        }
        
    }
    
    return arg_min;
}

template<typename T>
RMQ<T>::CartesianTree::CartesianTree(typename std::vector<T>::const_iterator begin,
                                     typename std::vector<T>::const_iterator end) {
    
    // initialize the nodes
    nodes.resize(end - begin);
    
    // we sweep over the nodes left-to-right, extending the prefix that is a valid
    // cartesian tree in each iteration
    root = 0;
    for (size_t i = 1; i < nodes.size(); ++i) {
        // find the highest node to the left that has a lower value
        size_t here = i - 1;
        while (here != -1 && *(begin + i) < *(begin + here)) {
            here = nodes[here].parent;
        }
        
        if (here != -1) {
            // the node here is the parent, and the new node takes the place of its right child
            nodes[i].left = nodes[here].right;
            if (nodes[i].left != -1) {
                nodes[nodes[i].left].parent = i;
            }
            nodes[here].right = i;
            nodes[i].parent = here;
        }
        else {
            // this entry is the min so far, it becomes the root
            nodes[i].left = root;
            nodes[root].parent = i;
            root = i;
        }
    }
}

template<typename T>
uint16_t RMQ<T>::CartesianTree::bit_code() const {
    uint16_t code = 0;
    uint32_t pos = 0;
    
    std::vector<size_t> stack(1, root);
    while (!stack.empty()) {
        size_t here = stack.back();
        stack.pop_back();
        
        if (nodes[here].left != -1) {
            code |= (1 << pos);
            stack.push_back(nodes[here].left);
        }
        ++pos;
        if (nodes[here].right != -1) {
            code |= (1 << pos);
            stack.push_back(nodes[here].right);
        }
        ++pos;
    }
    return code;
}


template<typename T>
void RMQ<T>::ExhaustiveRMQ::initialize(typename std::vector<T>::const_iterator begin,
                                       typename std::vector<T>::const_iterator end) {
    
    size_t size = end - begin;
    
    table.resize(size + 1);
    for (size_t i = 1; i < table.size(); ++i) {
        table[i].resize(i);
    }
    
    // seed the iteration with the trivial length 1 queries
    for (size_t i = 0; i < size; ++i) {
        table[i + 1][i] = i;
    }
    
    // fill out the rest of the table with queries of increasing length
    for (size_t len = 2; len <= size; ++len) {
        for (size_t b = 0, n = size - len; b <= n; ++b) {
            size_t e = b + len;
            size_t left_arg_min = table[e - 1][b]; // all but the last entry
            size_t right_arg_min = e - 1; // the last entry
            table[e][b] = (*(begin + left_arg_min) <= *(begin + right_arg_min)) ? left_arg_min : right_arg_min;
        }
    }
    if (debug_rmq) {
        std::cerr << "exhaustive table for vector of size " << size << ":\n";
        for (size_t i = 0; i < size; ++i) {
            for (size_t j = 0; j < table.size(); ++j) {
                if (j) {
                    std::cerr << ' ';
                }
                if (j <= i) {
                    std::cerr << '.';
                }
                else {
                    std::cerr << table[j][i];
                }
            }
            std::cerr << '\n';
        }
    }
}

template<typename T>
size_t RMQ<T>::ExhaustiveRMQ::range_arg_min(size_t begin, size_t end) const {
    return table[end][begin];
}

template<typename T>
bool RMQ<T>::ExhaustiveRMQ::initialized() const {
    return !table.empty();
}

template<typename T>
RMQ<T>::SparseTable::SparseTable(const std::vector<T>& arr) : arr(&arr) {
    
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
    
    if (debug_rmq) {
        std::cerr << "sparse table:\n";
        for (size_t k = 0; k < table.size(); ++k) {
            std::cerr << (1 << k) << ":\t";
            const auto& row = table[k];
            for (size_t i = 0; i < row.size(); ++i) {
                if (i) {
                    std::cerr << ' ';
                }
                std::cerr << row[i];
            }
            std::cerr << '\n';
        }
    }
}

template<typename T>
size_t RMQ<T>::SparseTable::range_arg_min(size_t begin, size_t end) const  {
    
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

#endif /* centrolign_range_min_query_hpp */
