#include "utility.hpp"

namespace centrolign {

using namespace std;

vector<size_t> range_vector(size_t size) {
    vector<size_t> range(size, 0);
    for (size_t i = 1; i < range.size(); ++i) {
        range[i] = i;
    }
    return range;
}

}
