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


// modified from sdsl-lite by simon gog
// using built-in method or
// 64-bit version of 32-bit proposal of
// http://www-graphics.stanford.edu/~seander/bithacks.html
const uint32_t lt_hi[] = {
    0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7
};
uint32_t hi_bit(uint64_t x)
{
#ifdef __GNUG__
    if (x == 0)
        return 0;
    return 63 - __builtin_clzll(x);
#else
    uint64_t t,tt; // temporaries
    if ((tt = x >> 32)) { // hi >= 32
        if ((t = tt >> 16)) { // hi >= 48
            return (tt = t >> 8) ? 56 + lt_hi[tt] : 48 + lt_hi[t];
        } else { // hi < 48
            return (t = tt >> 8) ? 40 + lt_hi[t] : 32 + lt_hi[tt];
        }
    } else { // hi < 32
        if ((t = x >> 16)) { // hi >= 16
            return (tt = t >> 8) ? 24 + lt_hi[tt] : 16 + lt_hi[t];
        } else { // hi < 16
            return (tt = x >> 8) ?  8 + lt_hi[tt] : lt_hi[x];
        }
    }
#endif
}

}
