#ifndef centrolign_score_function_hpp
#define centrolign_score_function_hpp

#include <cstdint>
#include <cmath>
#include <stdexcept>
#include <string>
#include <limits>

namespace centrolign {

/*
 * A shared anchor score function used by other modules
 */
class ScoreFunction {
public:
    
    ScoreFunction() = default;
    ~ScoreFunction() = default;
    
    enum AnchorScore {InverseCount = 0, LengthScaleInverseCount = 1, ConcaveLengthScaleInverseCount = 2, ConcaveLengthScaleCountDifference = 3};
    
    // weight of a full anchor
    inline double anchor_weight(size_t count1, size_t count2, size_t length) const;
    // weight of a partial anchor
    inline double anchor_weight(size_t count1, size_t count2, size_t length, size_t full_length) const;
    
    // how to score anchors
    AnchorScore anchor_score_function = ConcaveLengthScaleInverseCount;
    
    // if using inverse count weighting, the power to raise the count to in the denominator
    double pair_count_power = 0.5;
    // if using concave length scale, the maximum positive-scoring length
    double length_intercept = 1750.0;
    // if using concave length scale, the decay behavior of the concave term (higher delays decay until longer lengths)
    double length_decay_power = 3.0;
    
    // records the intrinsic scale of the scoring function on these sequences
    double score_scale = 0.303092; // ~ chr12 value
};


/**
 * Inline implementations
 */

inline double ScoreFunction::anchor_weight(size_t count1, size_t count2, size_t length) const {
    return anchor_weight(count1, count2, length, length);
}

inline double ScoreFunction::anchor_weight(size_t count1, size_t count2, size_t length, size_t full_length) const {
    
    assert(length <= full_length);
    double count = count1 * count2;
    double fraction = double(length) / double(full_length);
    
    switch (anchor_score_function) {
        case InverseCount:
            return fraction / pow(count, pair_count_power);
            
        case LengthScaleInverseCount:
            return fraction * length / pow(count, pair_count_power);
            
        case ConcaveLengthScaleInverseCount:
            return fraction * (length / pow(count, pair_count_power) - pow(length / length_intercept, length_decay_power) * length_intercept);
            
        case ConcaveLengthScaleCountDifference:
            return fraction * (length - count * pow(length / length_intercept, length_decay_power) * length_intercept);
            
        default:
            throw std::runtime_error("Unrecognized anchor scoring function " + std::to_string((int) anchor_score_function));
            break;
    }
    return -std::numeric_limits<double>::infinity();
}


}

#endif /* centrolign_score_function_hpp */
