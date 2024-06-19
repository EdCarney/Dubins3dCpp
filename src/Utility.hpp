#ifndef H_UTILITY
#define H_UTILITY

#include <cmath>
#include <numbers>

struct Utility {
    static double mod2pi(const double val) {
        constexpr double modVal = 2 * std::numbers::pi;
        const int num = static_cast<int>(std::floor(val / modVal));
        return val - num * modVal;
    }
};

#endif //H_UTILITY