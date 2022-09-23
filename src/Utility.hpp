#include <math.h>
#include <vector>

using namespace std;

#ifndef H_UTILITY
#define H_UTILITY

struct Utility
{
    static double mod2pi(double val)
    {
        double modVal = 2 * M_PI;
        int num = (int)floor(val / modVal);
        return val - num * modVal;
    }
};

#endif //H_UTILITY