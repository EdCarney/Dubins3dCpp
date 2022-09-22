#include <math.h>

struct Utility
{
    static double mod2pi(double val)
    {
        double modVal = 2 * M_PI;
        int num = (int)floor(val / modVal);
        return val - num * modVal;
    }
};