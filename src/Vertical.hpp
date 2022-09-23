#include <math.h>
#include <vector>
#include "DubinsManeuver2d.hpp"
#include "Utility.hpp"

using namespace std;

#ifndef H_VERTICAL
#define H_VERTICAL

class Vertical
{
    static DubinsStruct _lsl(DubinsManeuver2d maneuver);
    static DubinsStruct _rsr(DubinsManeuver2d maneuver);
    static DubinsStruct _lsr(DubinsManeuver2d maneuver, vector<double> pitchMax);
    static DubinsStruct _rsl(DubinsManeuver2d maneuver, vector<double> pitchMax);

    public:
        static DubinsManeuver2d createDubinsManeuver2D(vector<double> qi, vector<double> qf, double rhoMin, vector<double> pitchMax);
};

#endif //H_VERTICAL