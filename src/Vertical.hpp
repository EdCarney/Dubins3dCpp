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
    static DubinsStruct _lsr(DubinsManeuver2d maneuver, tuple<double, double> pitchMax);
    static DubinsStruct _rsl(DubinsManeuver2d maneuver, tuple<double, double> pitchMax);

    public:
        static DubinsManeuver2d createDubinsManeuver2D(State2d qi, State2d qf, double rhoMin, tuple<double, double> pitchMax);
};

#endif //H_VERTICAL