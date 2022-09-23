#include <math.h>
#include <string>
#include <vector>
#include "Geometry.hpp"
#include "Utility.hpp"
#include "Vertical.hpp"

using namespace std;

#ifndef DUBINS_MANEUVER_3D
#define DUBINS_MANEUVER_3D

class DubinsManeuver3d
{
    static vector<DubinsManeuver2d> _tryToConstruct(DubinsManeuver3d maneuver, double horizontalRadius);
    static DubinsManeuver3d _getLowerBound(State3d qi, State3d qf, double rhoMin = 1, vector<double> pitchLims = { -M_PI / 4.0, M_PI / 2.0 });
    static DubinsManeuver3d _getUpperBound(State3d qi, State3d qf, double rhoMin = 1, vector<double> pitchLims = { -M_PI / 4.0, M_PI / 2.0 });

    public:
        State3d _qi;
        State3d _qf;

        double _rhoMin;
        vector<double> _pitchLims;

        vector<DubinsManeuver2d> _path;
        double _length;

        DubinsManeuver3d(State3d qi, State3d qf, double rhoMin, vector<double> pitchLims);
        static DubinsManeuver3d createDubinsManeuver3d(State3d qi, State3d qf, double rhoMin, vector<double> pitchLims);
        static vector<State3d> computeSampling(DubinsManeuver3d maneuver, int numSamples = 1000);
};

#endif //DUBINS_MANEUVER_3D