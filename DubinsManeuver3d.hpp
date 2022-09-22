#include <math.h>
#include <string>
#include <vector>
#include "Utility.hpp"
#include "Vertical.hpp"

using namespace std;

#ifndef DUBINS_MANEUVER_3D
#define DUBINS_MANEUVER_3D

class DubinsManeuver3d
{
    static vector<vector<double>> _computeSampling(DubinsManeuver3d maneuver, int numSamples = 1000);
    static vector<DubinsManeuver2d> _tryToConstruct(DubinsManeuver3d maneuver, double horizontalRadius);
    static DubinsManeuver3d _getLowerBound(vector<double> qi, vector<double> qf, double rhoMin = 1, vector<double> pitchLims = { -M_PI / 4.0, M_PI / 2.0 });
    static DubinsManeuver3d _getUpperBound(vector<double> qi, vector<double> qf, double rhoMin = 1, vector<double> pitchLims = { -M_PI / 4.0, M_PI / 2.0 });

    public:
        vector<double> _qi;
        vector<double> _qf;

        double _rhoMin;
        vector<double> _pitchLims;

        vector<DubinsManeuver2d> _path;
        double _length;

        DubinsManeuver3d(vector<double> qi, vector<double> qf, double rhoMin, vector<double> pitchLims);
        static DubinsManeuver3d createDubinsManeuver3d(vector<double> qi, vector<double> qf, double rhoMin, vector<double> pitchLims);
};

#endif //DUBINS_MANEUVER_3D