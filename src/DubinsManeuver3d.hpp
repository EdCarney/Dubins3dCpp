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
    double _rhoMin, _length;
    State3d _qi, _qf;
    tuple<double, double> _pitchLims;
    vector<DubinsManeuver2d> _path;

    void _generateManeuver();
    vector<DubinsManeuver2d> _tryToConstruct(double horizontalRadius) const;
    DubinsManeuver3d _getLowerBound(double rhoMin = 1, tuple<double, double> pitchLims = { -M_PI / 4.0, M_PI / 2.0 }) const;
    DubinsManeuver3d _getUpperBound(double rhoMin = 1, tuple<double, double> pitchLims = { -M_PI / 4.0, M_PI / 2.0 }) const;

    public:
        double rhoMin() const;
        double length() const;
        const State3d& qi() const;
        const State3d& qf() const;
        double minPitch() const;
        double maxPitch() const;
        const vector<DubinsManeuver2d>& path() const;
        void setPath(vector<DubinsManeuver2d> path);

        vector<State3d> computeSampling(int numSamples = 1000) const;
        DubinsManeuver3d(State3d qi, State3d qf, double rhoMin, tuple<double, double> pitchLims);
};

#endif //DUBINS_MANEUVER_3D