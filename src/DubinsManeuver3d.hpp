#include <math.h>
#include <string>
#include <vector>
#include "Geometry.hpp"
#include "Utility.hpp"
#include "DubinsManeuver2d.hpp"

using namespace std;

#ifndef DUBINS_MANEUVER_3D
#define DUBINS_MANEUVER_3D

struct DubinsPath
{
    DubinsManeuver2d lat, lon;
    bool isEmpty() const { return &lat != NULL && &lon != NULL; }

    DubinsPath()
    {
        lat = *((DubinsManeuver2d*)NULL);
        lon = *((DubinsManeuver2d*)NULL);
    }
};

class DubinsManeuver3d
{
    double _rhoMin, _length;
    State3d _qi, _qf;
    tuple<double, double> _pitchLims;
    vector<DubinsManeuver2d> _path;

    void _generateManeuver();
    vector<DubinsManeuver2d> _tryToConstruct(double horizontalRadius) const;

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
        DubinsManeuver3d();
        DubinsManeuver3d(const State3d& qi, const State3d& qf, double rhoMin, const tuple<double, double>& pitchLims);
};

#endif //DUBINS_MANEUVER_3D