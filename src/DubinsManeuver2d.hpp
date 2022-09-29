#include <math.h>
#include <string>
#include <vector>
#include "Geometry.hpp"
#include "Utility.hpp"

using namespace std;

#ifndef DUBINS_MANEUVER_2D
#define DUBINS_MANEUVER_2D

struct DubinsStruct
{
    double t, p, q, length;
    string caseType;
};

struct DubinsParams
{
    double a, b, d, sa, ca, sb, cb;
};

class DubinsManeuver2d
{
    double _rhoMin;
    State2d _qi, _qf;
    DubinsStruct _maneuver;

    void _generateManeuver(double minLength, bool disableCCC);
    State2d _getPositionInSegment(double offset, State2d qi, char caseType) const;
    DubinsStruct _lsl(const DubinsParams& params) const;
    DubinsStruct _rsr(const DubinsParams& params) const;
    DubinsStruct _lsr(const DubinsParams& params) const;
    DubinsStruct _rsl(const DubinsParams& params) const;
    DubinsStruct _rlr(const DubinsParams& params) const;
    DubinsStruct _lrl(const DubinsParams& params) const;
    DubinsStruct _c() const;

    public:
        double rhoMin() const;
        const State2d& qi() const;
        const State2d& qf() const;
        const DubinsStruct& maneuver() const;
        void setManeuver(DubinsStruct maneuver);

        State2d getCoordinatesAt(double offset) const;
        vector<State2d> getSamplingPoints(double res = 0.1) const;
        DubinsManeuver2d();
        DubinsManeuver2d(const State2d& qi, const State2d& qf, double rhoMin = 1, double minLength = -1, bool disableCCC = false);
};

#endif //DUBINS_MANEUVER_2D