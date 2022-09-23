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

class DubinsManeuver2d
{
    double _rhoMin;
    State2d _qi, _qf;
    DubinsStruct _maneuver;

    void _generateManeuver(double minLength, bool disableCCC);
    State2d _getPositionInSegment(double offset, State2d qi, char caseType) const;
    DubinsStruct _lsl(double a, double b, double d, double sa, double ca, double sb, double cb) const;
    DubinsStruct _rsr(double a, double b, double d, double sa, double ca, double sb, double cb) const;
    DubinsStruct _lsr(double a, double b, double d, double sa, double ca, double sb, double cb) const;
    DubinsStruct _rsl(double a, double b, double d, double sa, double ca, double sb, double cb) const;
    DubinsStruct _rlr(double a, double b, double d, double sa, double ca, double sb, double cb) const;
    DubinsStruct _lrl(double a, double b, double d, double sa, double ca, double sb, double cb) const;
    DubinsStruct _c() const;

    public:
        double rhoMin() const;
        const State2d& qi() const;
        const State2d& qf() const;
        const DubinsStruct& maneuver() const;
        void setManeuver(DubinsStruct maneuver);

        State2d getCoordinatesAt(double offset) const;
        vector<State2d> getSamplingPoints(double res = 0.1) const;
        DubinsManeuver2d(State2d qi, State2d qf, double rhoMin = 1, double minLength = -1, bool disableCCC = false);
};

#endif //DUBINS_MANEUVER_2D