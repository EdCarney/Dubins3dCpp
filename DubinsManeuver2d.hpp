#include <math.h>
#include <string>
#include <vector>
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
    static DubinsStruct _lsl(DubinsManeuver2d maneuver, double a, double b, double d, double sa, double ca, double sb, double cb);
    static DubinsStruct _rsr(DubinsManeuver2d maneuver, double a, double b, double d, double sa, double ca, double sb, double cb);
    static DubinsStruct _lsr(DubinsManeuver2d maneuver, double a, double b, double d, double sa, double ca, double sb, double cb);
    static DubinsStruct _rsl(DubinsManeuver2d maneuver, double a, double b, double d, double sa, double ca, double sb, double cb);
    static DubinsStruct _rlr(DubinsManeuver2d maneuver, double a, double b, double d, double sa, double ca, double sb, double cb);
    static DubinsStruct _lrl(DubinsManeuver2d maneuver, double a, double b, double d, double sa, double ca, double sb, double cb);
    static DubinsStruct _c(DubinsManeuver2d maneuver);

    public:
        vector<double> _qi;
        vector<double> _qf;
        double _rhoMin;
        DubinsStruct _maneuver;

        DubinsManeuver2d(vector<double> qi, vector<double> qf, double rhoMin);
        static DubinsManeuver2d createDubinsManeuver2D(vector<double> qi, vector<double> qf, double rhoMin = 1, double minLength = -1, bool disableCCC = false);
        static vector<double> getCoordinatesAt(DubinsManeuver2d maneuver, double offset);
        static vector<double> getPositionInSegment(DubinsManeuver2d maneuver, double offset, vector<double> qi, char caseType);
        static vector<vector<double>> getSamplingPoints(DubinsManeuver2d maneuver, double res = 0.1);
};

#endif //DUBINS_MANEUVER_2D