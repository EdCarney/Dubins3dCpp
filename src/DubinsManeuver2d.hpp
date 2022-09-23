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
    double _rhoMin;
    vector<double> _qi;
    vector<double> _qf;
    DubinsStruct _maneuver;

    void _generateManeuver(double minLength, bool disableCCC);
    vector<double> _getPositionInSegment(double offset, vector<double> qi, char caseType) const;
    DubinsStruct _lsl(double a, double b, double d, double sa, double ca, double sb, double cb) const;
    DubinsStruct _rsr(double a, double b, double d, double sa, double ca, double sb, double cb) const;
    DubinsStruct _lsr(double a, double b, double d, double sa, double ca, double sb, double cb) const;
    DubinsStruct _rsl(double a, double b, double d, double sa, double ca, double sb, double cb) const;
    DubinsStruct _rlr(double a, double b, double d, double sa, double ca, double sb, double cb) const;
    DubinsStruct _lrl(double a, double b, double d, double sa, double ca, double sb, double cb) const;
    DubinsStruct _c() const;

    public:
        double rhoMin() const;
        const vector<double>& qi() const;
        const vector<double>& qf() const;
        double qi(int i) const;
        double qf(int i) const;
        const DubinsStruct& maneuver() const;
        void setManeuver(DubinsStruct maneuver);

        vector<double> getCoordinatesAt(double offset) const;
        vector<vector<double>> getSamplingPoints(double res = 0.1) const;
        DubinsManeuver2d(vector<double> qi, vector<double> qf, double rhoMin = 1, double minLength = -1, bool disableCCC = false);
};

#endif //DUBINS_MANEUVER_2D