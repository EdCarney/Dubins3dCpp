#include "DubinsManeuver2d.hpp"

void DubinsManeuver2d::_generateManeuver(double minLength, bool disableCCC)
{
    double dx = qf(0) - qi(0);
    double dy = qf(1) - qi(1);
    double D = sqrt(dx*dx + dy*dy);

    double d = D / rhoMin();

    double rotationAngle = Utility::mod2pi(atan2(dy, dx));
    double a = Utility::mod2pi(qi(2) - rotationAngle);
    double b = Utility::mod2pi(qf(2) - rotationAngle);

    double sa = sin(a);
    double ca = cos(a);
    double sb = sin(b);
    double cb = cos(b);

    DubinsStruct pathRLR, pathLRL, pathC;
    DubinsStruct pathLSL = _lsl(a, b, d, sa, ca, sb, cb);
    DubinsStruct pathRSR = _rsr(a, b, d, sa, ca, sb, cb);
    DubinsStruct pathLSR = _lsr(a, b, d, sa, ca, sb, cb);
    DubinsStruct pathRSL = _rsl(a, b, d, sa, ca, sb, cb);

    vector<DubinsStruct> paths;
    if (disableCCC)
    {
        paths = { pathLSL, pathRSR, pathLSR, pathRSL };
    }
    else
    {
        pathRLR = _rlr(a, b, d, sa, ca, sb, cb);
        pathLRL = _lrl(a, b, d, sa, ca, sb, cb);
        paths = { pathLSL, pathRSR, pathLSR, pathRSL, pathRLR, pathLRL };
    }

    double rhoCompare = rhoMin() * 0.00001;
    if (abs(d) < rhoCompare && abs(a) < rhoCompare && abs(b) < rhoCompare)
    {
        double dist2d = max(abs(qi(0) - qf(0)), abs(qi(1) - qf(1)));
        if (dist2d < rhoCompare)
        {
            pathC = _c();
            paths = { pathC };
        }
    }

    sort(paths.begin(), paths.end(),
        [](const DubinsStruct& a, const DubinsStruct& b) -> bool
        {
            return a.length < b.length;
        });

    if (minLength < 0)
    {
        _maneuver = paths[0];
    }
    else
    {
        for (auto p : paths)
        {
            if (p.length >= minLength)
            {
                _maneuver = p;
                break;
            }
        }
    }
}

vector<double> DubinsManeuver2d::_getPositionInSegment(double offset, vector<double> qi, char caseType) const
{
    vector<double> q = { 0, 0, 0 };
    switch (caseType)
    {
        case 'L':
            q[0] = qi[0] + sin(qi[2]+offset) - sin(qi[2]);
            q[1] = qi[1] - cos(qi[2]+offset) + cos(qi[2]);
            q[2] = qi[2] + offset;
            break;
        case 'R':
            q[0] = qi[0] - sin(qi[2]-offset) + sin(qi[2]);
            q[1] = qi[1] + cos(qi[2]-offset) - cos(qi[2]);
            q[2] = qi[2] - offset;
            break;
        case 'S':
            q[0] = qi[0] + cos(qi[2]) * offset;
            q[1] = qi[1] + sin(qi[2]) * offset;
            q[2] = qi[2];
            break;
    }
    return q;
}

DubinsStruct DubinsManeuver2d::_lsl(double a, double b, double d, double sa, double ca, double sb, double cb) const
{
    double aux, t, p, q, length;
    string caseType;
    aux = atan2(cb - ca, d + sa - sb);
    t = Utility::mod2pi(-a + aux);
    p = sqrt(2 + d*d - 2*cos(a - b) + 2*d*(sa - sb));
    q = Utility::mod2pi(b - aux);
    length = (t + p + q) * rhoMin();
    caseType = "LSL";
    return DubinsStruct {t, p, q, length, caseType};
}

DubinsStruct DubinsManeuver2d::_rsr(double a, double b, double d, double sa, double ca, double sb, double cb) const
{
    double aux, t, p, q, length;
    string caseType;
    aux = atan2(ca - cb, d - sa + sb);
    t = Utility::mod2pi(a - aux);
    p = sqrt(2 + d*d - 2*cos(a - b) + 2*d*(sb - sa));
    q = Utility::mod2pi(Utility::mod2pi(-b) + aux);
    length = (t + p + q) * rhoMin();
    caseType = "RSR";
    return DubinsStruct {t, p, q, length, caseType};
}

DubinsStruct DubinsManeuver2d::_lsr(double a, double b, double d, double sa, double ca, double sb, double cb) const
{
    double aux1, aux2, t, p, q, length;
    string caseType;
    aux1 = -2 + d*d + 2*cos(a - b) + 2*d*(sa + sb);
    if (aux1 > 0)
    {
        p = sqrt(aux1);
        aux2 = atan2(-ca-cb, d+sa+sb) - atan(-2/p);
        t = Utility::mod2pi(-a + aux2);
        q = Utility::mod2pi(-Utility::mod2pi(b) + aux2);
    }
    else
    {
        t = INFINITY;
        p = INFINITY;
        q = INFINITY;
    }
    length = (t+p+q) * rhoMin();
    caseType = "LSR";
    return DubinsStruct {t, p, q, length, caseType};
}

DubinsStruct DubinsManeuver2d::_rsl(double a, double b, double d, double sa, double ca, double sb, double cb) const
{
    double aux1, aux2, t, p, q, length;
    string caseType;
    aux1 = d*d - 2 + 2*cos(a-b) - 2*d*(sa+sb);
    if (aux1 > 0)
    {
        p = sqrt(aux1);
        aux2 = atan2(ca+cb, d-sa-sb) - atan(2/p);
        t = Utility::mod2pi(a - aux2);
        q = Utility::mod2pi(Utility::mod2pi(b) - aux2);
    }
    else
    {
        t = INFINITY;
        p = INFINITY;
        q = INFINITY;
    }
    length = (t+p+q) * rhoMin();
    caseType = "RSL";
    return DubinsStruct {t, p, q, length, caseType};
}

DubinsStruct DubinsManeuver2d::_rlr(double a, double b, double d, double sa, double ca, double sb, double cb) const
{
    double aux, t, p, q, length;
    string caseType;
    aux = (6 - d*d + 2*cos(a-b) + 2*d*(sa-sb))/8;
    if (abs(aux) <= 1)
    {
        p = Utility::mod2pi(-acos(aux));
        t = Utility::mod2pi(a - atan2(ca-cb, d-sa+sb) + p/2);
        q = Utility::mod2pi(a - b - t + p);
    }
    else
    {
        t = INFINITY;
        p = INFINITY;
        q = INFINITY;
    }
    length = (t+p+q) * rhoMin();
    caseType = "RLR";
    return DubinsStruct {t, p, q, length, caseType};
}

DubinsStruct DubinsManeuver2d::_lrl(double a, double b, double d, double sa, double ca, double sb, double cb) const
{
    double aux, t, p, q, length;
    string caseType;
    aux = (6 - d*d + 2*cos(a-b) + 2*d*(-sa+sb))/8;;
    if (abs(aux) <= 1)
    {
        p = Utility::mod2pi(-acos(aux));
        t = Utility::mod2pi(-a + atan2(-ca+cb, d+sa-sb) + p/2);
        q = Utility::mod2pi(b - a - t + p);
    }
    else
    {
        t = INFINITY;
        p = INFINITY;
        q = INFINITY;
    }
    length = (t+p+q) * rhoMin();
    caseType = "LRL";
    return DubinsStruct {t, p, q, length, caseType};
}

DubinsStruct DubinsManeuver2d::_c() const
{
    double t, p , q, length;
    string caseType;
    t = 0;
    p = 2*M_PI;
    q = 0;
    length = (t+p+q) * rhoMin();
    caseType = "RRR";
    return DubinsStruct {t, p, q, length, caseType};
}

double DubinsManeuver2d::rhoMin() const { return _rhoMin; }

const vector<double>& DubinsManeuver2d::qi() const { return _qi; }

const vector<double>& DubinsManeuver2d::qf() const { return _qf; }

double DubinsManeuver2d::qi(int i) const { return _qi.at(i); }

double DubinsManeuver2d::qf(int i) const { return _qf.at(i); }

const DubinsStruct& DubinsManeuver2d::maneuver() const { return _maneuver; }

void DubinsManeuver2d::setManeuver(DubinsStruct maneuver) { _maneuver = maneuver; }

DubinsManeuver2d::DubinsManeuver2d(vector<double> qi, vector<double> qf, double rhoMin, double minLength, bool disableCCC)
{
    _qi = qi;
    _qf = qf;
    _rhoMin = rhoMin;
    _maneuver = DubinsStruct {0, 0, 0, INFINITY, ""};
    _generateManeuver(minLength, disableCCC);
}

vector<double> DubinsManeuver2d::getCoordinatesAt(double offset) const
{
    double noffset = offset / rhoMin();
    vector<double> qir = { 0, 0, qi(2) };
    vector<double> q1, q2;

    double l1 = maneuver().t;
    double l2 = maneuver().p;

    vector<double> q;
    if (noffset < l1)
    {
        q = _getPositionInSegment(noffset, qir, maneuver().caseType.at(0));
    }
    else if(noffset < (l1 + l2))
    {
        q1 = _getPositionInSegment(l1, qir, maneuver().caseType.at(0));
        q = _getPositionInSegment(noffset - l1, q1, maneuver().caseType.at(1));
    }
    else
    {
        q1 = _getPositionInSegment(l1, qir, maneuver().caseType.at(0));
        q2 = _getPositionInSegment(l2, q1, maneuver().caseType.at(1));
        q = _getPositionInSegment(noffset - l1 - l2, q2, maneuver().caseType.at(2));
    }

    q[0] = q[0] * rhoMin() + qi(0);
    q[1] = q[1] * rhoMin() + qi(1);
    q[2] = Utility::mod2pi(q[2]);

    return q;
}

vector<vector<double>> DubinsManeuver2d::getSamplingPoints(double res) const
{
    int numPoints = (int)floor(maneuver().length / res) + 1;
    vector<vector<double>> points;
    points.reserve(numPoints);

    for (int i = 0; i < numPoints; i++)
    {
        points.push_back(getCoordinatesAt(i * res));
    }

    return points;
}