#include "DubinsManeuver2d.hpp"

void DubinsManeuver2d::_generateManeuver(double minLength, bool disableCCC)
{
    double dx = qf().x - qi().x;
    double dy = qf().y - qi().y;
    double D = sqrt(dx*dx + dy*dy);

    double d = D / rhoMin();

    double rotationAngle = Utility::mod2pi(atan2(dy, dx));
    double a = Utility::mod2pi(qi().theta - rotationAngle);
    double b = Utility::mod2pi(qf().theta - rotationAngle);

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
        double dist2d = max(abs(qi().x - qf().x), abs(qi().y - qf().y));
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

State2d DubinsManeuver2d::_getPositionInSegment(double offset, State2d qi, char caseType) const
{
    State2d q = { 0, 0, 0 };
    switch (caseType)
    {
        case 'L':
            q.x = qi.x + sin(qi.theta + offset) - sin(qi.theta);
            q.y = qi.y - cos(qi.theta + offset) + cos(qi.theta);
            q.theta = qi.theta + offset;
            break;
        case 'R':
            q.x = qi.x - sin(qi.theta - offset) + sin(qi.theta);
            q.y = qi.y + cos(qi.theta - offset) - cos(qi.theta);
            q.theta = qi.theta - offset;
            break;
        case 'S':
            q.x = qi.x + cos(qi.theta) * offset;
            q.y = qi.y + sin(qi.theta) * offset;
            q.theta = qi.theta;
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

const State2d& DubinsManeuver2d::qi() const { return _qi; }

const State2d& DubinsManeuver2d::qf() const { return _qf; }

const DubinsStruct& DubinsManeuver2d::maneuver() const { return _maneuver; }

void DubinsManeuver2d::setManeuver(DubinsStruct maneuver) { _maneuver = maneuver; }

DubinsManeuver2d::DubinsManeuver2d(State2d qi, State2d qf, double rhoMin, double minLength, bool disableCCC)
{
    _qi = qi;
    _qf = qf;
    _rhoMin = rhoMin;
    _maneuver = DubinsStruct {0, 0, 0, INFINITY, ""};
    _generateManeuver(minLength, disableCCC);
}

State2d DubinsManeuver2d::getCoordinatesAt(double offset) const
{
    double noffset = offset / rhoMin();
    State2d qir = { 0, 0, qi().theta };
    State2d q, q1, q2;

    double l1 = maneuver().t;
    double l2 = maneuver().p;

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

    q.x = q.x * rhoMin() + qi().x;
    q.y = q.y * rhoMin() + qi().y;
    q.theta = Utility::mod2pi(q.theta);

    return q;
}

vector<State2d> DubinsManeuver2d::getSamplingPoints(double res) const
{
    int numPoints = (int)floor(maneuver().length / res) + 1;
    vector<State2d> points;
    points.reserve(numPoints);

    for (int i = 0; i < numPoints; i++)
    {
        points.push_back(getCoordinatesAt(i * res));
    }

    return points;
}