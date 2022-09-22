#include "DubinsManeuver2d.hpp"

DubinsStruct DubinsManeuver2d::_lsl(DubinsManeuver2d maneuver, double a, double b, double d, double sa, double ca, double sb, double cb)
{
    double aux, t, p, q, length;
    string caseType;
    aux = atan2(cb - ca, d + sa - sb);
    t = Utility::mod2pi(-a + aux);
    p = sqrt(2 + d*d - 2*cos(a - b) + 2*d*(sa - sb));
    q = Utility::mod2pi(b - aux);
    length = (t + p + q) * maneuver._rhoMin;
    caseType = "LSL";
    return DubinsStruct {t, p, q, length, caseType};
}

DubinsStruct DubinsManeuver2d::_rsr(DubinsManeuver2d maneuver, double a, double b, double d, double sa, double ca, double sb, double cb)
{
    double aux, t, p, q, length;
    string caseType;
    aux = atan2(ca - cb, d - sa + sb);
    t = Utility::mod2pi(a - aux);
    p = sqrt(2 + d*d - 2*cos(a - b) + 2*d*(sb - sa));
    q = Utility::mod2pi(Utility::mod2pi(-b) + aux);
    length = (t + p + q) * maneuver._rhoMin;
    caseType = "RSR";
    return DubinsStruct {t, p, q, length, caseType};
}

DubinsStruct DubinsManeuver2d::_lsr(DubinsManeuver2d maneuver, double a, double b, double d, double sa, double ca, double sb, double cb)
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
    length = (t+p+q) * maneuver._rhoMin;
    caseType = "LSR";
    return DubinsStruct {t, p, q, length, caseType};
}

DubinsStruct DubinsManeuver2d::_rsl(DubinsManeuver2d maneuver, double a, double b, double d, double sa, double ca, double sb, double cb)
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
    length = (t+p+q) * maneuver._rhoMin;
    caseType = "RSL";
    return DubinsStruct {t, p, q, length, caseType};
}

DubinsStruct DubinsManeuver2d::_rlr(DubinsManeuver2d maneuver, double a, double b, double d, double sa, double ca, double sb, double cb)
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
    length = (t+p+q) * maneuver._rhoMin;
    caseType = "RLR";
    return DubinsStruct {t, p, q, length, caseType};
}

DubinsStruct DubinsManeuver2d::_lrl(DubinsManeuver2d maneuver, double a, double b, double d, double sa, double ca, double sb, double cb)
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
    length = (t+p+q) * maneuver._rhoMin;
    caseType = "LRL";
    return DubinsStruct {t, p, q, length, caseType};
}

DubinsStruct DubinsManeuver2d::_c(DubinsManeuver2d maneuver)
{
    double t, p , q, length;
    string caseType;
    t = 0;
    p = 2*M_PI;
    q = 0;
    length = (t+p+q) * maneuver._rhoMin;
    caseType = "RRR";
    return DubinsStruct {t, p, q, length, caseType};
}

DubinsManeuver2d DubinsManeuver2d::_buildBasicManeuver(vector<double> qi, vector<double> qf, double rhoMin)
{
    DubinsManeuver2d maneuver;
    maneuver._qi = qi;
    maneuver._qf = qf;
    maneuver._rhoMin = rhoMin;
    maneuver._maneuver = DubinsStruct {0, 0, 0, INFINITY, ""};
    return maneuver;
}

DubinsManeuver2d::DubinsManeuver2d() { }

DubinsManeuver2d DubinsManeuver2d::createDubinsManeuver2D(vector<double> qi, vector<double> qf, double rhoMin, double minLength, bool disableCCC)
{
    auto maneuver = _buildBasicManeuver(qi, qf, rhoMin);

    double dx = maneuver._qf[0] - maneuver._qi[0];
    double dy = maneuver._qf[1] - maneuver._qi[1];
    double D = sqrt(dx*dx + dy*dy);

    double d = D / maneuver._rhoMin;

    double rotationAngle = Utility::mod2pi(atan2(dy, dx));
    double a = Utility::mod2pi(maneuver._qi[2] - rotationAngle);
    double b = Utility::mod2pi(maneuver._qf[2] - rotationAngle);

    double sa = sin(a);
    double ca = cos(a);
    double sb = sin(b);
    double cb = cos(b);

    DubinsStruct pathRLR, pathLRL, pathC;
    DubinsStruct pathLSL = _lsl(maneuver, a, b, d, sa, ca, sb, cb);
    DubinsStruct pathRSR = _rsr(maneuver, a, b, d, sa, ca, sb, cb);
    DubinsStruct pathLSR = _lsr(maneuver, a, b, d, sa, ca, sb, cb);
    DubinsStruct pathRSL = _rsl(maneuver, a, b, d, sa, ca, sb, cb);

    vector<DubinsStruct> paths;
    if (disableCCC)
    {
        paths = { pathLSL, pathRSR, pathLSR, pathRSL };
    }
    else
    {
        pathRLR = _rlr(maneuver, a, b, d, sa, ca, sb, cb);
        pathLRL = _lrl(maneuver, a, b, d, sa, ca, sb, cb);
        paths = { pathLSL, pathRSR, pathLSR, pathRSL, pathRLR, pathLRL };
    }

    double rhoCompare = maneuver._rhoMin * 0.00001;
    if (abs(d) < rhoCompare && abs(a) < rhoCompare && abs(b) < rhoCompare)
    {
        double dist2d = max(abs(maneuver._qi[0] - maneuver._qf[0]), abs(maneuver._qi[1] - maneuver._qf[1]));
        if (dist2d < rhoCompare)
        {
            pathC = _c(maneuver);
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
        maneuver._maneuver = paths[0];
    }
    else
    {
        for (auto p : paths)
        {
            if (p.length >= minLength)
            {
                maneuver._maneuver = p;
                break;
            }
        }
    }

    return maneuver;
}

vector<double> DubinsManeuver2d::getCoordinatesAt(DubinsManeuver2d maneuver, double offset)
{
    double noffset = offset / maneuver._rhoMin;
    vector<double> qi = { 0, 0, maneuver._qi[2] };
    vector<double> q1, q2;

    double l1 = maneuver._maneuver.t;
    double l2 = maneuver._maneuver.p;

    vector<double> q;
    if (noffset < l1)
    {
        q = getPositionInSegment(maneuver, noffset, qi, maneuver._maneuver.caseType.at(0));
    }
    else if(noffset < (l1 + l2))
    {
        q1 = getPositionInSegment(maneuver, l1, qi, maneuver._maneuver.caseType.at(0));
        q = getPositionInSegment(maneuver, noffset - l1, q1, maneuver._maneuver.caseType.at(1));
    }
    else
    {
        q1 = getPositionInSegment(maneuver, l1, qi, maneuver._maneuver.caseType.at(0));
        q2 = getPositionInSegment(maneuver, l2, q1, maneuver._maneuver.caseType.at(1));
        q = getPositionInSegment(maneuver, noffset - l1 - l2, q2, maneuver._maneuver.caseType.at(2));
    }

    q[0] = q[0] * maneuver._rhoMin + maneuver._qi[0];
    q[1] = q[1] * maneuver._rhoMin + maneuver._qi[1];
    q[2] = Utility::mod2pi(q[2]);

    return q;
}

vector<double> DubinsManeuver2d::getPositionInSegment(DubinsManeuver2d maneuver, double offset, vector<double> qi, char caseType)
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

vector<vector<double>> DubinsManeuver2d::getSamplingPoints(DubinsManeuver2d maneuver, double res)
{
    int numPoints = (int)floor(maneuver._maneuver.length / res) + 1;
    vector<vector<double>> points;
    points.reserve(numPoints);

    for (int i = 0; i < numPoints; i++)
    {
        points.push_back(getCoordinatesAt(maneuver, i * res));
    }

    return points;
}