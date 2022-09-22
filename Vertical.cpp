#include "Vertical.hpp"

DubinsStruct Vertical::_lsl(DubinsManeuver2d maneuver)
{
    double theta1 = maneuver._qi[2];
    double theta2 = maneuver._qf[2];

    double t = INFINITY;
    double p = INFINITY;
    double q = INFINITY;

    if (theta1 <= theta2)
    {
        double radius = maneuver._rhoMin;

        vector<double> p1 = { maneuver._qi.at(0), maneuver._qi.at(1) };
        vector<double> p2 = { maneuver._qf.at(0), maneuver._qf.at(1) };

        double c1 = radius * cos(theta1);
        double s1 = radius * sin(theta1);
        double c2 = radius * cos(theta2);
        double s2 = radius * sin(theta2);

        vector<double> o1 = { p1.at(0) - s1, p1.at(1) + c1 };
        vector<double> o2 = { p2.at(0) - s2, p2.at(1) + c2 };

        vector<double> diff = { o2.at(0) - o1.at(0), o2.at(1) - o1.at(1) };
        double centerDistance = hypot(diff.at(0), diff.at(1));
        double centerAngle = atan2(diff.at(1), diff.at(0));

        t = Utility::mod2pi(centerAngle - theta1);
        p = centerDistance / radius;
        q = Utility::mod2pi(theta2 - centerAngle);

        double turnEndY, diffY, tol = 0.00001;

        if (t > M_PI)
        {
            t = 0;
            q = theta2 - theta1;
            turnEndY = o2.at(1) - radius*cos(theta1);
            diffY = turnEndY - p1.at(1);
            if (abs(theta1) > tol && (diffY < 0 == theta1 < 0))
            {
                p = (diffY / sin(theta1)) / radius;
            }
            else
            {
                t = INFINITY;
                p = INFINITY;
                q = INFINITY;
            }
        }
        if (q > M_PI)
        {
            t = theta2 - theta1;
            q = 0;
            turnEndY = o1.at(1) - radius*cos(theta2);
            diffY = p2.at(1) - turnEndY;
            if (abs(theta2) > tol && (diffY < 0 == theta2 < 0))
            {
                p = (diffY / sin(theta2)) / radius;
            }
            else
            {
                t = INFINITY;
                p = INFINITY;
                q = INFINITY;
            }
        }
    }

    double length = (t + p + q) * maneuver._rhoMin;
    string caseType = "LSL";

    return DubinsStruct { t, p, q, length, caseType };
}

DubinsStruct Vertical::_rsr(DubinsManeuver2d maneuver)
{
    double theta1 = maneuver._qi[2];
    double theta2 = maneuver._qf[2];

    double t = INFINITY;
    double p = INFINITY;
    double q = INFINITY;

    if (theta1 <= theta2)
    {
        double radius = maneuver._rhoMin;

        vector<double> p1 = { maneuver._qi.at(0), maneuver._qi.at(1) };
        vector<double> p2 = { maneuver._qf.at(0), maneuver._qf.at(1) };

        double c1 = radius * cos(theta1);
        double s1 = radius * sin(theta1);
        double c2 = radius * cos(theta2);
        double s2 = radius * sin(theta2);

        vector<double> o1 = { p1.at(0) + s1, p1.at(1) - c1 };
        vector<double> o2 = { p2.at(0) + s2, p2.at(1) - c2 };

        vector<double> diff = { o2.at(0) - o1.at(0), o2.at(1) - o1.at(1) };
        double centerDistance = hypot(diff.at(0), diff.at(1));
        double centerAngle = atan2(diff.at(1), diff.at(0));

        t = Utility::mod2pi(theta1 - centerAngle);
        p = centerDistance / radius;
        q = Utility::mod2pi(centerAngle - theta2 );

        double turnEndY, diffY, tol = 0.00001;

        if (t > M_PI)
        {
            t = 0;
            q = theta1 - theta2;
            turnEndY = o2.at(1) + radius*cos(theta1);
            diffY = turnEndY - p1.at(1);
            if (abs(theta1) > tol && (diffY < 0 == theta1 < 0))
            {
                p = (diffY / sin(theta1)) / radius;
            }
            else
            {
                t = INFINITY;
                p = INFINITY;
                q = INFINITY;
            }
        }
        if (q > M_PI)
        {
            t = theta1 - theta2;
            q = 0;
            turnEndY = o1.at(1) + radius*cos(theta2);
            diffY = p2.at(1) - turnEndY;
            if (abs(theta2) > tol && (diffY < 0 == theta2 < 0))
            {
                p = (diffY / sin(theta2)) / radius;
            }
            else
            {
                t = INFINITY;
                p = INFINITY;
                q = INFINITY;
            }
        }
    }

    double length = (t + p + q) * maneuver._rhoMin;
    string caseType = "RSR";

    return DubinsStruct { t, p, q, length, caseType };
}

DubinsStruct Vertical::_lsr(DubinsManeuver2d maneuver, vector<double> pitchMax)
{
    double theta1 = maneuver._qi[2];
    double theta2 = maneuver._qf[2];

    vector<double> p1 = { maneuver._qi.at(0), maneuver._qi.at(1) };
    vector<double> p2 = { maneuver._qf.at(0), maneuver._qf.at(1) };

    double radius = maneuver._rhoMin;

    double c1 = radius * cos(theta1);
    double s1 = radius * sin(theta1);
    double c2 = radius * cos(theta2);
    double s2 = radius * sin(theta2);

    vector<double> o1 = { p1.at(0) - s1, p1.at(1) + c1 };
    vector<double> o2 = { p2.at(0) + s2, p2.at(1) - c2 };

    vector<double> diff = { o2.at(0) - o1.at(0), o2.at(1) - o1.at(1) };
    double centerDistance = hypot(diff.at(0), diff.at(1));

    double alpha;
    if (centerDistance < 2.0 * radius)
    {
        diff[0] = sqrt(4.0 * radius * radius - diff.at(1) * diff.at(1));
        alpha = M_PI / 2.0;
    }
    else
    {
        alpha = asin(2.0 * radius / centerDistance);
    }

    double centerAngle = atan2(diff.at(1), diff.at(0)) + alpha;

    double t, p, q;
    if (centerAngle < pitchMax.at(1))
    {
        t = Utility::mod2pi(centerAngle - theta1);
        p = sqrt(max(0.0, centerDistance * centerDistance - 4.0 * radius * radius)) / radius;
        q = Utility::mod2pi(centerAngle - theta2);
    }
    else
    {
        centerAngle = pitchMax.at(1);
        t = Utility::mod2pi(centerAngle - theta1);
        q = Utility::mod2pi(centerAngle - theta2);

        c1 = radius * cos(centerAngle);
        s1 = radius * sin(centerAngle);
        c2 = radius * cos(centerAngle);
        s2 = radius * sin(centerAngle);

        vector<double> w1 = { o1.at(0) + s1, o1.at(1) - c1 };
        vector<double> w2 = { o2.at(0) - s2, o2.at(1) + c2 };

        p = ((w2.at(1) - w1.at(1)) / sin(centerAngle)) / radius;
    }

    double length = (t + p + q) * maneuver._rhoMin;
    string caseType = "LSR";

    return DubinsStruct { t, p, q, length, caseType };
}

DubinsStruct Vertical::_rsl(DubinsManeuver2d maneuver, vector<double> pitchMax)
{
    double theta1 = maneuver._qi[2];
    double theta2 = maneuver._qf[2];

    vector<double> p1 = { maneuver._qi.at(0), maneuver._qi.at(1) };
    vector<double> p2 = { maneuver._qf.at(0), maneuver._qf.at(1) };

    double radius = maneuver._rhoMin;

    double c1 = radius * cos(theta1);
    double s1 = radius * sin(theta1);
    double c2 = radius * cos(theta2);
    double s2 = radius * sin(theta2);

    vector<double> o1 = { p1.at(0) + s1, p1.at(1) - c1 };
    vector<double> o2 = { p2.at(0) - s2, p2.at(1) + c2 };

    vector<double> diff = { o2.at(0) - o1.at(0), o2.at(1) - o1.at(1) };
    double centerDistance = hypot(diff.at(0), diff.at(1));

    double alpha;
    if (centerDistance < 2.0 * radius)
    {
        diff[0] = sqrt(4.0 * radius * radius - diff.at(1) * diff.at(1));
        alpha = M_PI / 2.0;
    }
    else
    {
        alpha = asin(2.0 * radius / centerDistance);
    }

    double centerAngle = atan2(diff.at(1), diff.at(0)) - alpha;

    double t, p, q;
    if (centerAngle > pitchMax.at(0))
    {
        t = Utility::mod2pi(theta1 - centerAngle);
        p = sqrt(max(0.0, centerDistance * centerDistance - 4.0 * radius * radius)) / radius;
        q = Utility::mod2pi(theta2 - centerAngle);
    }
    else
    {
        centerAngle = pitchMax.at(0);
        t = Utility::mod2pi(theta1 - centerAngle);
        q = Utility::mod2pi(theta2 - centerAngle);

        c1 = radius * cos(centerAngle);
        s1 = radius * sin(centerAngle);
        c2 = radius * cos(centerAngle);
        s2 = radius * sin(centerAngle);

        vector<double> w1 = { o1.at(0) - s1, o1.at(1) + c1 };
        vector<double> w2 = { o2.at(0) + s2, o2.at(1) - c2 };

        p = ((w2.at(1) - w1.at(1)) / sin(centerAngle)) / radius;
    }

    double length = (t + p + q) * maneuver._rhoMin;
    string caseType = "RSL";

    return DubinsStruct { t, p, q, length, caseType };
}

DubinsManeuver2d Vertical::createDubinsManeuver2D(vector<double> qi, vector<double> qf, double rhoMin, vector<double> pitchMax)
{
    auto maneuver = DubinsManeuver2d(qi, qf, rhoMin);

    double dx = maneuver._qf[0] - maneuver._qi[0];
    double dy = maneuver._qf[1] - maneuver._qi[1];
    double D = sqrt(dx*dx + dy*dy);

    double d = D / maneuver._rhoMin;

    double rotationAngle = Utility::mod2pi(atan2(dy, dx));
    double a = Utility::mod2pi(maneuver._qi[2] - rotationAngle);
    double b = Utility::mod2pi(maneuver._qf[2] - rotationAngle);

    auto pathLSL = _lsl(maneuver);
    auto pathRSR = _rsr(maneuver);
    auto pathLSR = _lsr(maneuver, pitchMax);
    auto pathRSL = _rsl(maneuver, pitchMax);

    vector<DubinsStruct> paths = { pathLSR, pathLSL, pathRSR, pathRSL };

    sort(paths.begin(), paths.end(),
        [](const DubinsStruct& a, const DubinsStruct& b) -> bool
        {
            return a.length < b.length;
        });

    for (auto p : paths)
    {
        if (abs(p.t) < M_PI && abs(p.q) < M_PI)
        {
            double centerAngle = maneuver._qi[2] + (p.caseType.at(0) == 'L' ? p.t : -p.t);
            if (centerAngle > pitchMax[0] && centerAngle < pitchMax[1])
            {
                maneuver._maneuver = p;
                break;
            }
        }
    }

    return maneuver;
}