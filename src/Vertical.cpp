#include "Vertical.hpp"

DubinsStruct Vertical::_lsl(DubinsManeuver2d maneuver)
{
    double theta1 = maneuver.qi().theta;
    double theta2 = maneuver.qf().theta;

    double t = INFINITY;
    double p = INFINITY;
    double q = INFINITY;

    if (theta1 <= theta2)
    {
        double radius = maneuver.rhoMin();

        Point2d p1 = { maneuver.qi().x, maneuver.qi().y };
        Point2d p2 = { maneuver.qf().x, maneuver.qf().y };

        double c1 = radius * cos(theta1);
        double s1 = radius * sin(theta1);
        double c2 = radius * cos(theta2);
        double s2 = radius * sin(theta2);

        Point2d o1 = { p1.x - s1, p1.y + c1 };
        Point2d o2 = { p2.x - s2, p2.y + c2 };
        Point2d diff = o2 - o1;
        double centerDistance = o2.distanceTo(o1);
        double centerAngle = atan2(diff.y, diff.x);

        t = Utility::mod2pi(centerAngle - theta1);
        p = centerDistance / radius;
        q = Utility::mod2pi(theta2 - centerAngle);

        double turnEndY, diffY, tol = 0.00001;

        if (t > M_PI)
        {
            t = 0;
            q = theta2 - theta1;
            turnEndY = o2.y - radius * cos(theta1);
            diffY = turnEndY - p1.y;
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
            turnEndY = o1.y - radius*cos(theta2);
            diffY = p2.y - turnEndY;
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

    double length = (t + p + q) * maneuver.rhoMin();
    string caseType = "LSL";

    return DubinsStruct { t, p, q, length, caseType };
}

DubinsStruct Vertical::_rsr(DubinsManeuver2d maneuver)
{
    double theta1 = maneuver.qi().theta;
    double theta2 = maneuver.qf().theta;

    double t = INFINITY;
    double p = INFINITY;
    double q = INFINITY;

    if (theta1 <= theta2)
    {
        double radius = maneuver.rhoMin();

        Point2d p1 = { maneuver.qi().x, maneuver.qi().y };
        Point2d p2 = { maneuver.qf().x, maneuver.qf().y };

        double c1 = radius * cos(theta1);
        double s1 = radius * sin(theta1);
        double c2 = radius * cos(theta2);
        double s2 = radius * sin(theta2);

        Point2d o1 = { p1.x + s1, p1.y - c1 };
        Point2d o2 = { p2.x + s2, p2.y - c2 };
        Point2d diff = o2 - o1;
        double centerDistance = o2.distanceTo(o1);
        double centerAngle = atan2(diff.y, diff.x);

        t = Utility::mod2pi(theta1 - centerAngle);
        p = centerDistance / radius;
        q = Utility::mod2pi(centerAngle - theta2 );

        double turnEndY, diffY, tol = 0.00001;

        if (t > M_PI)
        {
            t = 0;
            q = theta1 - theta2;
            turnEndY = o2.y + radius*cos(theta1);
            diffY = turnEndY - p1.y;
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
            turnEndY = o1.y + radius*cos(theta2);
            diffY = p2.y - turnEndY;
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

    double length = (t + p + q) * maneuver.rhoMin();
    string caseType = "RSR";

    return DubinsStruct { t, p, q, length, caseType };
}

DubinsStruct Vertical::_lsr(DubinsManeuver2d maneuver, tuple<double, double> pitchMax)
{
    double theta1 = maneuver.qi().theta;
    double theta2 = maneuver.qf().theta;

    Point2d p1 = { maneuver.qi().x, maneuver.qi().y };
    Point2d p2 = { maneuver.qf().x, maneuver.qf().y };

    double radius = maneuver.rhoMin();

    double c1 = radius * cos(theta1);
    double s1 = radius * sin(theta1);
    double c2 = radius * cos(theta2);
    double s2 = radius * sin(theta2);

    Point2d o1 = { p1.x - s1, p1.y + c1 };
    Point2d o2 = { p2.x + s2, p2.y - c2 };
    Point2d diff = o2 - o1;
    double centerDistance = o2.distanceTo(o1);

    double alpha;
    if (centerDistance < 2.0 * radius)
    {
        diff.x = sqrt(4.0 * radius * radius - diff.y * diff.y);
        alpha = M_PI / 2.0;
    }
    else
    {
        alpha = asin(2.0 * radius / centerDistance);
    }

    double centerAngle = atan2(diff.y, diff.x) + alpha;

    double t, p, q;
    if (centerAngle < get<1>(pitchMax))
    {
        t = Utility::mod2pi(centerAngle - theta1);
        p = sqrt(max(0.0, centerDistance * centerDistance - 4.0 * radius * radius)) / radius;
        q = Utility::mod2pi(centerAngle - theta2);
    }
    else
    {
        centerAngle = get<1>(pitchMax);
        t = Utility::mod2pi(centerAngle - theta1);
        q = Utility::mod2pi(centerAngle - theta2);

        c1 = radius * cos(centerAngle);
        s1 = radius * sin(centerAngle);
        c2 = radius * cos(centerAngle);
        s2 = radius * sin(centerAngle);

        Point2d w1 = { o1.x + s1, o1.y - c1 };
        Point2d w2 = { o2.x - s2, o2.y + c2 };

        p = ((w2.y - w1.y) / sin(centerAngle)) / radius;
    }

    double length = (t + p + q) * maneuver.rhoMin();
    string caseType = "LSR";

    return DubinsStruct { t, p, q, length, caseType };
}

DubinsStruct Vertical::_rsl(DubinsManeuver2d maneuver, tuple<double, double> pitchMax)
{
    double theta1 = maneuver.qi().theta;
    double theta2 = maneuver.qf().theta;

    Point2d p1 = { maneuver.qi().x, maneuver.qi().y };
    Point2d p2 = { maneuver.qf().x, maneuver.qf().y };

    double radius = maneuver.rhoMin();

    double c1 = radius * cos(theta1);
    double s1 = radius * sin(theta1);
    double c2 = radius * cos(theta2);
    double s2 = radius * sin(theta2);

    Point2d o1 = { p1.x + s1, p1.y - c1 };
    Point2d o2 = { p2.x - s2, p2.y + c2 };

    Point2d diff = o2 - o1;
    double centerDistance = o2.distanceTo(o1);

    double alpha;
    if (centerDistance < 2.0 * radius)
    {
        diff.x = sqrt(4.0 * radius * radius - diff.y * diff.y);
        alpha = M_PI / 2.0;
    }
    else
    {
        alpha = asin(2.0 * radius / centerDistance);
    }

    double centerAngle = atan2(diff.y, diff.x) - alpha;

    double t, p, q;
    if (centerAngle > get<0>(pitchMax))
    {
        t = Utility::mod2pi(theta1 - centerAngle);
        p = sqrt(max(0.0, centerDistance * centerDistance - 4.0 * radius * radius)) / radius;
        q = Utility::mod2pi(theta2 - centerAngle);
    }
    else
    {
        centerAngle = get<0>(pitchMax);
        t = Utility::mod2pi(theta1 - centerAngle);
        q = Utility::mod2pi(theta2 - centerAngle);

        c1 = radius * cos(centerAngle);
        s1 = radius * sin(centerAngle);
        c2 = radius * cos(centerAngle);
        s2 = radius * sin(centerAngle);

        Point2d w1 = { o1.x - s1, o1.y + c1 };
        Point2d w2 = { o2.x + s2, o2.y - c2 };

        p = ((w2.y - w1.y) / sin(centerAngle)) / radius;
    }

    double length = (t + p + q) * maneuver.rhoMin();
    string caseType = "RSL";

    return DubinsStruct { t, p, q, length, caseType };
}

DubinsManeuver2d Vertical::createDubinsManeuver2D(State2d qi, State2d qf, double rhoMin, tuple<double, double> pitchMax)
{
    auto maneuver = DubinsManeuver2d(qi, qf, rhoMin);

    Point2d diff = maneuver.qf() - maneuver.qi();
    double D = maneuver.qf().distanceTo(maneuver.qi());

    double d = D / maneuver.rhoMin();

    double rotationAngle = Utility::mod2pi(atan2(diff.y, diff.x));
    double a = Utility::mod2pi(maneuver.qi().theta - rotationAngle);
    double b = Utility::mod2pi(maneuver.qf().theta - rotationAngle);

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
            double centerAngle = maneuver.qi().theta + (p.caseType.at(0) == 'L' ? p.t : -p.t);
            if (centerAngle > get<0>(pitchMax) && centerAngle < get<1>(pitchMax))
            {
                maneuver.setManeuver(p);
                break;
            }
        }
    }

    return maneuver;
}