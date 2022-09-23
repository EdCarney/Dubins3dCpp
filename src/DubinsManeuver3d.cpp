#include "DubinsManeuver3d.hpp"

vector<DubinsManeuver2d> DubinsManeuver3d::_tryToConstruct(DubinsManeuver3d maneuver, double horizontalRadius)
{
    State2d qi2d = { maneuver._qi.x, maneuver._qi.y, maneuver._qi.theta };
    State2d qf2d = { maneuver._qf.x, maneuver._qf.y, maneuver._qf.theta };

    auto Dlat = DubinsManeuver2d(qi2d, qf2d, horizontalRadius);

    State2d qi3d = { 0.0, maneuver._qi.z, maneuver._qi.gamma };
    State2d qf3d = { Dlat.maneuver().length, maneuver._qf.z, maneuver._qf.gamma };

    double tol = 0.00001;
    double verticalCurvature = sqrt((1.0 / pow(maneuver._rhoMin, 2)) - (1.0 / pow(horizontalRadius, 2)));

    if (verticalCurvature < tol)
        return vector<DubinsManeuver2d> { };

    double verticalRadius = 1.0 / verticalCurvature;

    auto Dlon = DubinsManeuver2d(qi3d, qf3d, verticalRadius);

    // TODO: error here in conditional? supposed to be Dlan in second case?
    if (Dlon.maneuver().caseType == "RLR" || Dlon.maneuver().caseType == "RLR")
        return vector<DubinsManeuver2d> { };

    if (Dlon.maneuver().caseType.at(0) == 'R')
    {
        if (maneuver._qi.gamma - Dlon.maneuver().t < maneuver._pitchLims.at(0))
        {
            return vector<DubinsManeuver2d> { };
        }
    }
    else
    {
        if (maneuver._qi.gamma + Dlon.maneuver().t > maneuver._pitchLims.at(1))
        {
            return vector<DubinsManeuver2d> { };
        }
    }

    return vector<DubinsManeuver2d> { Dlat, Dlon };
}

DubinsManeuver3d DubinsManeuver3d::_getLowerBound(State3d qi, State3d qf, double rhoMin, vector<double> pitchLims)
{
    auto maneuver = DubinsManeuver3d(qi, qf, rhoMin, pitchLims);

    double spiralRadius = rhoMin * (pow(cos(max(-pitchLims.at(0), pitchLims.at(1))), 2));

    State2d qi2d = { maneuver._qi.x, maneuver._qi.y, maneuver._qi.theta };
    State2d qf2d = { maneuver._qf.x, maneuver._qf.y, maneuver._qf.theta };
    auto Dlat = DubinsManeuver2d(qi2d, qf2d, spiralRadius);

    State2d qi3d = { 0.0, maneuver._qi.z, maneuver._qi.gamma };
    State2d qf3d = { Dlat.maneuver().length, maneuver._qf.z, maneuver._qf.gamma };
    auto Dlon = Vertical::createDubinsManeuver2D(qi3d, qf3d, maneuver._rhoMin, maneuver._pitchLims);

    if (Dlon.maneuver().caseType == "XXX")
    {
        maneuver._length = 0.0;
        return maneuver;
    }

    maneuver._path = { Dlat, Dlon };
    maneuver._length = Dlon.maneuver().length;
    return maneuver;
}

DubinsManeuver3d DubinsManeuver3d::_getUpperBound(State3d qi, State3d qf, double rhoMin, vector<double> pitchLims)
{
    auto maneuver = DubinsManeuver3d(qi, qf, rhoMin, pitchLims);

    double safeRadius = sqrt(2.0) * maneuver._rhoMin;

    Point2d p1 = { qi.x, qi.y };
    Point2d p2 = { qf.x, qf.y };
    Point2d diff = p2 - p1;
    double dist = p2.distanceTo(p1);

    if (dist < 4.0 * safeRadius)
    {
        maneuver._length = INFINITY;
        return maneuver;
    }

    State2d qi2d = { maneuver._qi.x, maneuver._qi.y, maneuver._qi.theta };
    State2d qf2d = { maneuver._qf.x, maneuver._qf.y, maneuver._qf.theta };
    auto Dlat = DubinsManeuver2d(qi2d, qf2d, safeRadius);

    State2d qi3d = { 0.0, maneuver._qi.z, maneuver._qi.gamma };
    State2d qf3d = { Dlat.maneuver().length, maneuver._qf.z, maneuver._qf.gamma };
    auto Dlon = Vertical::createDubinsManeuver2D(qi3d, qf3d, maneuver._rhoMin, maneuver._pitchLims);

    if (Dlon.maneuver().caseType == "XXX")
    {
        maneuver._length = INFINITY;
        return maneuver;
    }

    maneuver._path = { Dlat, Dlon };
    maneuver._length = Dlon.maneuver().length;
    return maneuver;
}

DubinsManeuver3d::DubinsManeuver3d(State3d qi, State3d qf, double rhoMin, vector<double> pitchLims)
{
    _qi = qi;
    _qf = qf;
    _rhoMin = rhoMin;
    _pitchLims = pitchLims;
    _length = -1;
}

DubinsManeuver3d DubinsManeuver3d::createDubinsManeuver3d(State3d qi, State3d qf, double rhoMin, vector<double> pitchLims)
{
    DubinsManeuver3d maneuver(qi, qf, rhoMin, pitchLims);

    double dz = qf.z - qi.z;

    double a = 1.0;
    double b = 1.0;

    vector<DubinsManeuver2d> fa, fb, fc;
    fa = _tryToConstruct(maneuver, maneuver._rhoMin * a);
    fb = _tryToConstruct(maneuver, maneuver._rhoMin * b);

    while (fb.size() < 2)
    {
        b *= 2.0;
        fb = _tryToConstruct(maneuver, maneuver._rhoMin * b);
    }

    if (fa.size() > 0)
    {
        maneuver._path = fa;
    }
    else
    {
        if (fb.size() < 2)
        {
            // ERROR NO MANEUVER!!
        }
    }

    double c, step = 0.1, tol = 0.00001;
    while (abs(step) > tol)
    {
        c = b + step;
        c = c < 1.0 ? 1.0 : c;

        fc = _tryToConstruct(maneuver, maneuver._rhoMin * c);
        if (fc.size() > 0)
        {
            if (fc.at(1).maneuver().length < fb.at(1).maneuver().length)
            {
                b = c;
                fb = fc;
                step *= 2.0;
                continue;
            }
        }
        step *= -0.1;
    }

    maneuver._path = fb;
    maneuver._length = fb.at(1).maneuver().length;
    return maneuver;
}

vector<State3d> DubinsManeuver3d::computeSampling(DubinsManeuver3d maneuver, int numSamples)
{
    auto Dlat = maneuver._path.at(0);
    auto Dlon = maneuver._path.at(1);

    vector<State3d> points(numSamples);
    double offsetLon, sampleLen = Dlon.maneuver().length / numSamples;

    for (int i = 0; i < numSamples; ++i)
    {
        offsetLon = sampleLen * i;
        auto qSZ = Dlon.getCoordinatesAt(offsetLon);
        auto qXY = Dlat.getCoordinatesAt(qSZ.x);
        points[i] = State3d { qXY.x, qXY.y, qSZ.y, qXY.theta, qSZ.theta };
    }

    return points;
}