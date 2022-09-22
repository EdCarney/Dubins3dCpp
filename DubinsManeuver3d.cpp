#include "DubinsManeuver3d.hpp"

vector<vector<double>> DubinsManeuver3d::_computeSampling(DubinsManeuver3d maneuver, int numSamples)
{
    auto Dlat = maneuver._path.at(0);
    auto Dlon = maneuver._path.at(1);

    vector<vector<double>> points(numSamples);
    double offsetLon, sampleLen = Dlon._maneuver.length / numSamples;

    for (int i = 0; i < numSamples; ++i)
    {
        offsetLon = sampleLen * i;
        auto qSZ = DubinsManeuver2d::getCoordinatesAt(Dlon, offsetLon);
        auto qXY = DubinsManeuver2d::getCoordinatesAt(Dlat, qSZ.at(0));
        points.push_back(vector<double> { qXY.at(0), qXY.at(1), qSZ.at(1), qXY.at(2), qSZ.at(2) });
    }

    return points;
}

vector<DubinsManeuver2d> DubinsManeuver3d::_tryToConstruct(DubinsManeuver3d maneuver, double horizontalRadius)
{
    vector<double> qi2d = { maneuver._qi.at(0), maneuver._qi.at(1), maneuver._qi.at(3) };
    vector<double> qf2d = { maneuver._qf.at(0), maneuver._qf.at(1), maneuver._qf.at(3) };

    auto Dlat = DubinsManeuver2d::createDubinsManeuver2d(qi2d, qf2d, horizontalRadius);

    vector<double> qi3d = { 0.0, maneuver._qi.at(2), maneuver._qi.at(4) };
    vector<double> qf3d = { Dlat._maneuver.length, maneuver._qf.at(2), maneuver._qf.at(4) };

    double tol = 0.00001;
    double verticalCurvature = sqrt((1.0 / pow(maneuver._rhoMin, 2)) - (1.0 / pow(horizontalRadius, 2)));

    if (verticalCurvature < tol)
        return vector<DubinsManeuver2d> { };

    double verticalRadius = 1.0 / verticalCurvature;

    auto Dlon = DubinsManeuver2d::createDubinsManeuver2d(qi3d, qf3d, verticalRadius);

    // TODO: error here in conditional? supposed to be Dlan in second case?
    if (Dlon._maneuver.caseType == "RLR" || Dlon._maneuver.caseType == "RLR")
        return vector<DubinsManeuver2d> { };

    if (Dlon._maneuver.caseType.at(0) == 'R')
    {
        if (maneuver._qi.at(4) - Dlon._maneuver.t < maneuver._pitchLims.at(0))
        {
            return vector<DubinsManeuver2d> { };
        }
    }
    else
    {
        if (maneuver._qi.at(4) + Dlon._maneuver.t > maneuver._pitchLims.at(1))
        {
            return vector<DubinsManeuver2d> { };
        }
    }

    return vector<DubinsManeuver2d> { Dlat, Dlon };
}

DubinsManeuver3d DubinsManeuver3d::_getLowerBound(vector<double> qi, vector<double> qf, double rhoMin, vector<double> pitchLims)
{
    auto maneuver = DubinsManeuver3d(qi, qf, rhoMin, pitchLims);

    double spiralRadius = rhoMin * (pow(cos(max(-pitchLims.at(0), pitchLims.at(1))), 2));

    vector<double> qi2d = { maneuver._qi.at(0), maneuver._qi.at(1), maneuver._qi.at(3) };
    vector<double> qf2d = { maneuver._qf.at(0), maneuver._qf.at(1), maneuver._qf.at(3) };
    auto Dlat = DubinsManeuver2d::createDubinsManeuver2d(qi2d, qf2d, spiralRadius);

    vector<double> qi3d = { 0.0, maneuver._qi.at(2), maneuver._qi.at(4) };
    vector<double> qf3d = { Dlat._maneuver.length, maneuver._qf.at(2), maneuver._qf.at(4) };
    auto Dlon = Vertical::createDubinsManeuver2D(qi3d, qf3d, maneuver._rhoMin, maneuver._pitchLims);

    if (Dlon._maneuver.caseType == "XXX")
    {
        maneuver._length = 0.0;
        return maneuver;
    }

    maneuver._path = { Dlat, Dlon };
    maneuver._length = Dlon._maneuver.length;
    return maneuver;
}

DubinsManeuver3d DubinsManeuver3d::_getUpperBound(vector<double> qi, vector<double> qf, double rhoMin, vector<double> pitchLims)
{
    auto maneuver = DubinsManeuver3d(qi, qf, rhoMin, pitchLims);

    double safeRadius = sqrt(2.0) * maneuver._rhoMin;

    vector<double> p1 = { qi.at(0), qi.at(1) };
    vector<double> p2 = { qf.at(0), qf.at(1) };
    vector<double> diff = { p2.at(0) - p1.at(0), p2.at(1) - p1.at(1) };
    double dist = hypot(diff.at(0), diff.at(1));

    if (dist < 4.0 * safeRadius)
    {
        maneuver._length = INFINITY;
        return maneuver;
    }

    vector<double> qi2d = { maneuver._qi.at(0), maneuver._qi.at(1), maneuver._qi.at(3) };
    vector<double> qf2d = { maneuver._qf.at(0), maneuver._qf.at(1), maneuver._qf.at(3) };
    auto Dlat = DubinsManeuver2d::createDubinsManeuver2d(qi2d, qf2d, safeRadius);

    vector<double> qi3d = { 0.0, maneuver._qi.at(2), maneuver._qi.at(4) };
    vector<double> qf3d = { Dlat._maneuver.length, maneuver._qf.at(2), maneuver._qf.at(4) };
    auto Dlon = Vertical::createDubinsManeuver2D(qi3d, qf3d, maneuver._rhoMin, maneuver._pitchLims);

    if (Dlon._maneuver.caseType == "XXX")
    {
        maneuver._length = INFINITY;
        return maneuver;
    }

    maneuver._path = { Dlat, Dlon };
    maneuver._length = Dlon._maneuver.length;
    return maneuver;
}

DubinsManeuver3d::DubinsManeuver3d(vector<double> qi, vector<double> qf, double rhoMin, vector<double> pitchLims)
{
    _qi = qi;
    _qf = qf;
    _rhoMin = rhoMin;
    _pitchLims = pitchLims;
    _length = -1;
}

DubinsManeuver3d DubinsManeuver3d::createDubinsManeuver3d(vector<double> qi, vector<double> qf, double rhoMin, vector<double> pitchLims)
{
    DubinsManeuver3d maneuver(qi, qf, rhoMin, pitchLims);

    double zi = maneuver._qi.at(2);
    double zf = maneuver._qf.at(2);
    double dz = zf - zi;

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
            if (fc.at(1)._maneuver.length < fb.at(1)._maneuver.length)
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
    maneuver._length = fb.at(1)._maneuver.length;
    return maneuver;
}