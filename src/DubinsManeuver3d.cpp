#include "DubinsManeuver3d.hpp"

void DubinsManeuver3d::_generateManeuver()
{
    double dz = _qf.z - _qi.z;

    double a = 1.0;
    double b = 1.0;

    vector<DubinsManeuver2d> fa, fb, fc;
    fa = _tryToConstruct(_rhoMin * a);
    fb = _tryToConstruct(_rhoMin * b);

    while (fb.size() < 2)
    {
        b *= 2.0;
        fb = _tryToConstruct(_rhoMin * b);
    }

    if (fa.size() > 0)
    {
        _path = fa;
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

        fc = _tryToConstruct(_rhoMin * c);
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

    _path = fb;
    _length = fb.at(1).maneuver().length;
}

vector<DubinsManeuver2d> DubinsManeuver3d::_tryToConstruct(double horizontalRadius)
{
    State2d qi2d = { _qi.x, _qi.y, _qi.theta };
    State2d qf2d = { _qf.x, _qf.y, _qf.theta };

    auto Dlat = DubinsManeuver2d(qi2d, qf2d, horizontalRadius);

    State2d qi3d = { 0.0, _qi.z, _qi.gamma };
    State2d qf3d = { Dlat.maneuver().length, _qf.z, _qf.gamma };

    double tol = 0.00001;
    double verticalCurvature = sqrt((1.0 / pow(_rhoMin, 2)) - (1.0 / pow(horizontalRadius, 2)));

    if (verticalCurvature < tol)
        return vector<DubinsManeuver2d> { };

    double verticalRadius = 1.0 / verticalCurvature;

    auto Dlon = DubinsManeuver2d(qi3d, qf3d, verticalRadius);

    // TODO: error here in conditional? supposed to be Dlan in second case?
    if (Dlon.maneuver().caseType == "RLR" || Dlon.maneuver().caseType == "RLR")
        return vector<DubinsManeuver2d> { };

    if (Dlon.maneuver().caseType.at(0) == 'R')
    {
        if (_qi.gamma - Dlon.maneuver().t < _pitchLims.at(0))
        {
            return vector<DubinsManeuver2d> { };
        }
    }
    else
    {
        if (_qi.gamma + Dlon.maneuver().t > _pitchLims.at(1))
        {
            return vector<DubinsManeuver2d> { };
        }
    }

    return vector<DubinsManeuver2d> { Dlat, Dlon };
}

DubinsManeuver3d DubinsManeuver3d::_getLowerBound(double rhoMin, vector<double> pitchLims)
{
    auto maneuver = DubinsManeuver3d(_qi, _qf, rhoMin, pitchLims);

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

DubinsManeuver3d DubinsManeuver3d::_getUpperBound(double rhoMin, vector<double> pitchLims)
{
    auto maneuver = DubinsManeuver3d(_qi, _qf, rhoMin, pitchLims);

    double safeRadius = sqrt(2.0) * maneuver._rhoMin;

    Point2d p1 = { _qi.x, _qi.y };
    Point2d p2 = { _qf.x, _qf.y };
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

double DubinsManeuver3d::rhoMin() const { return _rhoMin; }

double DubinsManeuver3d::length() const { return _length; }

const State3d& DubinsManeuver3d::qi() const { return _qi; }

const State3d& DubinsManeuver3d::qf() const { return _qf; }

double DubinsManeuver3d::minPitch() const { return _pitchLims.at(0); }

double DubinsManeuver3d::maxPitch() const { return _pitchLims.at(1); }

const vector<DubinsManeuver2d>& DubinsManeuver3d::path() const { return _path; }

void DubinsManeuver3d::setPath(vector<DubinsManeuver2d> path) { _path = path; }

DubinsManeuver3d::DubinsManeuver3d(State3d qi, State3d qf, double rhoMin, vector<double> pitchLims)
{
    _qi = qi;
    _qf = qf;
    _rhoMin = rhoMin;
    _pitchLims = pitchLims;
    _length = -1;
    _generateManeuver();
}

vector<State3d> DubinsManeuver3d::computeSampling(int numSamples)
{
    auto Dlat = _path.at(0);
    auto Dlon = _path.at(1);

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