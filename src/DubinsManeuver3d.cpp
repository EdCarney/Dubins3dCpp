#include <vector>
#include "DubinsManeuver3d.hpp"

void DubinsManeuver3d::_generateManeuver() {
    constexpr double a = 1.0;
    double b = 1.0;

    const auto fa = _tryToConstruct(_rhoMin * a);
    auto fb = _tryToConstruct(_rhoMin * b);

    while (fb.size() < 2) {
        b *= 2.0;
        fb = _tryToConstruct(_rhoMin * b);
    }

    if (!fa.empty()) {
        _path = fa;
    }
    else {
        if (fb.size() < 2) {
            // ERROR NO MANEUVER!!
            return;
        }
    }

    double step = 0.1;
    constexpr double tol = 0.00001;
    while (std::abs(step) > tol) {
        double c = b + step;
        c = c < 1.0 ? 1.0 : c;

        if (const auto fc = _tryToConstruct(_rhoMin * c); !fc.empty()) {
            if (fc.at(1).maneuver().length < fb.at(1).maneuver().length) {
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

std::vector<DubinsManeuver2d> DubinsManeuver3d::_tryToConstruct(const double horizontalRadius) const {
    const State2d qi2d = { _qi.x, _qi.y, _qi.theta };
    const State2d qf2d = { _qf.x, _qf.y, _qf.theta };

    const auto Dlat = DubinsManeuver2d(qi2d, qf2d, horizontalRadius);

    const State2d qi3d = { 0.0, _qi.z, _qi.gamma };
    const State2d qf3d = { Dlat.maneuver().length, _qf.z, _qf.gamma };

    const double verticalCurvature = sqrt((1.0 / pow(_rhoMin, 2)) - (1.0 / pow(horizontalRadius, 2)));

    if (constexpr double tol = 0.00001; verticalCurvature < tol)
        return std::vector<DubinsManeuver2d> { };

    const double verticalRadius = 1.0 / verticalCurvature;

    const auto Dlon = DubinsManeuver2d(qi3d, qf3d, verticalRadius);

    // TODO: error here in conditional? supposed to be Dlan in second case?
    if (Dlon.maneuver().caseType == "RLR" || Dlon.maneuver().caseType == "RLR")
        return std::vector<DubinsManeuver2d> { };

    if (Dlon.maneuver().caseType.at(0) == 'R') {
        if (_qi.gamma - Dlon.maneuver().t < minPitch()) {
            return std::vector<DubinsManeuver2d> { };
        }
    }
    else {
        if (_qi.gamma + Dlon.maneuver().t > maxPitch()) {
            return std::vector<DubinsManeuver2d> { };
        }
    }

    return std::vector { Dlat, Dlon };
}

std::vector<State3d> DubinsManeuver3d::computeSampling(const int numSamples) const {
    const auto Dlat = _path.at(0);
    const auto Dlon = _path.at(1);

    std::vector<State3d> points(numSamples);
    const double sampleLen = Dlon.maneuver().length / numSamples;

    for (int i = 0; i < numSamples; ++i)
    {
        const double offsetLon = sampleLen * i;
        const auto qSZ = Dlon.getCoordinatesAt(offsetLon);
        const auto qXY = Dlat.getCoordinatesAt(qSZ.x);
        points[i] = State3d { qXY.x, qXY.y, qSZ.y, qXY.theta, qSZ.theta };
    }

    return points;
}