#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>
#include "DubinsManeuver2d.hpp"
#include "Utility.hpp"

void DubinsManeuver2d::_generateManeuver(double minLength, bool disableCCC) {
    namespace ranges = std::ranges;

    const auto [x, y] = qf() - qi();
    const double D = qf().distanceTo(qi());

    const double d = D / rhoMin();

    const double rotationAngle = Utility::mod2pi(atan2(y, x));
    const double a = Utility::mod2pi(qi().theta - rotationAngle);
    const double b = Utility::mod2pi(qf().theta - rotationAngle);

    const double sa = std::sin(a);
    const double ca = std::cos(a);
    const double sb = std::sin(b);
    const double cb = std::cos(b);

    const DubinsParams params { a, b, d, sa, ca, sb, cb };

    const DubinsStruct pathLSL = _lsl(params);
    const DubinsStruct pathRSR = _rsr(params);
    const DubinsStruct pathLSR = _lsr(params);
    const DubinsStruct pathRSL = _rsl(params);

    std::vector<DubinsStruct> paths;
    if (disableCCC) {
        paths = { pathLSL, pathRSR, pathLSR, pathRSL };
    }
    else {
        DubinsStruct pathLRL;
        DubinsStruct pathRLR;
        pathRLR = _rlr(params);
        pathLRL = _lrl(params);
        paths = { pathLSL, pathRSR, pathLSR, pathRSL, pathRLR, pathLRL };
    }

    if (double rhoCompare = rhoMin() * 0.00001; std::abs(d) < rhoCompare && std::abs(a) < rhoCompare && std::abs(b) < rhoCompare) {
        if (double dist2d = std::max(std::abs(qi().x - qf().x), std::abs(qi().y - qf().y)); dist2d < rhoCompare) {
            DubinsStruct pathC;
            pathC = _c();
            paths = { pathC };
        }
    }

    ranges::sort(paths, [](const auto &x1, const auto &x2){ return x1.length < x2.length; });

    if (minLength < 0) {
        _maneuver = paths[0];
    }
    else if (const auto first = ranges::find_if(paths, [minLength](const auto &p){ return p.length >= minLength; }); first != paths.end()) {
        _maneuver = *first;
    }
}

State2d DubinsManeuver2d::_getPositionInSegment(const double offset, const State2d& qi, const char caseType) {
    State2d q = { 0, 0, 0 };
    switch (caseType) {
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
        default:
            throw std::runtime_error("Unknown segment case type");
    }
    return q;
}

DubinsStruct DubinsManeuver2d::_lsl(const DubinsParams& params) const {
    const double aux = atan2(params.cb - params.ca, params.d + params.sa - params.sb);
    const double t = Utility::mod2pi(-params.a + aux);
    const double p = sqrt(2 + params.d * params.d - 2 * cos(params.a - params.b) + 2 * params.d * (params.sa - params.sb));
    const double q = Utility::mod2pi(params.b - aux);
    const double length = (t + p + q) * rhoMin();
    const std::string caseType = "LSL";
    return DubinsStruct {t, p, q, length, caseType};
}

DubinsStruct DubinsManeuver2d::_rsr(const DubinsParams& params) const {
    const double aux = atan2(params.ca - params.cb, params.d - params.sa + params.sb);
    const double t = Utility::mod2pi(params.a - aux);
    const double p = sqrt(2 + params.d * params.d - 2 * cos(params.a - params.b) + 2 * params.d * (params.sb - params.sa));
    const double q = Utility::mod2pi(Utility::mod2pi(-params.b) + aux);
    const double length = (t + p + q) * rhoMin();
    const std::string caseType = "RSR";
    return DubinsStruct {t, p, q, length, caseType};
}

DubinsStruct DubinsManeuver2d::_lsr(const DubinsParams& params) const {
    double t, p, q;
    const double aux1 = -2 + params.d * params.d + 2 * cos(params.a - params.b) + 2 * params.d * (params.sa + params.sb);
    if (aux1 > 0) {
        p = sqrt(aux1);
        const double aux2 = atan2(-params.ca - params.cb, params.d + params.sa + params.sb) - atan(-2 / p);
        t = Utility::mod2pi(-params.a + aux2);
        q = Utility::mod2pi(-Utility::mod2pi(params.b) + aux2);
    }
    else {
        t = std::numeric_limits<double>::infinity();
        p = std::numeric_limits<double>::infinity();
        q = std::numeric_limits<double>::infinity();
    }
    const double length = (t + p + q) * rhoMin();
    const std::string caseType = "LSR";
    return DubinsStruct {t, p, q, length, caseType};
}

DubinsStruct DubinsManeuver2d::_rsl(const DubinsParams& params) const {
    double t, p, q;
    const double aux1 = params.d * params.d - 2 + 2 * cos(params.a - params.b) - 2 * params.d * (params.sa + params.sb);
    if (aux1 > 0) {
        p = sqrt(aux1);
        const double aux2 = atan2(params.ca + params.cb, params.d - params.sa - params.sb) - atan(2 / p);
        t = Utility::mod2pi(params.a - aux2);
        q = Utility::mod2pi(Utility::mod2pi(params.b) - aux2);
    }
    else {
        t = std::numeric_limits<double>::infinity();
        p = std::numeric_limits<double>::infinity();
        q = std::numeric_limits<double>::infinity();
    }
    const double length = (t + p + q) * rhoMin();
    const std::string caseType = "RSL";
    return DubinsStruct {t, p, q, length, caseType};
}

DubinsStruct DubinsManeuver2d::_rlr(const DubinsParams& params) const {
    double t, p, q;
    const double aux = (6 - params.d * params.d + 2 * cos(params.a - params.b) + 2 * params.d * (params.sa - params.sb)) / 8;
    if (std::abs(aux) <= 1) {
        p = Utility::mod2pi(-acos(aux));
        t = Utility::mod2pi(params.a - atan2(params.ca-params.cb, params.d-params.sa+params.sb) + p/2);
        q = Utility::mod2pi(params.a - params.b - t + p);
    }
    else {
        t = std::numeric_limits<double>::infinity();
        p = std::numeric_limits<double>::infinity();
        q = std::numeric_limits<double>::infinity();
    }
    const double length = (t + p + q) * rhoMin();
    const std::string caseType = "RLR";
    return DubinsStruct {t, p, q, length, caseType};
}

DubinsStruct DubinsManeuver2d::_lrl(const DubinsParams& params) const {
    double t, p, q;
    const double aux = (6 - params.d * params.d + 2 * cos(params.a - params.b) + 2 * params.d * (-params.sa + params.sb)) / 8;;
    if (std::abs(aux) <= 1) {
        p = Utility::mod2pi(-acos(aux));
        t = Utility::mod2pi(-params.a + atan2(-params.ca+params.cb, params.d+params.sa-params.sb) + p/2);
        q = Utility::mod2pi(params.b - params.a - t + p);
    }
    else {
        t = std::numeric_limits<double>::infinity();
        p = std::numeric_limits<double>::infinity();
        q = std::numeric_limits<double>::infinity();
    }
    const double length = (t + p + q) * rhoMin();
    const std::string caseType = "LRL";
    return DubinsStruct {t, p, q, length, caseType};
}

DubinsStruct DubinsManeuver2d::_c() const {
    const double t = 0;
    const double p = 2 * std::numbers::pi;
    const double q = 0;
    const double length = (t + p + q) * rhoMin();
    const std::string caseType = "RRR";
    return DubinsStruct {t, p, q, length, caseType};
}

State2d DubinsManeuver2d::getCoordinatesAt(const double offset) const {
    const double noffset = offset / rhoMin();
    const State2d qir = { 0, 0, qi().theta };
    State2d q, q1, q2;

    const double l1 = maneuver().t;
    const double l2 = maneuver().p;

    if (noffset < l1) {
        q = _getPositionInSegment(noffset, qir, maneuver().caseType.at(0));
    }
    else if(noffset < (l1 + l2)) {
        q1 = _getPositionInSegment(l1, qir, maneuver().caseType.at(0));
        q = _getPositionInSegment(noffset - l1, q1, maneuver().caseType.at(1));
    }
    else {
        q1 = _getPositionInSegment(l1, qir, maneuver().caseType.at(0));
        q2 = _getPositionInSegment(l2, q1, maneuver().caseType.at(1));
        q = _getPositionInSegment(noffset - l1 - l2, q2, maneuver().caseType.at(2));
    }

    q.x = q.x * rhoMin() + qi().x;
    q.y = q.y * rhoMin() + qi().y;
    q.theta = Utility::mod2pi(q.theta);

    return q;
}

std::vector<State2d> DubinsManeuver2d::getSamplingPoints(const double res) const {
    const int numPoints = static_cast<int>(std::floor(maneuver().length / res)) + 1;
    std::vector<State2d> points(numPoints);

    for (int i = 0; i < numPoints; i++) {
        points.push_back(getCoordinatesAt(i * res));
    }

    return points;
}