#ifndef DUBINS_MANEUVER_2D
#define DUBINS_MANEUVER_2D

#include <string>
#include <vector>
#include "Geometry.hpp"

struct DubinsStruct {
    double t, p, q, length;
    std::string caseType;
};

struct DubinsParams {
    double a, b, d, sa, ca, sb, cb;
};

class DubinsManeuver2d {
    State2d _qi, _qf;
    double _rhoMin;
    DubinsStruct _maneuver;

    void _generateManeuver(double minLength, bool disableCCC);
    [[nodiscard]] static State2d _getPositionInSegment(double offset, const State2d& qi, char caseType);
    [[nodiscard]] DubinsStruct _lsl(const DubinsParams& params) const;
    [[nodiscard]] DubinsStruct _rsr(const DubinsParams& params) const;
    [[nodiscard]] DubinsStruct _lsr(const DubinsParams& params) const;
    [[nodiscard]] DubinsStruct _rsl(const DubinsParams& params) const;
    [[nodiscard]] DubinsStruct _rlr(const DubinsParams& params) const;
    [[nodiscard]] DubinsStruct _lrl(const DubinsParams& params) const;
    [[nodiscard]] DubinsStruct _c() const;

    public:
        DubinsManeuver2d() = default;
        DubinsManeuver2d(const State2d& qi, const State2d& qf, const double rhoMin = 1, double minLength = -1, bool disableCCC = false)
            : _qi{qi}, _qf{qf}, _rhoMin{rhoMin}, _maneuver{0,0,0,std::numeric_limits<double>::infinity(), ""}
            { _generateManeuver(minLength, disableCCC); }
        [[nodiscard]] double rhoMin() const { return _rhoMin; }
        [[nodiscard]] const State2d& qi() const { return _qi; }
        [[nodiscard]] const State2d& qf() const { return _qf; };
        [[nodiscard]] const DubinsStruct& maneuver() const { return _maneuver; }
        void setManeuver(DubinsStruct maneuver) { _maneuver = _maneuver; }

        [[nodiscard]] State2d getCoordinatesAt(double offset) const;
        [[nodiscard]] std::vector<State2d> getSamplingPoints(double res = 0.1) const;
};

#endif //DUBINS_MANEUVER_2D