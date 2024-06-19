#ifndef DUBINS_MANEUVER_3D
#define DUBINS_MANEUVER_3D

#include <vector>
#include "Geometry.hpp"
#include "Utility.hpp"
#include "DubinsManeuver2d.hpp"

struct DubinsPath {
    DubinsManeuver2d lat, lon;
    bool isEmpty() const { return &lat != NULL && &lon != NULL; }

    DubinsPath() {
        lat = *((DubinsManeuver2d*)NULL);
        lon = *((DubinsManeuver2d*)NULL);
    }
};

class DubinsManeuver3d {
    State3d _qi, _qf;
    double _rhoMin, _length;
    std::tuple<double, double> _pitchLims;
    std::vector<DubinsManeuver2d> _path;

    void _generateManeuver();
    [[nodiscard]] std::vector<DubinsManeuver2d> _tryToConstruct(double horizontalRadius) const;

    public:
        DubinsManeuver3d() : _rhoMin{0}, _length{-1} { }
        DubinsManeuver3d(const State3d& qi, const State3d& qf, const double rhoMin, const std::tuple<double, double>& pitchLims)
            : _qi{qi}, _qf{qf}, _rhoMin{rhoMin}, _length{-1}, _pitchLims{pitchLims} { _generateManeuver(); }
        [[nodiscard]] double rhoMin() const { return _rhoMin; }
        [[nodiscard]] double length() const { return _length; }
        [[nodiscard]] const State3d& qi() const { return _qi; }
        [[nodiscard]] const State3d& qf() const { return _qf; }
        [[nodiscard]] double minPitch() const { return std::get<0>(_pitchLims); }
        [[nodiscard]] double maxPitch() const { return std::get<1>(_pitchLims); }
        [[nodiscard]] const std::vector<DubinsManeuver2d>& path() const { return _path; }
        [[nodiscard]] std::vector<State3d> computeSampling(int numSamples = 1000) const;
        void setPath(const std::vector<DubinsManeuver2d>& path) { _path = path; }
};

#endif //DUBINS_MANEUVER_3D