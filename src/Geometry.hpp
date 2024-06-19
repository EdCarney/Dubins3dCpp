#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <cmath>

struct Point2d {
    double x, y;

    Point2d operator+(const Point2d& p) const { return Point2d { x + p.x, y + p.y }; }
    Point2d operator-(const Point2d& p) const { return Point2d { x - p.x, y - p.y }; }
    [[nodiscard]] double distanceTo(const Point2d& p) const { return std::hypot(x - p.x, y - p.y); }
};

struct Point3d {
    double x, y, z;
}; 

struct State2d : Point2d {
    double theta;
};

struct State3d : Point3d {
    double theta, gamma;
};

#endif // GEOMETRY_H