#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <cmath>

struct Point2d {
    double x = 0, y = 0;

    Point2d operator+(const Point2d& p) const { return Point2d { x + p.x, y + p.y }; }
    Point2d operator-(const Point2d& p) const { return Point2d { x - p.x, y - p.y }; }
    [[nodiscard]] double distanceTo(const Point2d& p) const { return std::hypot(x - p.x, y - p.y); }
};

struct Point3d {
    double x = 0, y = 0, z = 0;
}; 

struct State2d : Point2d {
    double theta = 0;
};

struct State3d : Point3d {
    double theta = 0, gamma = 0;
};

#endif // GEOMETRY_H