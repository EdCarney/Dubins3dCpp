#include <math.h>

#ifndef GEOMETRY_H
#define GEOMETRY_H

struct Point2d
{
    double x, y;

    Point2d operator+(const Point2d& p) const { return Point2d { x + p.x, y + p.y }; }
    Point2d operator-(const Point2d& p) const { return Point2d { x - p.x, y - p.y }; }
    double distanceTo(const Point2d& p) const { return hypot(x - p.x, y - p.y); }
};

struct Point3d
{
    double x, y, z;

    Point3d operator+(const Point3d& p) const;
    Point3d operator-(const Point3d& p) const;
    double distanceTo(const Point3d& p) const;
}; 

struct State2d : public Point2d
{
    double theta;
};

struct State3d : public Point3d
{
    double theta, gamma;
};

#endif //