#ifndef GEOMETRY_H
#define GEOMETRY_H

struct State2d
{
    double x, y, theta; 
};

struct State3d : public State2d
{
    double z, gamma;
};

#endif //