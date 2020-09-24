#pragma once

#include "InteractingBody.h"
class ConstForce :
    public InteractingBody
{
public:
    ConstForce(double accelerationConstant);
    Eigen::Vector2d interact(InteractingBody& b);
private:
    double m_accelerationConstant;
};

