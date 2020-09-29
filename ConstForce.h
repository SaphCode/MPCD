#pragma once

#include "InteractingBody.h"
class ConstForce :
    public InteractingBody
{
public:
    ConstForce(Eigen::Vector2d acceleration);
    Eigen::Vector2d interact(InteractingBody& b);
private:
    Eigen::Vector2d m_acceleration;
};

