#pragma once

#include "InteractingBody.h"
class ConstForce :
    public InteractingBody
{
public:
    ConstForce(Eigen::Vector2d acceleration);
    void interact(InteractingBody& b) override;
private:
    Eigen::Vector2d m_force;
};

