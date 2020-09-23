#pragma once

#ifndef WALL_H
#define WALL_H

#include "ImmovableObstacle.h"
#include "Rectangle.h"
#include "PhysicalObject.h"


namespace MPCD {
    class Wall : public ImmovableObstacle
    {
    public:
        Wall(const double yPos, const ForceType type);

        bool isInBounds(const PhysicalObject& o) const override;
        Eigen::Vector2d getOvershoot(const PhysicalObject& o) const override;

    private:
        const double _yPos;
    };
}
#endif // !WALL_H