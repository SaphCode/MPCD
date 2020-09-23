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
        Wall(const Rectangle rect, const ForceType type);

        bool isInBounds(const PhysicalObject& o) const override;
        Eigen::Vector2d getCollisionPoint(const PhysicalObject& o) const override;

    private:
        const Rectangle _bounds;
    };
}
#endif // !WALL_H