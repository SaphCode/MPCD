#pragma once
#include "Obstacle.h"
#include <limits>
namespace MPCD {
    class ImmovableObstacle :
        public Obstacle
    {
    public:
        ImmovableObstacle(const Eigen::Vector2d pos, const ForceType type) : Obstacle(std::numeric_limits<double>::infinity(), pos, Eigen::Vector2d(0, 0), type) {

        }
        virtual bool Obstacle::isInBounds(const PhysicalObject& o) const = 0;
        virtual Eigen::Vector2d Obstacle::getCollisionPoint(const PhysicalObject& o) const = 0;
    };
}


