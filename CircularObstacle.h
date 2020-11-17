#pragma once

#ifndef CIRCULAR_OBSTCLE_H
#define CIRCULAR_OBSTCLE_H
#include "IObstacle.h"
#include "InteractingBody.h"

namespace MPCD {
    class CircularObstacle :
        public IObstacle,
        public InteractingBody
    {
    public:
        CircularObstacle(Eigen::Vector2d center, double radius);

        bool isInBounds(const Body& o) const override;
        Eigen::Vector2d getOvershoot(const Body& o) const;

        Eigen::Vector2d interact(InteractingBody& b) override;

        bool contains(Eigen::Vector2d point) const override;

        bool occupies(std::pair<int, int> index, double cell_dim) const override;

    private:
        Eigen::Vector2d m_center;
        const double m_radius;
    };
}
#endif // !CIRCULAR_OBSTCLE_H