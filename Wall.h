#pragma once

#ifndef WALL_H
#define WALL_H

#include "IObstacle.h"
#include "InteractingBody.h"
#include "Constants.h"

namespace MPCD {
    class Wall :
        public IObstacle,
        public InteractingBody
    {
    public:
        Wall(double yPos) :
            InteractingBody(
                std::numeric_limits<double>::infinity(), // infinite mass 
                Eigen::Vector2d((Constants::x_max - Constants::x_0) / 2, yPos), // pos is (center, y)
                Eigen::Vector2d(0, 0),
                BodyType::WALL), // vel is (0,0)
            m_yPos(yPos)
        {

        }

        ~Wall() {

        }

        bool isInBounds(const Body& o) const override;
        Eigen::Vector2d getOvershoot(const Body& o) const override;

        void interact(InteractingBody& b) override;

        bool contains(Eigen::Vector2d point) const override;
        bool occupies(std::pair<int, int> index, Eigen::Vector2d shift, double cell_dim) const override;

    private:
        double m_yPos;
    };
}
#endif // !WALL_H