#pragma once

#ifndef OBSTACLE_H
#define OBSTACLE_H
#include "PhysicalObject.h"
#include <Eigen/Dense>
#include "ForceType.h"
#include "Shape.h"
namespace MPCD {
	class Obstacle : public PhysicalObject
	{
	public:
		Obstacle(double mass, Eigen::Vector2d pos, Eigen::Vector2d vel, ForceType type) 
			: PhysicalObject(mass, pos, vel, type) {

		}
		//virtual std::map<std::pair<int, int>, bool> occupied();
		virtual bool isInBounds(const PhysicalObject& o) const = 0;
		virtual Eigen::Vector2d getOvershoot(const PhysicalObject& o) const = 0;
	private:
		//std::map<std::pair<int, int>, bool> _occupied;

	};
}
#endif // !OBSTACLE_H
