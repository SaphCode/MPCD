#pragma once

#ifndef OBSTACLE_H
#define OBSTACLE_H
#include "Body.h"
#include <Eigen/Dense>
#include "Shape.h"
namespace MPCD {
	class IObstacle
	{
	public:
		/*
		Obstacle(double mass, Eigen::Vector2d pos, Eigen::Vector2d vel, ForceType type) 
			: PhysicalObject(mass, pos, vel, type) {

		}
		*/
		//virtual std::map<std::pair<int, int>, bool> occupied();
		virtual bool isInBounds(const Body& o) const = 0;
		virtual Eigen::Vector2d getOvershoot(const Body& o) const = 0;
		virtual ~IObstacle() {}
		//std::map<std::pair<int, int>, bool> _occupied;

	};
}
#endif // !OBSTACLE_H
