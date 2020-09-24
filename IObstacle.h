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
		virtual bool isInBounds(const Body& o) const {
			return true;
		}
		virtual Eigen::Vector2d getOvershoot(const Body& o) const {
			return Eigen::Vector2d(0, 0);
		}
		//virtual ~IObstacle();
	private:
		//std::map<std::pair<int, int>, bool> _occupied;

	};
}
#endif // !OBSTACLE_H
