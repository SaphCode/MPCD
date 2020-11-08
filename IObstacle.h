#pragma once

#include "Body.h"
#include <Eigen/Dense>
namespace MPCD {
	class IObstacle
	{
	public:
		virtual bool isInBounds(const Body& o) const = 0;
		virtual Eigen::Vector2d getOvershoot(const Body& o) const = 0;
		virtual ~IObstacle() {}
	};
}
