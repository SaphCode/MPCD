#pragma once

#ifndef SHAPE_H
#define SHAPE_H

#include <Eigen/Dense>

class Shape
{
public:
	bool isInBounds();
	Eigen::Vector2d getCollisionPoint();
};
#endif // !SHAPE_H
