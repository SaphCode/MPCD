#pragma once

#ifndef CIRCLE_H
#define CIRCLE_H

#include "Shape.h"

class Circle : public Shape
{
public:
	Circle(const Eigen::Vector2d centerPosition, const double radius)
		: _centerPosition(centerPosition), _r(radius) {
	}
	double _r;
	Eigen::Vector2d _centerPosition;

};
#endif // !SHAPE_H
