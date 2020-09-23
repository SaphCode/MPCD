#pragma once

#ifndef RECTANGLE_H
#define RECTANGLE_H

#include "Shape.h"
#include "PhysicalObject.h"

class Rectangle : public Shape
{
public:
	Rectangle(const Eigen::Vector2d lowerLeft, const double width, const double height) 
		: x0(lowerLeft[0]), xMax(lowerLeft[0] + width), y0(lowerLeft[1]), yMax(lowerLeft[1] + height) {
	}
	const double x0;
	const double xMax;
	const double y0;
	const double yMax;
};
#endif // !RECTANGLE_H
