#pragma once

#include <Eigen/Dense>
#include "Force.h"
class PhysicalObject
{
public:
	double getEffect(PhysicalObject o, double timelapse);
private:
	double _mass;
	Eigen::Vector2d _pos;
	Eigen::Vector2d _vel;
	Force _field;
};

