#pragma once

#ifndef PHYSICALOBJECT_H
#define PHYSICALOBJECT_H


#include <Eigen/Dense>
#include "ForceType.h"
#include <algorithm>
class PhysicalObject
{
public:
	PhysicalObject(double mass, Eigen::Vector2d pos, Eigen::Vector2d vel, ForceType type) : _mass(mass), _pos(pos), _vel(vel), _f(type) {}
	PhysicalObject(double mass, Eigen::Vector2d pos, Eigen::Vector2d vel) : _mass(mass), _pos(pos), _vel(vel), _f(ForceType::NO_FORCE) {}
	PhysicalObject() : _mass(0), _pos(0,0), _vel(0,0), _f(ForceType::NO_FORCE) {}
	PhysicalObject(ForceType type) : _mass(0), _pos(0, 0), _vel(0, 0), _f(type) {}
	Eigen::Vector2d getPosition() const;
	Eigen::Vector2d getVelocity() const;
	void updateVelocity(double timelapse);
	void updatePosition(double timelapse);
	void registerObject(PhysicalObject& o);
	void unregisterObject(PhysicalObject& o);
protected:
	double _mass;
	Eigen::Vector2d _pos;
	Eigen::Vector2d _vel;
private:
	Eigen::Vector2d calculateForceOnThis() const;
	Eigen::Vector2d getEffect(PhysicalObject o) const;
	ForceType _f;
	std::vector<PhysicalObject*> _registeredObjects;
};

#endif // !PHYSICALOBJECT_H