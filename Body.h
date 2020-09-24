#pragma once

#ifndef IBODY_H
#define IBODY_H

#include <Eigen/Dense>

class Body {
public:
	Body(double mass, Eigen::Vector2d pos, Eigen::Vector2d vel);
	~Body();
	virtual void move(const double timelapse);
	virtual void correct(const Eigen::Vector2d newPos);
	Eigen::Vector2d getPosition() const;
	virtual Eigen::Vector2d getOldPosition(const double timelapse) const;
	Eigen::Vector2d getVelocity() const;
protected:
	double m_mass;
	Eigen::Vector2d m_pos;
	Eigen::Vector2d m_vel;
};

#endif