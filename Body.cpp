#include "Body.h"
#include "Constants.h"

#include <iostream>


Body::Body(double mass, Eigen::Vector2d pos, Eigen::Vector2d vel)
{
	m_mass = mass;
	m_pos = pos;
	m_vel = vel;
}

Body::~Body()
{
}

void Body::move(const double timelapse)
{
	m_pos += m_vel * timelapse;
}

void Body::correctPosition(const Eigen::Vector2d newPos)
{
	m_pos = newPos;
}

void Body::collided(const Eigen::Vector2d overshoot)
{
	m_pos -= 2 * overshoot;
	m_vel = -m_vel;
}

Eigen::Vector2d Body::getPosition() const
{
	return m_pos;
}

Eigen::Vector2d Body::getTorusPosition() const {
	double xMax = MPCD::Constants::x_max;
	return Eigen::Vector2d(std::fmod(m_pos[0], xMax), m_pos[1]);
}

Eigen::Vector2d Body::getOldTorusPosition() const {
	double xMax = MPCD::Constants::x_max;
	return Eigen::Vector2d(std::fmod(m_oldPosition[0], xMax), m_oldPosition[1]);
}

Eigen::Vector2d Body::getOldPosition() const
{
	return m_oldPosition;
}

Eigen::Vector2d Body::getVelocity() const
{
	return m_vel;
}

double Body::getMass() const {
	return m_mass;
}