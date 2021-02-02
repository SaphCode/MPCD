#include "Body.h"
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
	/*
	std::cout << "Body\n";
	std::cout << "Vel:\n" << m_vel;
	std::cout << "\nTimelapse:\n" << timelapse;
	*/
}

void Body::correctPosition(const Eigen::Vector2d newPos)
{
	m_pos = newPos;
}

void Body::collided(const Eigen::Vector2d overshoot)
{
	m_pos = m_pos - 2 * overshoot;
	m_vel = -m_vel;
}

Eigen::Vector2d Body::getPosition() const
{
	return m_pos;
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