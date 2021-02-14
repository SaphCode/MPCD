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
	if (_isnan(timelapse) || _isnan(m_pos[0]) || _isnan(m_pos[1]) || _isnan(m_vel[0]) || _isnan(m_vel[1])) {
		std::cout << "in correctPosition:\n";
		std::cout << "timelapse: " << timelapse << "\n";
		std::cout << "Vel: " << m_vel << "\n";
		std::cout << "Mass: " << getMass() << "\nPos: " << getPosition() << "\n";
		std::cout << "Old Pos: " << getOldPosition() << "\n";
		exit(-1);
	}
	m_pos += m_vel * timelapse;
}

void Body::correctPosition(const Eigen::Vector2d newPos)
{
	if (_isnan(newPos[0]) || _isnan(newPos[1])) {
		std::cout << "in correctPosition:\n";
		std::cout << "newpos: " << newPos << "\n";
		std::cout << "Vel: " << m_vel << "\n";
		std::cout << "Mass: " << getMass() << "\nPos: " << getPosition() << "\n";
		std::cout << "Old Pos: " << getOldPosition() << "\n";
		exit(-1);
	}
	m_pos = newPos;
}

void Body::collided(const Eigen::Vector2d overshoot)
{
	if (_isnan(overshoot[0]) || _isnan(overshoot[1])) {
		std::cout << "in collided:\n";
		std::cout << "overshoot: " << overshoot << "\n";
		std::cout << "Vel: " << m_vel << "\n";
		std::cout << "Mass: " << getMass() << "\nPos: " << getPosition() << "\n";
		std::cout << "Old Pos: " << getOldPosition() << "\n";
		exit(-1);
	}
	m_pos -= 2.0 * overshoot;
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