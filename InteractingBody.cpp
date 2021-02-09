#include "InteractingBody.h"
#include <iostream>

using namespace Eigen;

void InteractingBody::move(const double timelapse)
{
	std::cout << "Force: " << m_effect << "\n";
	std::cout << "Delta r: " << m_effect / m_mass * timelapse * timelapse / 2 << "\n";
	m_pos += m_effect / m_mass * timelapse * timelapse / 2;
	/*
	std::cout << "InteractingBody\n";
	std::cout << "Effect:\n" << m_effect;
	std::cout << "\nTimelapse:\n" << timelapse;
	std::cout << "\nMass:\n" << m_mass;
	std::cout << "\nPosition:\n" << m_pos;
	*/
	Body::move(timelapse);
}

void InteractingBody::updateVelocity(const double timelapse) {
	m_vel += timelapse * m_effect / m_mass;
}

void InteractingBody::addEffect(Eigen::Vector2d effect)
{
	m_effect += effect;
}
