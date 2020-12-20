#include "Monomer.h"
#include "Constants.h"

#include <iostream>

#define _USE_MATH_DEFINES
#include <math.h>

void Monomer::interact(InteractingBody& b)
{
	BodyType otherType = b.getType();

	Eigen::Vector2d pos = b.getPosition();
	Eigen::Vector2d rel = pos - m_pos;
	if (otherType == BodyType::MONOMER) {
		// monomer monomer interaction
		Eigen::Vector2d f = 1/2 * truncLennardJones(rel, MPCD::Constants::monomerMonomer_interaction_tuning, 2 * m_diameter); // 1/2 b/c of double counting
		m_effect += f;
		b.addEffect(-f);
	}
	else if (otherType == BodyType::PARTICLE) {
		// particle - monomer interaction
		Eigen::Vector2d f = truncLennardJones(rel, MPCD::Constants::solventMonomer_interaction_tuning, m_diameter);
		m_effect += f;
		b.addEffect(-f);
	}
	else if (otherType == BodyType::WALL) {
		pos = b.getPosition();
		double yWall = pos[1];
		double y = m_pos[1];
		rel = Eigen::Vector2d(0, yWall - y);
		Eigen::Vector2d f = truncLennardJones(rel, MPCD::Constants::obstacleMonomer_interaction_tuning, m_diameter);
		m_effect += f;
		// no effect on wall
	}
	else if (otherType == BodyType::OBSTACLE) {
		Eigen::Vector2d f = truncLennardJones(rel, MPCD::Constants::obstacleMonomer_interaction_tuning, (m_diameter + 2 * MPCD::Obstacles::radius));
		m_effect += f;
		// no effect on obstacle
	}
	else {
		throw std::exception();
	}
}

Eigen::Vector2d Monomer::truncLennardJones(Eigen::Vector2d rel, double tuning, double diameter) {
	double d = rel.stableNorm();
	Eigen::Vector2d f(0, 0);
	if (d < std::pow(2, 1 / 6) * diameter) {
		double f_abs = 4 * tuning * (-12 * std::pow(diameter, 12) / std::pow(d, 13) + 6 * std::pow(diameter, 6) / std::pow(d, 7));
		f[0] = rel[0] / d * f_abs;
		f[1] = rel[1] / d * f_abs;
	}
	if (f.stableNorm() > 1000.0) {
		std::cout << f << std::endl;
		f[0] = f[0] > 0 ? std::min(f[0], 1000.0) : std::max(f[0], -1000.0);
		f[1] = f[1] > 0 ? std::min(f[1], 1000.0) : std::max(f[1], -1000.0);
		std::cout << f << std::endl;
	}
	return f;
}

void Monomer::move(const double timelapse)
{
	// position update logic
	m_oldPosition = m_pos;
	m_pos += timelapse * m_vel + 1 / 2 * timelapse * timelapse * m_effect / m_mass;
}

void Monomer::collide(Eigen::Vector2d mean_cell_velocity, double rotationAngle, double temperatureScalingFactor) {
	Eigen::Matrix2d rotationMatrix;
	rotationMatrix(0, 0) = std::cos(rotationAngle);
	rotationMatrix(1, 0) = std::sin(rotationAngle);
	rotationMatrix(0, 1) = -std::sin(rotationAngle);
	rotationMatrix(1, 1) = std::cos(rotationAngle);

	m_vel = mean_cell_velocity + temperatureScalingFactor * rotationMatrix * (m_vel - mean_cell_velocity);

}
