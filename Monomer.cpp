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
		// monomer - monomer interaction
		Eigen::Vector2d f = 1/2 * truncshiftedLennardJones(rel, MPCD::Constants::monomerMonomer_interaction_tuning, m_diameter); // 1/2 b/c of double counting
		m_effect += f;
		b.addEffect(-f);
	}
	else if (otherType == BodyType::PARTICLE) {
		// no interaction with particles
		/*Eigen::Vector2d f = truncLennardJones(rel, MPCD::Constants::solventMonomer_interaction_tuning, m_diameter);
		m_effect += f;
		b.addEffect(-f);*/
	}
	else if (otherType == BodyType::WALL) {
		pos = b.getPosition();
		double yWall = pos[1];
		double y = m_pos[1];
		rel = Eigen::Vector2d(0, yWall - y);
		Eigen::Vector2d f = truncLennardJonesWall(rel, MPCD::Constants::monomerMonomer_interaction_tuning, m_diameter);
		m_effect += f;
		// no effect on wall
	}
	else if (otherType == BodyType::OBSTACLE) {
		Eigen::Vector2d f = truncshiftedLennardJones(rel, MPCD::Constants::monomerMonomer_interaction_tuning, (m_diameter + 2 * MPCD::Obstacles::radius) / 2);
		m_effect += f;

		// no effect on obstacle
	}
	else {
		throw std::exception();
	}
}

Eigen::Vector2d Monomer::getRelPositionTorus(Eigen::Vector2d otherPos)
{
	Eigen::Vector2d naivePos = otherPos - m_pos;
	double xSpan = MPCD::Constants::x_max - MPCD::Constants::x_0;
	if (m_pos[0] - otherPos[0] >= MPCD::Constants::x_max / 2) {
		naivePos[0] = MPCD::Constants::x_max - m_pos[0] + otherPos[0];
	}
	else if (otherPos[0] - m_pos[0] >= MPCD::Constants::x_max) {
		naivePos[0] = -1 * (MPCD::Constants::x_max - otherPos[0] + m_pos[0]);
	}
	
	return naivePos;
}


void Monomer::monomerInteraction(Eigen::Vector2d rel, double tuning, double diameter) {
	m_effect += truncshiftedLennardJones(rel, tuning, diameter);
}

void Monomer::nonlinearSpring(Eigen::Vector2d rel, double tuning, double diameter) {
	const double d = rel.stableNorm();
	const double R0 = MPCD::Constants::monomer_bond_length;
	const double k = MPCD::Constants::monomer_spring_constant;
	if (d < R0) {
		m_effect += k * tuning/(diameter*diameter) * R0*R0 * d / (1 - d * d / (R0 * R0)) * rel.normalized();
	}
	else {
		m_effect += 1000 * k * d * rel.normalized();
	}
}


Eigen::Vector2d Monomer::truncshiftedLennardJones(Eigen::Vector2d rel, double tuning, double diameter) {
	const double r = rel.stableNorm();
	const double r_end = std::pow(2, 1 / 6) * diameter;
	Eigen::Vector2d f(0, 0);
	if (r <= r_end) {
		Eigen::Vector2d lj = lennardJones(rel, tuning, diameter);
		Eigen::Vector2d r_hat = rel/r;
		//Eigen::Vector2d shift = lennardJones(r_hat * r_end, tuning, diameter);
		f = lj; //- shift
	}
	return f;
}

Eigen::Vector2d Monomer::lennardJones(Eigen::Vector2d rel, double tuning, double diameter) {
	double r = rel.stableNorm();
	Eigen::Vector2d r_hat = rel / r;

	double f_abs = 4 * tuning * (-12 * std::pow(diameter, 12) / std::pow(r, 13) - -6 * std::pow(diameter, 6) / std::pow(r, 7));
	Eigen::Vector2d f = f_abs * r_hat;
	return f;
}

Eigen::Vector2d Monomer::truncLennardJonesWall(Eigen::Vector2d rel, double tuning, double diameter) {
	const double r = rel.stableNorm();
	const double r_end = std::pow(2 / 5, 1 / 6) * diameter;
	Eigen::Vector2d r_hat = rel / r;
	Eigen::Vector2d f(0, 0);
	if (r <= r_end) {
		double f_abs = tuning * (2 / 15 * -9 * std::pow(diameter, 9) / std::pow(r, 10) - -3 * std::pow(diameter, 3) / std::pow(r, 4));
		f = f_abs * r_hat;
	}
	return f;
}

void Monomer::move(const double timelapse)
{
	// position update logic
	m_oldPosition = m_pos;
	InteractingBody::move(timelapse);
}

void Monomer::collide(Eigen::Vector2d mean_cell_velocity, double rotationAngle, double temperatureScalingFactor) {
	Eigen::Matrix2d rotationMatrix;
	rotationMatrix(0, 0) = std::cos(rotationAngle);
	rotationMatrix(1, 0) = std::sin(rotationAngle);
	rotationMatrix(0, 1) = -std::sin(rotationAngle);
	rotationMatrix(1, 1) = std::cos(rotationAngle);

	m_vel = mean_cell_velocity + temperatureScalingFactor * rotationMatrix * (m_vel - mean_cell_velocity);

}
