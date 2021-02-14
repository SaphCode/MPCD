#include "Monomer.h"
#include "Constants.h"

#include <iostream>

#define _USE_MATH_DEFINES
#include <math.h>

void Monomer::interact(InteractingBody& b)
{
	BodyType otherType = b.getType();

	Eigen::Vector2d pos = b.getTorusPosition();
	if (otherType == BodyType::WALL) {
		double yWall = pos[1];
		double y = getTorusPosition()[1];
		Eigen::Vector2d rel = Eigen::Vector2d(0, yWall - y);
		Eigen::Vector2d f = truncLennardJonesWall(rel, MPCD::Constants::monomerMonomer_interaction_tuning, m_diameter);
		m_effect += capForce(f);
		// no effect on wall
	}
	else if (otherType == BodyType::OBSTACLE) {
		Eigen::Vector2d rel = pos - getTorusPosition();
		Eigen::Vector2d f = truncshiftedLennardJones(rel, MPCD::Constants::monomerMonomer_interaction_tuning, (m_diameter + 2 * MPCD::Obstacles::radius) / 2);
		m_effect += capForce(f);
		// no effect on obstacle
	}
	else {
		throw std::exception();
	}
}

Eigen::Vector2d Monomer::getRelPositionTorus(Eigen::Vector2d otherPos)
{
	Eigen::Vector2d naivePos = otherPos - m_pos;
	/*double xSpan = MPCD::Constants::x_max - MPCD::Constants::x_0;
	if (m_pos[0] - otherPos[0] >= MPCD::Constants::x_max / 2) {
		naivePos[0] = MPCD::Constants::x_max - m_pos[0] + otherPos[0];
	}
	else if (otherPos[0] - m_pos[0] >= MPCD::Constants::x_max/2) {
		naivePos[0] = -1 * (MPCD::Constants::x_max - otherPos[0] + m_pos[0]);
	}*/
	return naivePos;
}


Eigen::Vector2d Monomer::capForce(Eigen::Vector2d f)
{
	double max_cell_movement = 0.1; // 0.05 cell at one md timestep is allowed
	double maxForce = 2.0 * max_cell_movement * m_mass / (MPCD::Constants::md_timestep * MPCD::Constants::md_timestep);
	double force = f.stableNorm();
	Eigen::Vector2d cappedForce = force <= maxForce ? f : maxForce * f/force;
	return cappedForce;
}

void Monomer::monomerInteraction(Eigen::Vector2d rel, double tuning, double diameter) {
	m_effect += capForce(truncshiftedLennardJones(rel, tuning, diameter));
}

void Monomer::nonlinearSpring(Eigen::Vector2d rel, double tuning, double diameter) {
	const double d = rel.stableNorm();
	const double R0 = MPCD::Constants::monomer_bond_length;
	const double k = MPCD::Constants::monomer_spring_constant;
	if (d < R0) {
		Eigen::Vector2d f = k * tuning / (diameter * diameter) * R0 * R0 * d / (1 - d * d / (R0 * R0)) * rel/d;
		m_effect += capForce(f);
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
	/*const double max_f_abs = 20 * 0.1 / (MPCD::Constants::md_timestep * MPCD::Constants::md_timestep);
	if (f_abs > max_f_abs) { // 20 * c / delta_t^2
		f_abs = max_f_abs;
	}*/
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
		/*const double max_f_abs = 20 * 0.1 / (MPCD::Constants::md_timestep * MPCD::Constants::md_timestep);
		if (f_abs > max_f_abs) { // 20 * c / delta_t^2
			f_abs = max_f_abs;
		}*/
	}
	return f;
}

void Monomer::move(const double timelapse)
{
	// position update logic
	m_oldPosition = m_pos;
	InteractingBody::move(timelapse);
}

bool Monomer::testMove(const double timelapse)
{
	
	return false;
}

void Monomer::collide(Eigen::Vector2d mean_cell_velocity, double rotationAngle, double temperatureScalingFactor) {
	Eigen::Matrix2d rotationMatrix;
	rotationMatrix(0, 0) = std::cos(rotationAngle);
	rotationMatrix(1, 0) = std::sin(rotationAngle);
	rotationMatrix(0, 1) = -std::sin(rotationAngle);
	rotationMatrix(1, 1) = std::cos(rotationAngle);

	m_vel = mean_cell_velocity + temperatureScalingFactor * rotationMatrix * (m_vel - mean_cell_velocity);

}
