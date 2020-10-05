#include "Particle.h"
#include <iostream>
#include <cmath>
#include "Constants.h"
//#include "Grid.h"

using namespace MPCD;

/*
void Particle::shift(Eigen::Vector2d amount) {
	_position += amount;
}
*/

void Particle::collide(Eigen::Vector2d mean_cell_velocity, double rotationAngle, double temperatureScalingFactor) {
	Eigen::Matrix2d rotationMatrix;
	rotationMatrix(0, 0) = std::cos(rotationAngle);
	rotationMatrix(1, 0) = std::sin(rotationAngle);
	rotationMatrix(0, 1) = -std::sin(rotationAngle);
	rotationMatrix(1, 1) = std::cos(rotationAngle);
	
	double abs_m_vel = m_vel.stableNorm();
	double m_vel_x = m_vel[0];
	double m_vel_y = m_vel[1];

	m_vel = mean_cell_velocity + temperatureScalingFactor * rotationMatrix * (m_vel - mean_cell_velocity);
	
	double abs_m_vel_after = m_vel.stableNorm();
	double m_vel_x_after = m_vel[0];
	double m_vel_y_after = m_vel[1];
	//if (abs_m_vel >= 10) {

	//}
}

void Particle::move(double timelapse)
{
	InteractingBody::move(timelapse);
}

/*
bool MPCD::operator==(const Particle& lhs, const Particle& rhs) {
	return (lhs._pos == rhs._pos) && (lhs._vel == rhs._vel);
}
*/
