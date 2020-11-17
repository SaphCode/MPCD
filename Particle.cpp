#include "Particle.h"
#include <iostream>
#include <cmath>
#include "Constants.h"

using namespace MPCD;

void Particle::collide(Eigen::Vector2d mean_cell_velocity, double rotationAngle, double temperatureScalingFactor) {
	Eigen::Matrix2d rotationMatrix;
	rotationMatrix(0, 0) = std::cos(rotationAngle);
	rotationMatrix(1, 0) = std::sin(rotationAngle);
	rotationMatrix(0, 1) = -std::sin(rotationAngle);
	rotationMatrix(1, 1) = std::cos(rotationAngle);

	m_vel = m_vel + temperatureScalingFactor * rotationMatrix * (m_vel - mean_cell_velocity); // TODO: woops
}

void Particle::move(double timelapse)
{
	m_oldPosition = m_pos;
	InteractingBody::move(timelapse);
}