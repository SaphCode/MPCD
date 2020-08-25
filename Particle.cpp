#include "Particle.h"
#include <iostream>
#include "MPCD.h"
#include <cmath>
//#include "Grid.h"

using namespace MPCD;

Particle::Particle(Eigen::Vector2d position, Eigen::Vector2d velocity) {
	_position = position;
	_velocity = velocity;
}

Particle::Particle() {}

Particle::~Particle() {}

/* If you make the variables PRIVATE */
Eigen::Vector2d Particle::getPosition() {
	return _position;
}

Eigen::Vector2d Particle::getVelocity() {
	return _velocity;
}

/*
Eigen::Vector2i Particle::getCellIndex() {
	return _cell_index;
}
*/
/*
void Particle::setPosition(Eigen::Vector2d newPosition) {
	_position = newPosition;
}

void Particle::setVelocity(Eigen::Vector2d newVelocity) {
	_velocity = newVelocity;
}
*/

/* If you make the variables PRIVATE */

void Particle::move(double time_step) {
	_position += time_step * _velocity;
}

void Particle::updateVelocity(Eigen::Vector2d mean_cell_velocity, double rotationAngle) {
	Eigen::Matrix2d rotationMatrix;
	rotationMatrix << cos(rotationAngle), -sin(rotationAngle),
					sin(rotationAngle), cos(rotationAngle);
	_velocity = mean_cell_velocity + rotationMatrix * (_velocity - mean_cell_velocity);
}

Eigen::Vector2i Particle::getCellIndex(Eigen::Vector2d shiftedPosition, double cell_dim) {
	int i = std::round(std::floor(shiftedPosition(0) / cell_dim));
	int j = std::round(std::floor(shiftedPosition(1) / cell_dim));
	Eigen::Vector2i cell_index(i, j);
	return cell_index;
	
	/*
	std::cout << "cell_index: (" << _cell_index(0) << ", " << _cell_index(1) << ")" << std::endl;
	std::cout << "cell_dim: " << cell_dim << std::endl;
	std::cout << "(update) shiftedPos: (" << shiftedPosition(0) << "," << shiftedPosition(1) << ")" << std::endl;
	*/
}

Eigen::Vector2i Particle::shift(Eigen::Vector2d amount, double cell_dim) {
	Eigen::Vector2d shiftedPosition = _position + amount;
	//std::cout << "shiftedPos: (" << shiftedPosition(0) << "," << shiftedPosition(1) << ")" << std::endl;
	return getCellIndex(shiftedPosition, cell_dim);
}


std::ostream& MPCD::operator<<(std::ostream& output, const Particle& p) {
	output << "Particle:\nLocation(\n" << p._position << "\n)\nVelocity(\n" << p._velocity << "\n)\n" << std::endl;
	return output;
}
