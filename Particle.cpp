#include "Particle.h"
#include <iostream>
#include <cmath>
#include "Constants.h"
//#include "Grid.h"

using namespace MPCD;

Particle::Particle(Eigen::Vector2d position, Eigen::Vector2d velocity, double mass) : PhysicalObject(mass, position, velocity) {

}

Particle::Particle() : PhysicalObject() {}


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

void Particle::stream(double timeLapse) {
	PhysicalObject::updatePosition(timeLapse);
}

void MPCD::Particle::correctPosition(Eigen::Vector2d newPos)
{
	_pos = newPos;
}

/*
void Particle::shift(Eigen::Vector2d amount) {
	_position += amount;
}
*/

void Particle::updateVelocity(double timelapse, Eigen::Vector2d mean_cell_velocity, double rotationAngle) {
	PhysicalObject::updateVelocity(timelapse);
	Eigen::Matrix2d rotationMatrix;
	rotationMatrix << cos(rotationAngle), -sin(rotationAngle),
					sin(rotationAngle), cos(rotationAngle);
	_vel = mean_cell_velocity + rotationMatrix * (_vel - mean_cell_velocity);
}

/*
void Particle::shift(Eigen::Vector2d amount) {
	_position += amount;
}
*/
bool MPCD::operator==(const Particle& lhs, const Particle& rhs) {
	return (lhs._pos == rhs._pos) && (lhs._vel == rhs._vel);
}

/*
void Particle::_updateCellIndex() {
	double cell_dim = MPCD::Constants::Grid::cell_dim;
	double grid_shifted_x = _position(0) - MPCD::Constants::Grid::grid_x_shift;
	double grid_shifted_y = _position(1) - MPCD::Constants::Grid::grid_y_shift;
	int i = std::round(std::floor(grid_shifted_y / cell_dim));
	int j = std::round(std::floor(grid_shifted_x / cell_dim));
	Eigen::Vector2i cell_index(i, j);
	_cell_index = cell_index;
}

/*
bool Particle::operator<(const Particle& p) {
	Eigen::Vector2i other_index = p._cell_index;
	int size = _cell_index.size();
	assert(size == other_index.size());
	for (int i = 0; i < size; i++) {
		if (_cell_index[i] < other_index[i]) {
			return true;
		}
		else if (_cell_index[i] > other_index[i]) {
			return false;
		}
	}
	return false;
}
*/
