// MPCD.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <Eigen/Dense>
#include "Particle.h"
#include "Simulation.h"
#include "Grid.h"
#include "Constants.h"
#include <filesystem>
#include <fstream>
#include <map>
#include <cmath>

#define _USE_MATH_DEFINES
#include <math.h>

using namespace Eigen;
using namespace MPCD;
using namespace std::chrono;

MPCD::Simulation::Simulation(std::vector<Particle>& particles) {
	_time_lapse = MPCD::Constants::time_lapse;
	_cell_dim = MPCD::Constants::Grid::cell_dim;
	_particles = particles;
	Xoshiro rg_shift_x(-MPCD::Constants::Grid::max_shift, MPCD::Constants::Grid::max_shift);
	Xoshiro rg_shift_y(-MPCD::Constants::Grid::max_shift, MPCD::Constants::Grid::max_shift);
	Xoshiro rg_angle(0.0, 2 * M_PI);
	Xoshiro rg_sign(-1, 1);
	_rg_shift_x = rg_shift_x;
	_rg_shift_y = rg_shift_y;
	_rg_angle = rg_angle;
	_rg_sign = rg_sign;
	init();
}

std::map<int, Vector2d> MPCD::Simulation::timestep()
{
	reset();
	Vector2d shift(_rg_shift_x.next(), _rg_shift_y.next());
	moveShiftPrepare(shift);
	calculateCellQuantities(_totalCellVelocities, _numCellParticles);
	updateVelocity(shift);
	
	return _meanCellVelocities;
}

/* O(N) */
void MPCD::Simulation::moveShiftPrepare(Eigen::Vector2d shift) {
	for (auto it = _particles.begin(); it != _particles.end(); ++it) {
		it->move(); // move const
		it->shift(shift); // shift // const
		Vector2d velocity = it->getVelocity();
		int linear_index = Grid::convertToLinearIndex(it->getCellIndex());
		_totalCellVelocities[linear_index] += velocity;
		_numCellParticles[linear_index] += 1;
		/*//it->getCellIndex();
		Eigen::Vector2d pos = it->getPosition();
		if (pos[0] > max[0]) {
			max[0] = pos[0];
		}
		else if (pos[0] < min[0]) {
			min[0] = pos[0];
		}
		if (pos[1] > max[1]) {
			max[1] = pos[1];
		}
		else if (pos[1] < min[1]) {
			min[1] = pos[1];
		}
	}
	double cell_dim = MPCD::Constants::Grid::cell_dim;

	int min_col = std::floor(min[0] / cell_dim);
	int max_col = std::floor(max[0] / cell_dim);
	int cols = max_col - min_col;
	assert(cols > 0);
	int min_row = std::floor(min[1] / cell_dim);
	int max_row = std::floor(max[1] / cell_dim);
	int rows = max_row - min_row;
	assert(rows > 0);*/
	}
}



/* O(N) (because num cells should be proportional to particles) */
std::tuple<std::map<int, Eigen::Vector2d>, std::map<int, double>> MPCD::Simulation::calculateCellQuantities(std::map<int, Eigen::Vector2d> totalCellVelocities, std::map<int, int> numCellParticles) {
	Vector2d zero(0, 0);
	for (auto const& [key, val] : totalCellVelocities) {
		Vector2d totalCellVelocity = totalCellVelocities[key];
		if (totalCellVelocity == zero) {
			_meanCellVelocities[key] = zero;
		}
		else {
			int numParticles = numCellParticles[key];
			assert(numParticles != 0);
			_meanCellVelocities[key] = totalCellVelocity / numParticles;
		}
		_cellRotationAngle[key] = _rg_angle.next();
	}
	return std::make_pair(_meanCellVelocities, _cellRotationAngle);
}

/* O(N) */
void MPCD::Simulation::updateVelocity(Eigen::Vector2d shift) {
	std::cout << "--Updating velocities .." << "\n";
	for (auto it = _particles.begin(); it != _particles.end(); ++it) {
		int linear_index = Grid::convertToLinearIndex(it->getCellIndex()); // const
		double s = _rg_sign.next();
		if (s < 0) {
			s = -1;
		}
		else if (s >= 0) {
			s = 1;
		}
		it->updateVelocity(_meanCellVelocities[linear_index], s * _cellRotationAngle[linear_index]); // const
		it->shift(-shift);
	}
}

void MPCD::Simulation::init() {
	Vector2d zero(0, 0);
	for (int i = 0; i < MPCD::Constants::Grid::max_rows; i++) {
		for (int j = 0; j < MPCD::Constants::Grid::max_cols; j++) {
			Vector2i index(i, j);
			int linearIndex = Grid::convertToLinearIndex(index);
			_totalCellVelocities.insert(std::make_pair(linearIndex, zero));
			_meanCellVelocities.insert(std::make_pair(linearIndex, zero));
		}
	}
}

void MPCD::Simulation::reset() {
	Vector2d reset_vel(0, 0);
	for (int i = 0; i < MPCD::Constants::Grid::max_rows; i++) {
		for (int j = 0; j < MPCD::Constants::Grid::max_cols; j++) {
			Vector2i index(i, j);
			int linearIndex = Grid::convertToLinearIndex(index);
			_totalCellVelocities[linearIndex] = reset_vel;
			_meanCellVelocities[linearIndex] = reset_vel;
			_numCellParticles[linearIndex] = 0;
		}
	}
}
