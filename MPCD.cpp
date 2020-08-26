// MPCD.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <Eigen/Dense>
#include "Particle.h"
#include "MPCD.h"
#include "Grid.h"
#include "Constants.h"
#include <filesystem>
#include <fstream>
#include <map>

#define _USE_MATH_DEFINES
#include <math.h>

using namespace Eigen;
using namespace MPCD;

int main()
{
	std::vector<Particle> particles;
	particles.reserve(MPCD::Constants::number);

	/* dont worry the numbers are just seeds */
	Xoshiro xs_xpos(0.0, 1.0);
	Xoshiro xs_ypos(0.0, 1.0);
	Xoshiro xs_xvel(-0.01, 0.01);
	Xoshiro xs_yvel(-0.01, 0.01);

	for (int i = 0; i < MPCD::Constants::number; i++) {
		double xs_x = xs_xpos.next();
		double xs_y = xs_ypos.next();
		double xs_vx = xs_xvel.next();
		double xs_vy = xs_yvel.next();

		Vector2d pos(xs_x, xs_y);
		Vector2d vel(xs_vx, xs_vy);

		Particle xs_p(pos, vel);

		particles.push_back(xs_p);
	}

	double time_lapse = MPCD::Constants::time_lapse;

	Xoshiro rg_angle(0.0, 2 * M_PI);
	double max_shift = MPCD::Constants::Grid::max_shift;
	Xoshiro rg_shift_x(-max_shift, max_shift);
	Xoshiro rg_shift_y(-max_shift, max_shift);

	int timesteps = MPCD::Constants::timesteps;

	for (int t = 0; t < timesteps; t++) {
		timestep(particles, rg_shift_x, rg_shift_y, rg_angle);
	}

}

void MPCD::timestep(std::vector<Particle>& particles, Xoshiro& rg_shift_x, Xoshiro& rg_shift_y, Xoshiro& rg_angle)
{
	std::map<int, Vector2d> meanCellVelocities;
	std::map<int, int> numParticles;
	std::map<int, double> rotationAngles;
	std::map<int, bool> calculationDone;

	double cell_dim = MPCD::Constants::Grid::cell_dim;

	moveAndShift(particles, meanCellVelocities, numParticles, rg_shift_x, rg_shift_y);
	calculateCellQuantities(meanCellVelocities, numParticles, rotationAngles, rg_angle);
	updateVelocity(particles, meanCellVelocities, rotationAngles);

}

/* O(N) */
void MPCD::moveAndShift(std::vector<Particle>& particles, std::map<int, Vector2d> &meanCellVelocities, std::map<int, int> &numParticles, Xoshiro& rg_shift_x, Xoshiro& rg_shift_y) {
	for (auto it = particles.begin(); it != particles.end(); ++it) {
		it->move(); // move const

		Vector2d shift(rg_shift_x.next(), rg_shift_y.next());
		Vector2i cell_index = it->shift(shift); // shift // const
		int linear_index = MPCD::Grid::convertToLinearIndex(cell_index);
		Vector2d vel = it->getVelocity();
		meanCellVelocities[linear_index] += vel; // add vel to total cell velocity
		numParticles[linear_index] += 1; // add 1 to particle of cell
	}
}

/* O(N) (because num cells should be proportional to particles) */
void MPCD::calculateCellQuantities(std::map<int, Vector2d>& meanCellVelocities, std::map<int, int>& numParticles, std::map<int, double>& rotationAngles, Xoshiro& rg_angle) {
	for (auto const& [key, val] : meanCellVelocities) {
		meanCellVelocities[key] = val / numParticles[key];
		rotationAngles[key] = rg_angle.next();
	}
}

/* O(N) */
void MPCD::updateVelocity(std::vector<Particle>& particles, std::map<int, Vector2d> meanCellVelocities, std::map<int, double> rotationAngles) {
	for (auto it = particles.begin(); it != particles.end(); ++it) {
		Vector2i index = it->getCellIndex(it->getPosition()); // const
		int linearIndex = MPCD::Grid::convertToLinearIndex(index); // const
		it->updateVelocity(meanCellVelocities[linearIndex], rotationAngles[linearIndex]); // const
	}
}
