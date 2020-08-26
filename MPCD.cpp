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

#define _USE_MATH_DEFINES
#include <math.h>

using namespace Eigen;
using namespace MPCD;

int main()
{
	Grid g(MPCD::Constants::min_particles_per_cell);
	std::vector<Particle> particles;
	particles.reserve(MPCD::Constants::number);

	/* dont worry the numbers are just seeds */
	Xoshiro xs_xpos(0.0, 1.0);
	Xoshiro xs_ypos(0.0, 1.0);
	Xoshiro xs_xvel(-0.01, 0.01);
	Xoshiro xs_yvel(-0.01, 0.01);
	Xoshiro angles(0.0, 2 * M_PI);
	Xoshiro rg_angle = angles;

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
	double max_shift = g.getMaxShift();
	Xoshiro rg_shift_x(-max_shift, max_shift);
	Xoshiro rg_shift_y(-max_shift, max_shift);

	int timesteps = MPCD::Constants::timesteps;

	for (int t = 0; t < timesteps; t++) {
		timestep(particles, g, time_lapse, rg_shift_x, rg_shift_y, rg_angle);
	}

	
}

void MPCD::timestep(std::vector<Particle>& particles, Grid g, double time_lapse, Xoshiro& rg_shift_x, Xoshiro& rg_shift_y, Xoshiro& rg_angle)
{
	std::map<int, Vector2d> meanCellVelocities;
	std::map<int, int> numParticles;
	std::map<int, double> rotationAngles;
	std::map<int, bool> calculationDone;

	double cell_dim = g.getCellDim();

	for (auto it = particles.begin(); it != particles.end(); ++it) {
		it->move(time_lapse);

		Vector2d shift(rg_shift_x.next(), rg_shift_y.next());
		Vector2i cell_index = it->shift(shift, cell_dim);
		int linear_index = g.convertToLinearIndex(cell_index);
		Vector2d vel = it->getVelocity();
		meanCellVelocities[linear_index] += vel;
		numParticles[linear_index] += 1;
	}

	for (auto const& [key, val] : meanCellVelocities) {
		meanCellVelocities[key] =  val / numParticles[key];
		rotationAngles[key] = rg_angle.next();
	}

	for (auto it = particles.begin(); it != particles.end(); ++it) {
		Vector2i index = it->getCellIndex(it->getPosition(), cell_dim);
		int linearIndex = g.convertToLinearIndex(index);
		it->updateVelocity(meanCellVelocities[linearIndex], rotationAngles[linearIndex]);
	}

}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file


