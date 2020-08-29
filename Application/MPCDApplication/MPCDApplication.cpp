// MPCDApplication.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <Eigen/Dense>
#include "Particle.h"
#include "Constants.h"
#include "Xoshiro.h"
#include "MPCD.h"

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

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
