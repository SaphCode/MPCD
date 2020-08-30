// MPCDApplication.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <Eigen/Dense>
#include "Particle.h"
#include "Constants.h"
#include "Xoshiro.h"
#include "MPCD.h"
#include "Grid.h"

#define _USE_MATH_DEFINES
#include <math.h>

using namespace Eigen;
using namespace MPCD;

int main()
{
	double xvel = 1;
	double yvel = 1;
	Vector2d vel(xvel, yvel);

	std::vector<Particle> particles;
	int num = MPCD::Constants::number;
	particles.reserve(MPCD::Constants::Grid::rows * MPCD::Constants::Grid::cols);

	double cell_dim = MPCD::Constants::Grid::cell_dim;
	Vector2d startPos(cell_dim / 10, cell_dim / 10);
	std::map<int, double> totalCellVelocitiesX;
	std::map<int, double> totalCellVelocitiesY;
	std::map<int, int> cellParticles;

	Xoshiro rg_angle(0.0, 2 * 3.141); // not important in this test

	/*
	* This block should place 1 particle in each cell. The mean of every cell should then be the velocity of the particle itself.
	*/
	for (int i = 0; i < MPCD::Constants::Grid::rows; i++) {
		for (int j = 0; j < MPCD::Constants::Grid::cols; j++) {
			Vector2d offset(j * cell_dim, i * cell_dim);
			Vector2d pos = startPos + offset;
			Particle p(pos, vel);
			particles.push_back(p);

			int linearIndex = MPCD::Grid::convertToLinearIndex(p.getCellIndex(p.getPosition()), MPCD::Constants::Grid::cols);

			Vector2d vel = p.getVelocity();
			std::cout << "Adding " << vel[0] << " to " << totalCellVelocitiesX[linearIndex] << std::endl;
			totalCellVelocitiesX[linearIndex] += vel[0];
			std::cout << "This gives: " << totalCellVelocitiesX[linearIndex] << std::endl;
			std::cout << "Adding " << vel[1] << " to " << totalCellVelocitiesY[linearIndex] << std::endl;
			totalCellVelocitiesY[linearIndex] += vel[1];
			std::cout << "This gives: " << totalCellVelocitiesY[linearIndex] << std::endl;
			cellParticles[linearIndex] += 1;
		}
	}

	//ASSERT_EQ(totalCellVelocitiesX.size(), totalCellVelocitiesY.size());
	for (auto const& [key, val] : totalCellVelocitiesX) {
		std::cout << "(" << key << ", " << val << ")" << std::endl;
		std::cout << "Particles in Cell: " << cellParticles[key] << std::endl;
		totalCellVelocitiesX[key] /= cellParticles[key];
		totalCellVelocitiesY[key] /= cellParticles[key];
	}

	//updateVelocity(particles, meanCellVelocities, cellRotationAngles);
	std::tuple <std::map<int, double>, std::map<int, double>, std::map<int, double >> cellMeanVelocityAndRotationAngle = MPCD::calculateCellQuantities(totalCellVelocitiesX, totalCellVelocitiesY, cellParticles, rg_angle);
	std::map<int, double> meanCellVelocityX = std::get<0>(cellMeanVelocityAndRotationAngle);
	std::map<int, double> meanCellVelocityY = std::get<1>(cellMeanVelocityAndRotationAngle);
	std::map<int, double> cellRotationAngles = std::get<2>(cellMeanVelocityAndRotationAngle);

	for (auto const& [key, val] : meanCellVelocityX) {
		double epsilon = 0.0001;
		//std::cout << "(" << key << ", " << val << ")" << std::endl;
		//ASSERT_LT(std::abs(val[0] - xvel) / xvel, epsilon);
		//ASSERT_LT(std::abs(val[1] - yvel) / yvel, epsilon);
	}


	
	/*


	std::vector<Particle> particles;
	particles.reserve(MPCD::Constants::number);

	/* dont worry the numbers are just seeds 
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
	*/
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
