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


void MPCD::timestep(std::vector<Particle>& particles, Xoshiro& rg_shift_x, Xoshiro& rg_shift_y, Xoshiro& rg_angle)
{
	double cell_dim = MPCD::Constants::Grid::cell_dim;

	Vector2d shift(rg_shift_x.next(), rg_shift_y.next());
	Xoshiro sign(-1, 1);

	std::tuple<std::map<int, double>, std::map<int, double>, std::map<int, int>> totalCellVelocitiesAndParticles = moveAndPrepare(particles, shift);
	std::map<int, double> totalCellVelocityX = std::get<0>(totalCellVelocitiesAndParticles);
	std::map<int, double> totalCellVelocityY = std::get<1>(totalCellVelocitiesAndParticles);
	std::map<int, int> particlesPerCell = std::get<2>(totalCellVelocitiesAndParticles);

	std::tuple<std::map<int, double>, std::map<int, double>, std::map<int, double>> cellMeanVelocityAndRotationAngle = calculateCellQuantities(totalCellVelocityX, totalCellVelocityY, particlesPerCell, rg_angle);
	std::map<int, double> meanCellVelocityX = std::get<0>(cellMeanVelocityAndRotationAngle);
	std::map<int, double> meanCellVelocityY = std::get<1>(cellMeanVelocityAndRotationAngle);
	std::map<int, double> cellRotationAngles = std::get<2>(cellMeanVelocityAndRotationAngle);

	updateVelocity(particles, shift, meanCellVelocityX, meanCellVelocityY, cellRotationAngles, sign);
}

/* O(N) */
std::tuple<std::map<int, double>, std::map<int, double>, std::map<int, int>> MPCD::moveAndPrepare(std::vector<Particle>& particles, Eigen::Vector2d shift) {
	std::map<int, double> totalCellVelocityX;
	std::map<int, double> totalCellVelocityY;
	std::map<int, int> numParticles;
	int cols = MPCD::Constants::Grid::cols;

	for (auto it = particles.begin(); it != particles.end(); ++it) {
		it->move(); // move const

		Vector2i cell_index = it->shift(shift); // shift // const
		int linear_index = MPCD::Grid::convertToLinearIndex(cell_index, cols);
		Vector2d vel = it->getVelocity();
		totalCellVelocityX[linear_index] += vel[0]; // add vel to total cell velocity
		totalCellVelocityY[linear_index] += vel[1]; // add vel to total cell velocity
		numParticles[linear_index] += 1; // add 1 to particle of cell
	}

	return std::make_tuple(totalCellVelocityX, totalCellVelocityY, numParticles);
}

/* O(N) (because num cells should be proportional to particles) */
std::tuple<std::map<int, double>, std::map<int, double>, std::map<int, double>> MPCD::calculateCellQuantities(std::map<int, double> totalCellVelocityX, std::map<int, double> totalCellVelocityY, std::map<int, int> numParticles, Xoshiro& rg_angle) {
	std::map<int, double> rotationAngles;
	for (auto const& [key, val] : totalCellVelocityX) {
		totalCellVelocityX[key] /= numParticles[key];
		totalCellVelocityY[key] /= numParticles[key];
		rotationAngles[key] = rg_angle.next();
	}
	return std::make_tuple(totalCellVelocityX, totalCellVelocityY, rotationAngles);
}

/* O(N) */
void MPCD::updateVelocity(std::vector<Particle>& particles, Eigen::Vector2d shift, std::map<int, double> meanCellVelocityX, std::map<int, double> meanCellVelocityY, std::map<int, double> rotationAngles, Xoshiro & sign) {
	int cols = MPCD::Constants::Grid::cols;
	for (auto it = particles.begin(); it != particles.end(); ++it) {
		Vector2i index = it->getCellIndex(it->getPosition()); // const
		int linearIndex = MPCD::Grid::convertToLinearIndex(index, cols); // const
		Vector2d meanCellVelocity(meanCellVelocityX[linearIndex], meanCellVelocityY[linearIndex]);
		double s = sign.next();
		if (s < 0) {
			s = -1;
		}
		else if (s >= 0) {
			s = 1;
		}
		it->updateVelocity(meanCellVelocity, s * rotationAngles[linearIndex]); // const
		it->shift(-shift);
	}
}
