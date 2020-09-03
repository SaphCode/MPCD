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
#include <chrono>

#define _USE_MATH_DEFINES
#include <math.h>

using namespace Eigen;
using namespace MPCD;
using namespace std::chrono;

vectorMap MPCD::timestep(std::vector<Particle>& particles, Xoshiro& rg_shift_x, Xoshiro& rg_shift_y, Xoshiro& rg_angle)
{
	double cell_dim = MPCD::Constants::Grid::cell_dim;

	Vector2d shift(rg_shift_x.next(), rg_shift_y.next());
	Xoshiro sign(-1, 1);

	std::tuple<vectorMap, std::map<Eigen::Vector2i, int>> totalCellVelocitiesAndParticles = moveAndPrepare(particles, shift);
	vectorMap totalCellVelocities = std::get<0>(totalCellVelocitiesAndParticles);
	std::map<Eigen::Vector2i, int> particlesPerCell = std::get<1>(totalCellVelocitiesAndParticles);

	std::tuple<vectorMap, std::map<Eigen::Vector2i, double>> cellMeanVelocityAndRotationAngle = calculateCellQuantities(totalCellVelocities, particlesPerCell, rg_angle);
	vectorMap meanCellVelocities = std::get<0>(cellMeanVelocityAndRotationAngle);
	std::map<Eigen::Vector2i, double> cellRotationAngles = std::get<1>(cellMeanVelocityAndRotationAngle);

	updateVelocity(particles, shift, meanCellVelocities, cellRotationAngles, sign);
	return meanCellVelocities;
}

/* O(N) */
std::tuple<vectorMap, std::map<Eigen::Vector2i, int>> MPCD::moveAndPrepare(std::vector<Particle>& particles, Eigen::Vector2d shift) {
	vectorMap totalCellVelocities;
	std::map<int, double> totalCellVelocityY;
	std::map<Eigen::Vector2i, int> numParticles;

	std::cout << "--Moving Particles .." << "\n";

	auto start_loop = high_resolution_clock::now();
	for (auto it = particles.begin(); it != particles.end(); ++it) {
		it->move(); // move const

		Vector2i cell_index = it->shift(shift); // shift // const
		Vector2d vel = it->getVelocity();
		auto start_timing_map = high_resolution_clock::now();
		bool newlyInserted = totalCellVelocities.insert(std::make_pair(cell_index, vel)).second;
		if (!newlyInserted) {
			totalCellVelocities[cell_index] = totalCellVelocities[cell_index] + vel;
		}
		numParticles[cell_index] += 1;
		auto finish_timing_map = high_resolution_clock::now();

		auto ms_map = duration_cast<microseconds>(finish_timing_map - start_timing_map);
		//std::cout << "Took " << ms_map.count() << " mus for map insert/calculate. Inserted: " << newlyInserted << "\n";
	}
	auto end_loop = high_resolution_clock::now();
	auto ms_loop = duration_cast<microseconds>(end_loop - start_loop);
	int particles_size = particles.size();
	std::cout << "Took " << ms_loop.count() << " mus for the whole loop. This means about " << ms_loop.count() / particles_size << " mus per particle." << "\n";

	return std::make_tuple(totalCellVelocities, numParticles);
}

/* O(N) (because num cells should be proportional to particles) */
std::tuple<vectorMap, std::map<Eigen::Vector2i, double>> MPCD::calculateCellQuantities(vectorMap totalCellVelocities, std::map<Eigen::Vector2i, int> numParticles, Xoshiro& rg_angle) {

	std::cout << "--Calculating Means .." << "\n";
	std::map<Eigen::Vector2i, double> rotationAngles;
	for (auto const& [key, val] : totalCellVelocities) {
		totalCellVelocities[key] /= numParticles[key];
		rotationAngles[key] = rg_angle.next();
	}
	return std::make_tuple(totalCellVelocities, rotationAngles);
}

/* O(N) */
void MPCD::updateVelocity(std::vector<Particle>& particles, Eigen::Vector2d shift, vectorMap meanCellVelocities, std::map<Eigen::Vector2i, double> rotationAngles, Xoshiro & sign) {
	std::cout << "--Updating velocities .." << "\n";
	for (auto it = particles.begin(); it != particles.end(); ++it) {
		Vector2i index = it->getCellIndex(it->getPosition()); // const
		double s = sign.next();
		if (s < 0) {
			s = -1;
		}
		else if (s >= 0) {
			s = 1;
		}
		it->updateVelocity(meanCellVelocities[index], s * rotationAngles[index]); // const
		it->shift(-shift);
	}
}