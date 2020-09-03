#pragma once

#ifndef MPCD_H
#define MPCD_H

#include "Particle.h"
#include "Grid.h"
#include "Xoshiro.h"
#include <tuple>
#include <Eigen/Dense>
#include "VectorCompare.h"

#include <map>

namespace MPCD {
	vectorMap initMap();

	/* One timestep */
	vectorMap timestep(std::vector<Particle>& particles, Xoshiro& rg_shift_x, Xoshiro& rg_shift_y, Xoshiro& rg_angle);
	/* O(N)
	moves the particles, and shifts them. the total velocities of cells are figured out with the shifted particle position.*/
	std::tuple<vectorMap, std::map<Eigen::Vector2i, int>> moveAndPrepare(std::vector<Particle>& particles, Eigen::Vector2d shift);
	/* O(N) (because num cells should be proportional to particles)
	calulates cell means and rotation angles. */
	std::tuple<vectorMap, std::map<Eigen::Vector2i, double>> calculateCellQuantities(vectorMap totalCellVelocities, std::map<Eigen::Vector2i, int> numParticles, Xoshiro& rg_angle);
	/* O(N)
	updates the particles velocities according to their shifted position, and flips the sign of the rotation angle in 1/2 of cases. After
	the collision the particles are shifted back to their original positions. */
	void updateVelocity(std::vector<Particle>& particles, Eigen::Vector2d shift, vectorMap meanCellVelocities, std::map<Eigen::Vector2i, double> rotationAngles, Xoshiro& sign);
}

#endif // !MPCD_H
