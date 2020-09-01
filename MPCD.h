#pragma once

#ifndef MPCD_H
#define MPCD_H

#include "Particle.h"
#include "Grid.h"
#include "Xoshiro.h"
#include <tuple>

#include <map>

namespace MPCD {
	/* One timestep */
	void timestep(std::vector<Particle>& particles, Xoshiro& rg_shift_x, Xoshiro& rg_shift_y, Xoshiro& rg_angle);
	/* O(N)
	moves the particles, and shifts them. the total velocities of cells are figured out with the shifted particle position.*/
	std::tuple<std::map<int, double>,std::map<int, double>, std::map<int, int>> moveAndPrepare(std::vector<Particle>& particles, Eigen::Vector2d shift);
	/* O(N) (because num cells should be proportional to particles)
	calulates cell means and rotation angles. */
	std::tuple<std::map<int, double>, std::map<int, double>, std::map<int, double>> calculateCellQuantities(std::map<int, double> totalCellVelocityX, std::map<int, double> totalCellVelocityY, std::map<int, int> numParticles, Xoshiro& rg_angle);
	/* O(N)
	updates the particles velocities according to their shifted position, and flips the sign of the rotation angle in 1/2 of cases. After
	the collision the particles are shifted back to their original positions. */
	void updateVelocity(std::vector<Particle>& particles, Eigen::Vector2d shift, std::map<int, double> meanCellVelocityX, std::map<int, double> meanCellVelocityY, std::map<int, double> rotationAngles, Xoshiro & sign);
}

#endif // !MPCD_H
