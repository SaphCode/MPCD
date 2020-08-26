#pragma once

#ifndef MPCD_H
#define MPCD_H

#include "Particle.h"
#include "Grid.h"
#include "Xoshiro.h"

#include <map>

namespace MPCD {
	/* One timestep */
	void timestep(std::vector<Particle>& particles, Xoshiro& rg_shift_x, Xoshiro& rg_shift_y, Xoshiro& rg_angle);
	/* O(N) */
	void moveAndShift(std::vector<Particle>& particles, std::map<int, Eigen::Vector2d>& meanCellVelocities, std::map<int, int>& numParticles, Xoshiro& rg_shift_x, Xoshiro& rg_shift_y);
	/* O(N) (because num cells should be proportional to particles) */
	void calculateCellQuantities(std::map<int, Eigen::Vector2d>& meanCellVelocities, std::map<int, int>& numParticles, std::map<int, double>& rotationAngles, Xoshiro& rg_angle);
	/* O(N) */
	void updateVelocity(std::vector<Particle>& particles, std::map<int, Eigen::Vector2d> meanCellVelocities, std::map<int, double> rotationAngles);
}

#endif // !MPCD_H
