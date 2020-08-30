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
	/* O(N) */
	std::tuple<std::map<int, double>,std::map<int, double>, std::map<int, int>> moveAndPrepare(std::vector<Particle>& particles, Xoshiro& rg_shift_x, Xoshiro& rg_shift_y);
	/* O(N) (because num cells should be proportional to particles) */
	std::tuple<std::map<int, double>, std::map<int, double>, std::map<int, double>> calculateCellQuantities(std::map<int, double> totalCellVelocityX, std::map<int, double> totalCellVelocityY, std::map<int, int> numParticles, Xoshiro& rg_angle);
	/* O(N) */
	void updateVelocity(std::vector<Particle>& particles, std::map<int, double> meanCellVelocityX, std::map<int, double> meanCellVelocityY, std::map<int, double> rotationAngles);
}

#endif // !MPCD_H
