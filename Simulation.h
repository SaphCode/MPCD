#pragma once

#ifndef MPCD_H
#define MPCD_H

#include "Particle.h"
#include "Xoshiro.h"
#include <Eigen/Dense>

#include <map>

namespace MPCD {
	class Simulation {
	private:
		std::map<int, Eigen::Vector2d> _totalCellVelocities;
		std::map<int, Eigen::Vector2d> _meanCellVelocities;
		std::map<int, int> _numCellParticles;
		std::map<int, double> _cellRotationAngle;
		std::vector<Particle> _particles;
		double _cell_dim;
		double _time_lapse;
		Xoshiro _rg_shift_x;
		Xoshiro _rg_shift_y;
		Xoshiro _rg_angle;
		Xoshiro _rg_sign;

		void init();
		void reset();
	public:
		Simulation(std::vector<Particle>& particles);
		~Simulation() {}
		/* One timestep */
		std::map<int, Eigen::Vector2d> timestep();
		/* O(N)
		moves the particles, and shifts them. the total velocities of cells are figured out with the shifted particle position.*/
		void moveShiftPrepare(Eigen::Vector2d shift);
		/* O(N) (because num cells should be proportional to particles)
		calulates cell means and rotation angles. */
		std::tuple<std::map<int, Eigen::Vector2d>, std::map<int, double>> calculateCellQuantities(std::map<int, Eigen::Vector2d> totalCellVelocities, std::map<int, int> numCellParticles);
		/* O(N)
		updates the particles velocities according to their shifted position, and flips the sign of the rotation angle in 1/2 of cases. After
		the collision the particles are shifted back to their original positions. */
		void updateVelocity(Eigen::Vector2d shift);
	};
	
}

#endif // !MPCD_H
