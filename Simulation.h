#pragma once

#ifndef MPCD_H
#define MPCD_H

#include "Particle.h"
#include "Xoshiro.h"
#include <Eigen/Dense>
#include <boost/unordered/unordered_map.hpp>

#include <map>

namespace MPCD {
	class Simulation {
	private:
		boost::unordered::unordered_map<std::pair<int, int>, Eigen::Vector2d> _totalCellVelocities;
		boost::unordered::unordered_map<std::pair<int, int>, Eigen::Vector2d> _meanCellVelocities;
		boost::unordered::unordered_map<std::pair<int, int>, int> _numCellParticles;
		boost::unordered::unordered_map<std::pair<int, int>, double> _cellRotationAngle;
		bool _draw = false;
		std::vector<Particle> _particles;
		double _cell_dim;
		double _time_lapse;
		Xoshiro _rg_shift_x;
		Xoshiro _rg_shift_y;
		Xoshiro _rg_angle;
		Xoshiro _rg_sign;
		int _w;

		void reset(int t);
	public:
		Simulation(std::vector<Particle>& particles, bool draw);
		~Simulation() {}
		/* One timestep */
		void timestep(int t);
		/* O(N)
		moves the particles, and shifts them. the total velocities of cells are figured out with the shifted particle position.*/
		void moveShiftPrepare(Eigen::Vector2d shift, int t);
		/* O(N) (because num cells should be proportional to particles)
		calulates cell means and rotation angles. */
		void calculateCellQuantities();
		/* O(N)
		updates the particles velocities according to their shifted position, and flips the sign of the rotation angle in 1/2 of cases. After
		the collision the particles are shifted back to their original positions. */
		void updateVelocity(Eigen::Vector2d shift);
	};
	
}

#endif // !MPCD_H
