#pragma once
#include <Eigen/Dense>
#include <map>
#include "MPCD.h"
#include "Xoshiro.h"
#include "Particle.h"

#ifndef GRID_H
#define GRID_H

namespace MPCD {
	class Grid
	{
	public:
		Grid(int minNumberPerCell);
		~Grid();

		/* Calculates the mean cell velocity, rotation angle and total number of particles (shifted) per cell.
		* @return 2 maps: mean velocity and rotation angle.
		*/
		std::tuple<std::map<Eigen::Vector2i, Eigen::Vector2d>, std::map<Eigen::Vector2i, double>> calculateCellValues(std::vector<MPCD::Particle> particles);

		/* Converts 2d indexes into linear indexes.
			Example:
			A = [[1, 2, 3]
				 [4, 5, 6]
				 [7, 8, 9]]
			A[0,0], index = (0,0)
			returns 0
			A[1,1], index = (1,1)
			returns 4
			@param index needs to be 2d!
			@returns linear index
		*/
		//static int convertToLinearIndex(Eigen::Vector2i index);

	private:
		int min_particles_per_cell;
		int min_num_cells = number / min_particles_per_cell;
		int wanted_num_cells = min_num_cells * 2;
		double cell_dim; // cutting dim by x has x^2 effect on area, and therefore on expected particle number.
		double max_shift;
		int rows;
		int cols;
		Xoshiro rg_angle;
		Xoshiro rg_shift_x;
		Xoshiro rg_shift_y;
	};
}
#endif