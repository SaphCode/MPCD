#pragma once

#ifndef GRID_H
#define GRID_H

#include <Eigen/Dense>
#include <map>

namespace MPCD {
	namespace Grid {
		static std::map<int, Eigen::Vector2d> totalVelocityPerCell;
		static std::map<int, int> totalParticlesPerCell;
		/* Calculates the mean cell velocity, rotation angle and total number of particles (shifted) per cell.
		* @return 2 maps: mean velocity and rotation angle.

		//std::tuple<std::map<int, Eigen::Vector2d>, std::map<int, double>> calculateCellValues(std::vector<MPCD::Particle> particles);

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
			@param cols: the number of cols of your object
			@returns linear index
		*/
		int convertToLinearIndex(Eigen::Vector2i index, int cols);
	}
}
#endif
/*
namespace MPCD {
	class Grid
	{
	public:
		Grid(int minNumberPerCell);
		~Grid();

		
		
		
		
		double getCellDim();
		double getMaxShift();
		int getMinParticlesPerCell();
		int getMinNumCells();
		int getWantedNumCells();
		int getRows();
		int getCols();

	private:
		int min_particles_per_cell;
		int min_num_cells = MPCD::Constants::number / min_particles_per_cell;
		int wanted_num_cells = min_num_cells * 2;
		double cell_dim; // cutting dim by x has x^2 effect on area, and therefore on expected particle number.
		double max_shift;
		int rows;
		int cols;
	};
}
#endif
*/