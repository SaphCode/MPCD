#pragma once

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>

namespace MPCD {
	namespace Constants {
		constexpr int number = 1000000;
		constexpr int seed = 3453453453;
		constexpr double time_lapse = 0.1;
		constexpr int timesteps = 1000;
		namespace Pipe {
			constexpr double x_0 = -1.0;
			constexpr double y_0 = -1.0;
			constexpr double x_max = 1.0;
			constexpr double y_max = 1.0;
			constexpr double width = x_max - x_0;
			constexpr double height = y_max - y_0;
			constexpr double aspect_ratio = width / height;
		}
		namespace Grid {
			constexpr int average_particles_per_cell = 20;
			constexpr int num_cells = MPCD::Constants::number / average_particles_per_cell;
			//constexpr int wanted_num_cells = min_num_cells * 2;
			const double cell_dim = std::sqrt(MPCD::Constants::Pipe::width * MPCD::Constants::Pipe::height / num_cells );
			const double max_shift = cell_dim / 2;
			constexpr double grid_x_shift = MPCD::Constants::Pipe::x_0;
			constexpr double grid_y_shift = MPCD::Constants::Pipe::y_0;
		}
	}
		/*
		const int x_cells = std::round(sqrt(number * 0.9 / min_particles_per_cell)); // at least 5 particles per cell. 90% to make sure. sqrt gives 1 dimension, round gives int
		const int y_cells = std::round(x_cells / Pipe::aspect_ratio);
		const int num_cells = x_cells * y_cells; // how many cells does grid have
		const double grid_width = Pipe::width * (1 + 1 / x_cells) * 1.05; // at least * (1 + 1/x_cells), should be more
		const double grid_height = Pipe::height * (1 + 1 / y_cells) * 1.05; // at least * (1 + 1/y_cells), should be more
		const double cell_dim_x = grid_width / x_cells; // this is prob not rightst
		const double cell_dim = grid_width / x_cells; // this is prob not rightst
		const double cell_size_y = grid_height / y_cells;
		const double x_0_offset = Pipe::x_0 - cell_dim / 2;
		const double y_0_offset = Pipe::y_0 - cell_dim / 2;
		*/
}

/*
namespace Detail
{
	double constexpr sqrtNewtonRaphson(double x, double curr, double prev)
	{
		return curr == prev
			? curr
			: sqrtNewtonRaphson(x, 0.5 * (curr + x / curr), curr);
	}
}

/*
* Constexpr version of the square root
* Return value:
*   - For a finite and non-negative value of "x", returns an approximation for the square root of "x"
*   - Otherwise, returns NaN

double constexpr sqrt(double x)
{
	return x >= 0 && x < std::numeric_limits<double>::infinity()
		? Detail::sqrtNewtonRaphson(x, x, 0)
		: std::numeric_limits<double>::quiet_NaN();
}
*/
#endif