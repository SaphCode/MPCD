#pragma once

#ifndef MPCD_H
#define MPCD_H

#include <cmath>

namespace MPCD {
	constexpr int number = 10000;
	constexpr int seed = 34534;
	namespace Pipe {
		constexpr double x_0 = 0;
		constexpr double y_0 = 0;
		constexpr double x_max = 1.0;
		constexpr double y_max = 1.0;
		constexpr double width = x_max - x_0;
		constexpr double height = y_max - y_0;
		constexpr double aspect_ratio = width / height;
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
#endif

