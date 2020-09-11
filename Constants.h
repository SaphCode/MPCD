#pragma once

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>

namespace MPCD {
	namespace Constants {
		constexpr int number = Grid::average_particles_per_cell * (Pipe::x_max / Grid::cell_dim) * (Pipe::y_max / Grid::cell_dim); // will be a func of Grid
		constexpr int seed = 34534211;
		constexpr double time_lapse = 0.1; // will be a func of diff parameters
		constexpr int timesteps = 1000;
		namespace Pipe {
			constexpr double x_0 = 0;
			constexpr double y_0 = 0;
			constexpr double x_max = 400;
			constexpr double y_max = 20;
			constexpr double width = x_max - x_0;
			constexpr double height = y_max - y_0;
			constexpr double aspect_ratio = width / height;
		}
		namespace Grid {
			constexpr int average_particles_per_cell = 10;
			constexpr int num_cells = MPCD::Constants::number / average_particles_per_cell;
			//constexpr int wanted_num_cells = min_num_cells * 2;
			constexpr double cell_dim = 1;
			constexpr double max_shift = cell_dim / 2;
		}
	}
}
#endif