#pragma once

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>

namespace MPCD {
	namespace Constants {
		constexpr int number = 500000;
		constexpr int seed = 34534211;
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
}
#endif