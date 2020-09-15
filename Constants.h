#pragma once

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>

namespace MPCD {
	namespace Constants {
		//constexpr int number = Grid::average_particles_per_cell * (Pipe::x_max / Grid::cell_dim) * (Pipe::y_max / Grid::cell_dim); // will be a func of Grid
		constexpr int seed = 34534211;
		
		constexpr int timesteps = 100;
		constexpr int average_particles_per_cell = 10;
		constexpr double cell_dim = 1;
		constexpr double k_boltzmann = 1.38064852e-23;
		constexpr double temperature = 309.15; // = 36Celsius
		constexpr double particle_mass = 2.988e-26;
		const double unit_of_time = std::sqrt((particle_mass * std::pow(cell_dim, 2) / (k_boltzmann * temperature)));
		const double time_lapse = 0.1 * unit_of_time; // will be a func of diff parameters
		constexpr double x_0 = 0;
		constexpr double y_0 = 0;
		constexpr double x_max = 400;
		constexpr double y_max = 20;
		/*namespace Pipe {
			
			constexpr double width = x_max - x_0;
			constexpr double height = y_max - y_0;
			constexpr double aspect_ratio = width / height;
		}*/
		/*namespace Grid {
			
			constexpr int cell_dim = 1;
			constexpr int num_rows = (Pipe::y_max - Pipe::y_0) / cell_dim;
			constexpr int num_cols = (Pipe::x_max - Pipe::x_0) / cell_dim;
			constexpr int num_cells = MPCD::Constants::number / average_particles_per_cell;
			//constexpr int wanted_num_cells = min_num_cells * 2;
			constexpr double max_shift = cell_dim / 2;
		}*/
	}
}
#endif