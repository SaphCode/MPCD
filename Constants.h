#pragma once

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>

namespace MPCD {
	namespace Constants {
		//constexpr int number = Grid::average_particles_per_cell * (Pipe::x_max / Grid::cell_dim) * (Pipe::y_max / Grid::cell_dim); // will be a func of Grid
		constexpr int seed = 5646542;
		
		constexpr int timesteps = 200;
		constexpr int average_particles_per_cell = 10;
		constexpr double cell_dim = 1;
		constexpr double k_boltzmann = 1;
		constexpr double temperature = 2; // = 36Celsius
		constexpr double viscosity = 1.6735;
		constexpr double particle_mass = 1;
		const double unit_of_time = std::sqrt((particle_mass * std::pow(cell_dim, 2) / (k_boltzmann * temperature)));
		const double force_const = cell_dim / (unit_of_time * unit_of_time); // kg (particle mass) * m (1) / s^2 (unit of time) about 10^-20, so a = 10^3 seems reasonable, but its really small force
		const double time_lapse = 0.1;// * unit_of_time; // will be a func of diff parameters
		constexpr double x_0 = 0;
		constexpr double y_0 = 0;
		constexpr double x_max = 40;
		constexpr double y_max = 5;
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