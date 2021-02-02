#pragma once

#include <cmath>

namespace MPCD {
	namespace Constants {
		//constexpr int number = Grid::average_particles_per_cell * (Pipe::x_max / Grid::cell_dim) * (Pipe::y_max / Grid::cell_dim); // will be a func of Grid
		constexpr int seed = 234234;
		
		constexpr int timesteps = 3000;
		constexpr int average_particles_per_cell = 10;
		constexpr double cell_dim = 1;
		constexpr double k_boltzmann = 1;
		constexpr double temperature = 1; // = 36Celsius
		constexpr double monomerMonomer_interaction_tuning = 0.5 * k_boltzmann * temperature;
		constexpr double monomer_diameter = cell_dim;
		constexpr double monomer_bond_length = cell_dim;
		constexpr double monomer_spring_constant = 5.0 * 10.0 * 10.0 * 10.0 * k_boltzmann * temperature / (cell_dim * cell_dim);
		constexpr int num_monomers = 50;
		constexpr int num_md_timesteps = 50;
		constexpr double particle_mass = 1;
		constexpr double monomer_mass = average_particles_per_cell * particle_mass; // avg particle mass per cell
		const double unit_of_time = std::sqrt((particle_mass * std::pow(cell_dim, 2) / (k_boltzmann * temperature)));
		const double const_force = 0.01; // 1 at the moment, kg (particle mass) * m (1) / s^2 (unit of time) about 10^-20, so a = 10^3 seems reasonable, but its really small force
		const double time_lapse = 0.1 * unit_of_time;// * unit_of_time; // will be a func of diff parameters
		const double md_timestep = 1.0 / (double)num_md_timesteps * time_lapse; // TODO: is this ok?
		constexpr double x_0 = 0;
		constexpr double y_0 = 0;
		constexpr double x_max = 400;
		constexpr double y_max = 20;
	}
	namespace Obstacles {
		constexpr double radius = 2.5;
		constexpr double x_start = 100.0;
		constexpr double x_end = 200.0;
		constexpr double x_space = x_end - x_start;
		constexpr int num_per_row = (int)(x_space / (2.0 * 2.0 * radius));
		constexpr double center_to_center_spacing = 2 * radius + 2 * radius; // (->) + -- + (<-)
		constexpr double y_center_lower = Constants::y_0 + radius;
		constexpr double y_center_upper = Constants::y_max - radius;
	}
}