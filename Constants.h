#pragma once

#include <cmath>

namespace MPCD {
	namespace Constants {
		constexpr int average_particles_per_cell = 10;
		constexpr double cell_dim = 1;
		constexpr double k_boltzmann = 1;
		constexpr double temperature = 1;
		constexpr double monomerMonomer_interaction_tuning = 0.5 * k_boltzmann * temperature;
		constexpr double monomer_diameter = cell_dim;
		constexpr double monomer_bond_length = 1.5 * monomer_diameter;
		constexpr double monomer_spring_constant = 30.0 * monomerMonomer_interaction_tuning / (monomer_diameter * monomer_diameter);//30.0
		constexpr int num_monomers = 50;
		constexpr int num_md_timesteps = 400;
		constexpr double particle_mass = 1.0;
		constexpr double monomer_mass = average_particles_per_cell * particle_mass; // avg particle mass per cell
		constexpr double unit_of_time = 1.0; // a * mass_solvent/sqrt(kB T) // std::sqrt((particle_mass* std::pow(cell_dim, 2) / (k_boltzmann * temperature)))
		constexpr double const_force = 0.03 * particle_mass * cell_dim * unit_of_time * unit_of_time;
		constexpr double time_lapse = 0.1 * unit_of_time;
		constexpr double md_timestep = 1.0 / (double)num_md_timesteps * time_lapse;
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
		constexpr double y_center_lower = Constants::y_0 + 2.5 + radius;
		constexpr double y_center_upper = Constants::y_max - 2.5 - radius;
	}
}