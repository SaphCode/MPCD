#include "pch.h"
#include "Grid.h"
#include "Constants.h"
#include "Out.h"
#include <filesystem>
#include <fstream>
#include "Locations.h"
#include <iomanip>
#include "Xoshiro.h"
#include <iostream>
#include <stdexcept>
#include "Simulation.h"

#define _USE_MATH_DEFINES
#include <math.h>

using namespace MPCD;
using namespace Eigen;

TEST(MPCD, Timestep) {
	std::vector<Particle> particles;
	particles.reserve(MPCD::Constants::number);

	/* dont worry the numbers are just seeds */
	Xoshiro xs_xpos(MPCD::Constants::Pipe::x_0, MPCD::Constants::Pipe::x_max);
	Xoshiro xs_ypos(MPCD::Constants::Pipe::y_0, MPCD::Constants::Pipe::y_max);
	Xoshiro xs_xvel(-MPCD::Constants::Pipe::width / 100, MPCD::Constants::Pipe::width / 100);
	Xoshiro xs_yvel(-MPCD::Constants::Pipe::height / 100, MPCD::Constants::Pipe::height / 100);

	for (int i = 0; i < MPCD::Constants::number; i++) {
		double xs_x = xs_xpos.next();
		double xs_y = xs_ypos.next();
		double xs_vx = xs_xvel.next();
		double xs_vy = xs_yvel.next();

		Vector2d pos(xs_x, xs_y);
		Vector2d vel(xs_vx, xs_vy);

		Particle xs_p(pos, vel);

		particles.push_back(xs_p);
	}

	double time_lapse = MPCD::Constants::time_lapse;

	std::filesystem::path cwd = std::filesystem::current_path();
	
	Out out(cwd.string() + l_data);

	int timesteps = MPCD::Constants::timesteps;
	ASSERT_GT(timesteps, 0);

	int w = -1;
	if (timesteps < 99) {
		w = 2;
	}
	else if (timesteps < 999) {
		w = 3;
	}
	else if (timesteps < 9999) {
		w = 4;
	}
	else if (timesteps < 99999) {
		w = 5;
	}
	else {
		throw std::exception("timesteps too large or too short: " + timesteps);
	}

	bool drawParticles = true;
	bool drawCells = true;
	
	Simulation sim(particles);

	for (int t = 0; t < timesteps; t++) {
		std::map<int, Eigen::Vector2d> meanCellVelocities = sim.timestep();
		if (drawParticles) {
			std::stringstream s;
			std::stringstream av;
			s << std::setfill('0') << std::setw(w) << t;
			av << "_av" << MPCD::Constants::Grid::average_particles_per_cell;
			std::string filename = "particles_timestep_" + s.str() + av.str() + ".csv";
			out.writeToOut(particles, filename);
		}
		if (drawCells) {
			std::stringstream s;
			std::stringstream av;
			s << std::setfill('0') << std::setw(w) << t;
			av << "_av" << MPCD::Constants::Grid::average_particles_per_cell;
			std::string filename = "cellvalues_timestep_" + s.str() + av.str() + ".csv";
			out.writeToOut(meanCellVelocities, filename, "i,j,vx,vy");
		}
	}
	
	

	/*

	Grid g(min_particles_per_cell);

	std::tuple<std::map<int, Vector2d>, std::map<int, double>> origMaps = g.calculateCellValues(particles);
	std::map<int, Vector2d> origMeanVelocities = std::get<0>(origMaps);
	std::map<int, double> origRotationAngles = std::get<1>(origMaps);

	ASSERT_EQ(origMeanVelocities.size(), origRotationAngles.size());

	int steps = 1;
	//for (int step = 0; step < steps; step++) {
	for (auto it = particles.begin(); it != particles.end(); ++it) {
		it->move(time_step);
	}

	std::tuple<std::map<int, Vector2d>, std::map<int, double>> maps = g.calculateCellValues(particles);
	std::map<int, Vector2d> meanCellVelocities = std::get<0>(maps);
	std::map<int, double> rotationAngles = std::get<1>(maps);

	double cell_dim = g.getCellDim();
	for (auto it = particles.begin(); it != particles.end(); ++it) {
		Vector2i cellIndex = it->getCellIndex(it->getPosition(), cell_dim);
		int linearIndex = g.convertToLinearIndex(cellIndex);
		it->updateVelocity(meanCellVelocities[linearIndex], rotationAngles[linearIndex]);
	}

	ASSERT_EQ(true, false); // this is here to remind yourself to fix this test. in mpcd.cpp the new timestep method was built, test that instead.
	//}	

	std::tuple<std::map<int, Vector2d>, std::map<int, double>> maps_new = g.calculateCellValues(particles);
	std::map<int, Vector2d> meanCellVelocities_new = std::get<0>(maps_new);
	std::map<int, double> rotationAngles_new = std::get<1>(maps_new);

	ASSERT_EQ(meanCellVelocities_new.size(), rotationAngles_new.size());

	std::filesystem::path cwd = std::filesystem::current_path();

	Out out(cwd.string() + l_data);
	out.writeToOut(origMeanVelocities, "before_collision.csv", "cell_vx_b,cell_vy_b");
	out.writeToOut(meanCellVelocities_new, "after_collision.csv", "cell_vx_a,cell_vy_a");
	*/
}