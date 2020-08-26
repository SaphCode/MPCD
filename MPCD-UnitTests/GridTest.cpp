#include "pch.h"
#include "Constants.h"
#include <Eigen/Dense>
#include "Xoshiro.h"
#include "seeder.h"
#include "Particle.h"
#include <filesystem>
#include <fstream>
#include "Grid.h"
#include "MPCD.h"

using namespace MPCD;
using namespace Eigen;

const int min_per_cell = MPCD::Constants::Grid::min_particles_per_cell;

TEST(Grid, Constants) {
	double cell_dim = MPCD::Constants::Grid::cell_dim;
	int num_hypothetical_x_cells = std::round(MPCD::Constants::Pipe::width / cell_dim);
	int num_hypothetical_y_cells = std::round(MPCD::Constants::Pipe::height / cell_dim);

	if (MPCD::Constants::Pipe::width > MPCD::Constants::Pipe::height) {
		ASSERT_GT(num_hypothetical_x_cells, num_hypothetical_y_cells);
	}
	else {
		ASSERT_LE(num_hypothetical_x_cells, num_hypothetical_y_cells);
	}

	int num_hypothetical_cells = num_hypothetical_x_cells * num_hypothetical_y_cells;

	double epsilon = 0.05;

	int wanted_num_cells = MPCD::Constants::Grid::wanted_num_cells;
	int min_num_cells = MPCD::Constants::Grid::min_num_cells;
	ASSERT_LT(std::abs(num_hypothetical_cells - wanted_num_cells) / wanted_num_cells, epsilon);

	ASSERT_GT(wanted_num_cells, min_num_cells);

	std::filesystem::path cwd = std::filesystem::current_path();

	std::ofstream outFile(cwd.string() + "//Data//constants.csv");
	outFile << "cell_dim,pipe_width,pipe_height" << "\n"; // header columns
	outFile << cell_dim << "," << MPCD::Constants::Pipe::width << "," << MPCD::Constants::Pipe::height << std::endl;
}

TEST(Grid, MinNumberPerCell) {
	const double cell_dim = MPCD::Constants::Grid::cell_dim;
	const double time_step = MPCD::Constants::time_lapse;
	const double aspect_ratio = MPCD::Constants::Pipe::width / MPCD::Constants::Pipe::height;
	const double max_x_position = MPCD::Constants::Pipe::width;
	const double max_y_position = MPCD::Constants::Pipe::height;
	const double max_x_velocity = std::max(MPCD::Constants::Pipe::width, MPCD::Constants::Pipe::height) / 100.0; // also in RandomGenerator.cpp
	const double max_y_velocity = max_x_velocity / aspect_ratio;

	Xoshiro xs_xpos(0.0, max_x_position);
	Xoshiro xs_ypos(0.0, max_y_position);
	Xoshiro xs_xvel(-max_x_velocity, max_x_velocity);
	Xoshiro xs_yvel(-max_y_velocity, max_y_velocity);

	const int rows = std::ceil(MPCD::Constants::Pipe::height / cell_dim);
	const int cols = std::ceil(MPCD::Constants::Pipe::width / cell_dim);
	std::vector<std::vector<int>> frequencies;
	frequencies.resize(rows, std::vector<int>(cols, 0));

	const Vector2d zero(0, 0);

	for (int i = 0; i < MPCD::Constants::number; i++) {
		double xs_x = xs_xpos.next();
		double xs_y = xs_ypos.next();
		double xs_vx = xs_xvel.next();
		double xs_vy = xs_yvel.next();

		Vector2d pos(xs_x, xs_y);
		Vector2d vel(xs_vx, xs_vy);
		Particle p(pos, vel);
		
		Vector2i cell_index = p.shift(zero);

		ASSERT_LT(cell_index(0), rows);
		ASSERT_LT(cell_index(1), cols);

		frequencies[cell_index(0)][cell_index(1)] += 1;

		//particles.push_back(p);
	}

	int num_smaller = 0;

	std::filesystem::path cwd = std::filesystem::current_path();

	std::ofstream outFile(cwd.string() + "//Data//cell_frequencies.csv");
	outFile << "i,j,n" << "\n"; // header columns

	int row_ind = 0;
	int col_ind = 0;
	for (auto row = frequencies.begin(); row != frequencies.end(); ++row) {
		for (auto col = row->begin(); col != row->end(); ++col) {
			//EXPECT_GE(*col, Grid::min_particles_per_cell);
			outFile << row_ind << "," << col_ind << "," << *col << "\n";
			col_ind++;
			if (*col < MPCD::Constants::Grid::min_particles_per_cell) {
				num_smaller++;
			}
		}
		row_ind++;
	}

	outFile.close();
	
	double epsilon = 0.01;

	EXPECT_LT(num_smaller / MPCD::Constants::Grid::wanted_num_cells, epsilon);

	/*
	Test that every cell has about 5 particles
	*/
}

TEST(Grid, MeanVelocity) {
	
	double xvel = 1;
	double yvel = 1;
	Vector2d vel(xvel, yvel);

	std::vector<Particle> particles;
	int num = MPCD::Constants::number;
	particles.reserve(num);

	double cell_dim = MPCD::Constants::Grid::cell_dim;
	Vector2d startPos(cell_dim / 10, cell_dim / 10);
	std::map<int, Vector2d> totalCellVelocities;
	std::map<int, int> cellParticles;
	std::map<int, double> cellRotationAngles;

	Xoshiro rg_angle(0.0, 2 * 3.141); // not important in this test

	/*
	* This block should place 1 particle in each cell. The mean of every cell should then be the velocity of the particle itself.
	*/
	for (int i = 0; i < MPCD::Constants::Grid::rows; i++) {
		for (int j = 0; j < MPCD::Constants::Grid::cols; j++) {
			Vector2d offset(j * cell_dim, i * cell_dim);
			Vector2d pos = startPos + offset;
			Particle p(pos, vel);
			particles.push_back(p);
			
			int linearIndex = MPCD::Grid::convertToLinearIndex(p.getCellIndex(p.getPosition()));

			Vector2d vel = p.getVelocity();
			totalCellVelocities[linearIndex] += vel;
			cellParticles[linearIndex] += 1;
		}
	}

	MPCD::calculateCellQuantities(totalCellVelocities, cellParticles, cellRotationAngles, rg_angle);
	
	std::map<int, Vector2d> meanCellVelocities = totalCellVelocities;
	for (auto const& [key, val] : meanCellVelocities) {
		double epsilon = 0.0001;
		ASSERT_LT(std::abs(val[0] - xvel) / xvel, epsilon);
		ASSERT_LT(std::abs(val[1] - yvel) / yvel, epsilon);
	}

	/* ------------------------------------ */

	std::vector<Particle> particles2;
	particles2.reserve(num);

	Vector2d pos2(3.1 * cell_dim, 3.1 * cell_dim);
	double vx2 = 1.1111;
	double vy2 = -12;
	Vector2d vel2(vx2, vy2);

	std::map<int, Vector2d> totalCellVelocities2;
	std::map<int, int> cellParticles2;
	std::map<int, double> cellRotationAngles2;

	for (int i = 0; i < num; i++) {
		Particle p(pos2, vel2);
		particles2.push_back(p);
		int linearIndex = MPCD::Grid::convertToLinearIndex(p.getCellIndex(pos2));
		totalCellVelocities2[linearIndex] += vel2;
		cellParticles2[linearIndex] += 1;
	}

	calculateCellQuantities(totalCellVelocities2, cellParticles2, cellRotationAngles2, rg_angle);
	std::map<int, Vector2d> meanCellVelocities2 = totalCellVelocities2;

	for (auto const& [key, val] : meanCellVelocities2) {
		double epsilon = 0.0001;
		ASSERT_LT(std::abs(val[0] - vx2) / vx2, epsilon);
		ASSERT_LT(std::abs(val[1] - vy2) / vy2, epsilon);
	}
}