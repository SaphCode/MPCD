#include "pch.h"
#include "Constants.h"
#include <Eigen/Dense>
#include "Xoshiro.h"
#include "seeder.h"
#include "Particle.h"
#include <filesystem>
#include <fstream>
#include "Grid.h"

using namespace MPCD;
using namespace Eigen;

const int min_per_cell = MPCD::Constants::min_particles_per_cell;

TEST(Grid, Constants) {
	Grid g(min_per_cell);
	int num_hypothetical_x_cells = std::round(MPCD::Constants::Pipe::width / g.getCellDim());
	int num_hypothetical_y_cells = std::round(MPCD::Constants::Pipe::height / g.getCellDim());

	if (MPCD::Constants::Pipe::width > MPCD::Constants::Pipe::height) {
		ASSERT_GT(num_hypothetical_x_cells, num_hypothetical_y_cells);
	}
	else {
		ASSERT_LE(num_hypothetical_x_cells, num_hypothetical_y_cells);
	}

	int num_hypothetical_cells = num_hypothetical_x_cells * num_hypothetical_y_cells;

	double epsilon = 0.05;

	ASSERT_LT(std::abs(num_hypothetical_cells - g.getWantedNumCells()) / g.getWantedNumCells(), epsilon);

	ASSERT_GT(g.getWantedNumCells(), g.getMinNumCells());

	std::filesystem::path cwd = std::filesystem::current_path();

	std::ofstream outFile(cwd.string() + "//Data//constants.csv");
	outFile << "cell_dim,pipe_width,pipe_height" << "\n"; // header columns
	outFile << g.getCellDim() << "," << MPCD::Constants::Pipe::width << "," << MPCD::Constants::Pipe::height << std::endl;
}

TEST(Grid, MinNumberPerCell) {
	Grid g(min_per_cell);
	const double time_step = 1.0;
	const double aspect_ratio = MPCD::Constants::Pipe::width / MPCD::Constants::Pipe::height;
	const double max_x_position = MPCD::Constants::Pipe::width;
	const double max_y_position = MPCD::Constants::Pipe::height;
	const double max_x_velocity = std::max(MPCD::Constants::Pipe::width, MPCD::Constants::Pipe::height) / 100.0; // also in RandomGenerator.cpp
	const double max_y_velocity = max_x_velocity / aspect_ratio;

	Xoshiro xs_xpos(0.0, max_x_position);
	Xoshiro xs_ypos(0.0, max_y_position);
	Xoshiro xs_xvel(-max_x_velocity, max_x_velocity);
	Xoshiro xs_yvel(-max_y_velocity, max_y_velocity);

	const int rows = std::ceil(MPCD::Constants::Pipe::height / g.getCellDim());
	const int cols = std::ceil(MPCD::Constants::Pipe::width / g.getCellDim());
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
		
		Vector2i cell_index = p.shift(zero, g.getCellDim());

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
			if (*col < g.getMinParticlesPerCell()) {
				num_smaller++;
			}
		}
		row_ind++;
	}

	outFile.close();
	
	double epsilon = 0.01;

	EXPECT_LT(num_smaller / g.getWantedNumCells(), epsilon);

	/*
	Test that every cell has about 5 particles
	*/
}