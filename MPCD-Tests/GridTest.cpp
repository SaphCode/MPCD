#include "pch.h"
#include "MPCD.h"
#include <Eigen/Dense>
#include "Xoshiro.h"
#include "seeder.h"
#include "Particle.h"
#include <filesystem>
#include <fstream>

using namespace MPCD;
using namespace Eigen;

TEST(Grid, Constants) {
	int num_hypothetical_x_cells = std::round(Pipe::width / Grid::cell_dim);
	int num_hypothetical_y_cells = std::round(Pipe::height / Grid::cell_dim);

	if (Pipe::width > Pipe::height) {
		ASSERT_GT(num_hypothetical_x_cells, num_hypothetical_y_cells);
	}
	else {
		ASSERT_LE(num_hypothetical_x_cells, num_hypothetical_y_cells);
	}

	int num_hypothetical_cells = num_hypothetical_x_cells * num_hypothetical_y_cells;

	double epsilon = 0.05;

	ASSERT_LT(std::abs(num_hypothetical_cells - Grid::wanted_num_cells) / Grid::wanted_num_cells, epsilon);

	ASSERT_GT(Grid::wanted_num_cells, Grid::min_num_cells);

	std::filesystem::path cwd = std::filesystem::current_path();

	std::ofstream outFile(cwd.string() + "//Data//constants.csv");
	outFile << "cell_dim,pipe_width,pipe_height" << "\n"; // header columns
	outFile << Grid::cell_dim << "," << Pipe::width << "," << Pipe::height << std::endl;
}

TEST(Grid, MinNumberPerCell) {
	const double time_step = 1.0;
	const double aspect_ratio = Pipe::width / Pipe::height;
	const double max_x_position = Pipe::width;
	const double max_y_position = Pipe::height;
	const double max_x_velocity = std::max(Pipe::width, Pipe::height) / 100.0; // also in RandomGenerator.cpp
	const double max_y_velocity = max_x_velocity / aspect_ratio;

	Xoshiro xs_xpos(0.0, max_x_position);
	Xoshiro xs_ypos(0.0, max_y_position);
	Xoshiro xs_xvel(-max_x_velocity, max_x_velocity);
	Xoshiro xs_yvel(-max_y_velocity, max_y_velocity);

	const int rows = std::ceil(Pipe::height / Grid::cell_dim);
	const int cols = std::ceil(Pipe::width / Grid::cell_dim);
	std::vector<std::vector<int>> frequencies;
	frequencies.resize(rows, std::vector<int>(cols, 0));

	const Vector2d zero(0, 0);

	for (int i = 0; i < number; i++) {
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
			if (*col < Grid::min_particles_per_cell) {
				num_smaller++;
			}
		}
		row_ind++;
	}

	outFile.close();
	
	double epsilon = 0.01;

	EXPECT_LT(num_smaller / Grid::wanted_num_cells, epsilon);

	/*
	Test that every cell has about 5 particles
	*/
}

TEST(Grid, MeanVelocity) {
	Vector2d pos(0, 0);
	double xvel = 1;
	double yvel = 1;
	Vector2d vel(xvel, yvel);
	std::vector<Particle> particles;
	int num = 5;
	particles.reserve(num);

	for (int i = 0; i < num; i++) {
		Particle p(pos, vel);
		particles.push_back(p);
	}
	Grid g(5);

	g.calculateCellValues(particles);
	ASSERT_EQ(meanCellVelocity[0], xvel);
	ASSERT_EQ(meanCellVelocity[1], yvel);

	double mxvel = -xvel;
	double myvel = -yvel;
	Vector2d mvel(mxvel, myvel);

	for (int i = 0; i < num; i++) {
		Particle p(pos, mvel);
		particles.push_back(p);
	}

	Grid g2(10);
	g2.calculateCellValues();
	
	ASSERT_EQ(meanCellVelocity[0], num*(xvel + mxvel)/(2*num));
	ASSERT_EQ(meanCellVelocity[1], num * (yvel + myvel) / (2 * num));
}