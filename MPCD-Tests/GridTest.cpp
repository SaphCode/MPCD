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

/* Converts 2d indexes into linear indexes.
	Example:
	A = [[1, 2, 3]
			[4, 5, 6]
			[7, 8, 9]]
	A[0,0], index = (0,0)
	returns 0
	A[1,1], index = (1,1)
	returns 4
	@param index needs to be 2d!
	@returns linear index
*/
TEST(Grid, ConvertToLinearIndex_3x3) {
	Eigen::Matrix3i A;
	A << 0, 1, 2,
		3, 4, 5,
		6, 7, 8;

	int cols = 3;
	
	ASSERT_EQ(A(0, 0), 0);
	ASSERT_EQ(A(0, 1), 1);
	ASSERT_EQ(A(0, 2), 2);
	ASSERT_EQ(A(1, 0), 3);
	ASSERT_EQ(A(1, 1), 4);
	ASSERT_EQ(A(1, 2), 5);
	ASSERT_EQ(A(2, 0), 6);
	ASSERT_EQ(A(2, 1), 7);
	ASSERT_EQ(A(2, 2), 8);

	Eigen::Vector2i index_00(0, 0);
	Eigen::Vector2i index_01(0, 1);
	Eigen::Vector2i index_02(0, 2);
	Eigen::Vector2i index_10(1, 0);
	Eigen::Vector2i index_11(1, 1);
	Eigen::Vector2i index_12(1, 2);
	Eigen::Vector2i index_20(2, 0);
	Eigen::Vector2i index_21(2, 1);
	Eigen::Vector2i index_22(2, 2);

	ASSERT_EQ(Grid::convertToLinearIndex(index_00, cols), A(0, 0));
	ASSERT_EQ(Grid::convertToLinearIndex(index_01, cols), A(0, 1));
	ASSERT_EQ(Grid::convertToLinearIndex(index_02, cols), A(0, 2));
	ASSERT_EQ(Grid::convertToLinearIndex(index_10, cols), A(1, 0));
	ASSERT_EQ(Grid::convertToLinearIndex(index_11, cols), A(1, 1));
	ASSERT_EQ(Grid::convertToLinearIndex(index_12, cols), A(1, 2));
	ASSERT_EQ(Grid::convertToLinearIndex(index_20, cols), A(2, 0));
	ASSERT_EQ(Grid::convertToLinearIndex(index_21, cols), A(2, 1));
	ASSERT_EQ(Grid::convertToLinearIndex(index_22, cols), A(2, 2));

}

TEST(Grid, ConvertToLinearIndex_WholeGrid) {
	const double xvel = 1;
	const double yvel = 1;
	Vector2d vel(xvel, yvel);

	std::vector<Particle> particles;
	const int num = MPCD::Constants::number;
	particles.reserve(num);

	const double cell_dim = MPCD::Constants::Grid::cell_dim;
	Vector2d startPos(cell_dim / 20, cell_dim / 20);
	std::map<int, double> totalCellVelocityX;
	std::map<int, double> totalCellVelocityY;
	std::map<int, int> cellParticles;

	const int rows = MPCD::Constants::Grid::rows;
	const int cols= MPCD::Constants::Grid::cols;

	Xoshiro rg_angle(0.0, 2 * 3.141); // not important in this test

	/*
	* This block should place 1 particle in each cell. The mean of every cell should then be the velocity of the particle itself.
	*/
	Particle init(startPos, vel);
	Eigen::Vector2i startIndex = init.getCellIndex(init.getPosition());
	for (int i = startIndex[0]; i < rows; i++) {
		int linearIndex_before = Grid::convertToLinearIndex(startIndex, cols) - 1 + i * cols;
		for (int j = startIndex[1]; j < MPCD::Constants::Grid::cols; j++) {
			Vector2d offset(j * cell_dim, i * cell_dim);
			Vector2d pos = startPos + offset;
			Particle p(pos, vel);
			particles.push_back(p);

			int linearIndex = MPCD::Grid::convertToLinearIndex(p.getCellIndex(p.getPosition()), cols);
			ASSERT_EQ(p.getCellIndex(p.getPosition())[0], i);
			ASSERT_EQ(p.getCellIndex(p.getPosition())[1], j);
			Vector2i index(i, j);
			ASSERT_EQ(linearIndex, MPCD::Grid::convertToLinearIndex(index, rows));
			ASSERT_EQ(linearIndex, linearIndex_before + 1);
			linearIndex_before = linearIndex;
		}
	}
}



TEST(Grid, MeanVelocity_OneParticlePerCell) {
	double xvel = 1;
	double yvel = 1;
	Vector2d vel(xvel, yvel);

	std::vector<Particle> particles;
	int num = MPCD::Constants::number;
	particles.reserve((double)MPCD::Constants::Grid::rows * MPCD::Constants::Grid::cols);

	double cell_dim = MPCD::Constants::Grid::cell_dim;
	Vector2d startPos(cell_dim / 10, cell_dim / 10);
	std::map<int, double> totalCellVelocityX;
	std::map<int, double> totalCellVelocityY;
	std::map<int, int> cellParticles;

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
			
			int linearIndex = MPCD::Grid::convertToLinearIndex(p.getCellIndex(p.getPosition()), MPCD::Constants::Grid::cols);

			Vector2d vel = p.getVelocity();
			totalCellVelocityX[linearIndex] += vel[0];
			totalCellVelocityY[linearIndex] += vel[1];
			cellParticles[linearIndex] += 1;
		}
	}

	//updateVelocity(particles, meanCellVelocities, cellRotationAngles);
	std::tuple<std::map<int, double>, std::map<int, double>, std::map<int, double>> cellMeanVelocityAndRotationAngle = MPCD::calculateCellQuantities(totalCellVelocityX, totalCellVelocityY, cellParticles, rg_angle);
	std::map<int, double> meanCellVelocityX = std::get<0>(cellMeanVelocityAndRotationAngle);
	std::map<int, double> meanCellVelocityY = std::get<1>(cellMeanVelocityAndRotationAngle);
	std::map<int, double> cellRotationAngles = std::get<2>(cellMeanVelocityAndRotationAngle);
	
	for (auto const& [key, val] : meanCellVelocityX) {
		double epsilon = 0.0001;
		ASSERT_LT(std::abs(meanCellVelocityX[key] - xvel) / xvel, epsilon);
		ASSERT_LT(std::abs(meanCellVelocityY[key] - yvel) / yvel, epsilon);
	}

}

TEST(Grid, MeanVelocity_AllParticlesInOneCell) {
	std::vector<Particle> particles;
	const int num = MPCD::Constants::number;
	particles.reserve(num);
	const double cell_dim = MPCD::Constants::Grid::cell_dim;
	const int cols = MPCD::Constants::Grid::cols;

	Xoshiro rg_angle(0, 2 * 3.141); // not importatnt here

	Vector2d pos(3.1 * cell_dim, 3.1 * cell_dim);
	double vx = 1.1111;
	double vy = -12;
	Vector2d vel(vx, vy);

	std::map<int, double> totalCellVelocityX;
	std::map<int, double> totalCellVelocityY;
	std::map<int, int> cellParticles;

	for (int i = 0; i < num; i++) {
		Particle p(pos, vel);
		particles.push_back(p);
		int linearIndex = MPCD::Grid::convertToLinearIndex(p.getCellIndex(pos), cols);
		totalCellVelocityX[linearIndex] += vel[0];
		totalCellVelocityY[linearIndex] += vel[1];
		cellParticles[linearIndex] += 1;
	}

	std::tuple<std::map<int, double>, std::map<int, double>, std::map<int, double>> cellMeanVelocityAndRotationAngle = MPCD::calculateCellQuantities(totalCellVelocityX, totalCellVelocityY, cellParticles, rg_angle);
	std::map<int, double> meanCellVelocityX = std::get<0>(cellMeanVelocityAndRotationAngle);
	std::map<int, double> meanCellVelocityY = std::get<0>(cellMeanVelocityAndRotationAngle);
	std::map<int, double> cellRotationAngles = std::get<1>(cellMeanVelocityAndRotationAngle);

	for (auto const& [key, val] : meanCellVelocityX) {
		double epsilon = 0.0001;
		ASSERT_LT(std::abs(meanCellVelocityX[key] - vx) / vx, epsilon);
		ASSERT_LT(std::abs(meanCellVelocityY[key] - vy) / vy, epsilon);
	}
}