#include "pch.h"
#include "Constants.h"
#include <Eigen/Dense>
#include "Xoshiro.h"
#include "seeder.h"
#include "Particle.h"
#include <filesystem>
#include <fstream>
#include "Grid.h"
#include "Simulation.h"
#include "Compare.h"

using namespace MPCD;
using namespace Eigen;

const int average_particles_per_cell = MPCD::Constants::Grid::average_particles_per_cell;

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

	int num_cells = MPCD::Constants::Grid::num_cells;
	ASSERT_LT(std::abs(num_hypothetical_cells - num_cells) / num_cells, epsilon);

	std::filesystem::path cwd = std::filesystem::current_path();

	std::string filename("//Data//constants_");
	std::string num(std::to_string(average_particles_per_cell));
	std::string csv(".csv");

	std::ofstream outFile(cwd.string() + filename + num + csv);
	outFile << "timesteps,time_lapse,cell_dim,x_0,y_0,x_max,y_max,average_particles_per_cell,total_number_of_particles" << "\n"; // header columns
	outFile << MPCD::Constants::timesteps << "," << MPCD::Constants::time_lapse << "," << cell_dim << "," << MPCD::Constants::Pipe::x_0 << "," << MPCD::Constants::Pipe::y_0 << "," <<
		MPCD::Constants::Pipe::x_max << "," << MPCD::Constants::Pipe::y_max << "," << MPCD::Constants::Grid::average_particles_per_cell << "," << MPCD::Constants::number << std::endl;
	outFile.close();
}

TEST(Grid, GridOffsetRight) {
	Vector2d pos_lower_left_corner(MPCD::Constants::Pipe::x_0, MPCD::Constants::Pipe::y_0);
	Vector2d dummy_vel(0, 0);
	Particle p_lower_left_corner(pos_lower_left_corner, dummy_vel);
	Vector2i index_lower_left_corner = p_lower_left_corner.getCellIndex();
	Vector2i supposedIndex_lower_left_corner(0, 0);
	areVectorsEqual(index_lower_left_corner, supposedIndex_lower_left_corner);

	Vector2d pos_lower_right_corner(MPCD::Constants::Pipe::x_max, MPCD::Constants::Pipe::y_0);
	Particle p_lower_right_corner(pos_lower_right_corner, dummy_vel);
	Vector2i index_lower_right_corner = p_lower_right_corner.getCellIndex();
	Vector2i supposedIndex_lower_right_corner(0, MPCD::Constants::Grid::max_cols);
	areVectorsEqual(index_lower_right_corner, supposedIndex_lower_right_corner);

	Vector2d pos_upper_left_corner(MPCD::Constants::Pipe::x_0, MPCD::Constants::Pipe::y_max);
	Particle p_upper_left_corner(pos_upper_left_corner, dummy_vel);
	Vector2i index_upper_left_corner = p_upper_left_corner.getCellIndex();
	Vector2i supposedIndex_upper_left_corner(MPCD::Constants::Grid::max_rows, 0);
	areVectorsEqual(index_upper_left_corner, supposedIndex_upper_left_corner);

	Vector2d pos_upper_right_corner(MPCD::Constants::Pipe::x_max, MPCD::Constants::Pipe::y_max);
	Particle p_upper_right_corner(pos_upper_right_corner, dummy_vel);
	Vector2i index_upper_right_corner = p_upper_right_corner.getCellIndex();
	Vector2i supposedIndex_upper_right_corner(MPCD::Constants::Grid::max_rows, 0);
	areVectorsEqual(index_upper_right_corner, supposedIndex_upper_right_corner);
}

TEST(Grid, AverageNumberPerCell) {
	const double cell_dim = MPCD::Constants::Grid::cell_dim;
	const double time_step = MPCD::Constants::time_lapse;
	const double aspect_ratio = MPCD::Constants::Pipe::width / MPCD::Constants::Pipe::height;
	const double min_x_position = MPCD::Constants::Pipe::x_0;
	const double max_x_position = MPCD::Constants::Pipe::x_max;
	const double min_y_position = MPCD::Constants::Pipe::y_0;
	const double max_y_position = MPCD::Constants::Pipe::y_max;
	const double max_x_velocity = std::max(MPCD::Constants::Pipe::width, MPCD::Constants::Pipe::height) / 100.0; // also in RandomGenerator.cpp
	const double max_y_velocity = max_x_velocity / aspect_ratio;

	Xoshiro xs_xpos(min_x_position, max_x_position);
	Xoshiro xs_ypos(min_y_position, max_y_position);
	Xoshiro xs_xvel(-max_x_velocity, max_x_velocity);
	Xoshiro xs_yvel(-max_y_velocity, max_y_velocity);

	const int rows = MPCD::Constants::Grid::max_rows;
	const int cols = MPCD::Constants::Grid::max_cols;
	std::map<int, int> frequencies;

	const Vector2d zero(0, 0);

	for (int i = 0; i < MPCD::Constants::number; i++) {
		double xs_x = xs_xpos.next();
		double xs_y = xs_ypos.next();
		double xs_vx = xs_xvel.next();
		double xs_vy = xs_yvel.next();

		Vector2d pos(xs_x, xs_y);
		Vector2d vel(xs_vx, xs_vy);
		Particle p(pos, vel);

		Vector2i cell_index = p.getCellIndex();
		int linearIndex = Grid::convertToLinearIndex(cell_index);

		ASSERT_LE(cell_index(0), rows);
		ASSERT_LE(cell_index(1), cols);

		frequencies[linearIndex] += 1;

		//particles.push_back(p);
	}

	int num_smaller = 0;

	std::filesystem::path cwd = std::filesystem::current_path();

	std::string filename("//Data//cell_frequencies");
	std::string num("_av_" + std::to_string(average_particles_per_cell));
	std::string csv(".csv");

	std::ofstream outFile(cwd.string() + filename + num + csv);
	outFile << "i,j,n" << "\n"; // header columns

	int row_ind = 0;
	int col_ind = 0;
	int total = 0;
	for (auto const& [key, value] : frequencies) {
		Vector2i index = Grid::convertToIndex(key);
		outFile << index[0] << "," << index[1] << "," << value << "\n";
		if (value < average_particles_per_cell) {
			num_smaller++;
		}
		total += value;
	}
	outFile.close();
	
	double threshhold = 0.55; // THIS TEST IS BROKEN

	ASSERT_LT((double)num_smaller / MPCD::Constants::Grid::num_cells, threshhold);
	ASSERT_GE(total / MPCD::Constants::Grid::num_cells, average_particles_per_cell);

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

	ASSERT_EQ(A(0, 0), 0);
	ASSERT_EQ(A(0, 1), 1);
	ASSERT_EQ(A(0, 2), 2);
	ASSERT_EQ(A(1, 0), 3);
	ASSERT_EQ(A(1, 1), 4);
	ASSERT_EQ(A(1, 2), 5);
	ASSERT_EQ(A(2, 0), 6);
	ASSERT_EQ(A(2, 1), 7);
	ASSERT_EQ(A(2, 2), 8);
	
	int cols = 3;

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
	Vector2d startPos(MPCD::Constants::Grid::grid_x_shift + cell_dim / 10, MPCD::Constants::Grid::grid_y_shift + cell_dim / 10);

	const int rows = MPCD::Constants::Grid::max_rows;
	const int cols = MPCD::Constants::Grid::max_cols;

	Xoshiro rg_angle(0.0, 2 * 3.141); // not important in this test

	/*
	* This block should place 1 particle in each cell. The mean of every cell should then be the velocity of the particle itself.
	*/
	Particle init(startPos, vel);
	Eigen::Vector2i startIndex = init.getCellIndex();
	for (int i = startIndex[0]; i < rows; i++) {
		int linearIndex_before = Grid::convertToLinearIndex(startIndex) - 1 + i * cols;
		for (int j = startIndex[1]; j < cols; j++) {
			Vector2d offset(j * cell_dim, i * cell_dim);
			Vector2d pos = startPos + offset;
			Particle p(pos, vel);
			particles.push_back(p);

			int linearIndex = MPCD::Grid::convertToLinearIndex(p.getCellIndex());
			ASSERT_EQ(p.getCellIndex()[0], i);
			ASSERT_EQ(p.getCellIndex()[1], j);
			Vector2i index(i, j);
			ASSERT_EQ(linearIndex, MPCD::Grid::convertToLinearIndex(index));
			ASSERT_EQ(linearIndex, linearIndex_before + 1);
			linearIndex_before = linearIndex;
		}
	}
}



TEST(Grid, MeanVelocity_OneParticlePerCell) {
	double xvel = 1;
	double yvel = 1;
	Vector2d startVel(xvel, yvel);
	Vector2d zero(0, 0);

	const int rows = MPCD::Constants::Grid::max_rows;
	const int cols = MPCD::Constants::Grid::max_cols;

	std::vector<Particle> particles;
	int num = MPCD::Constants::number;
	particles.reserve((double)rows * cols);

	double cell_dim = MPCD::Constants::Grid::cell_dim;
	Vector2d startPos(MPCD::Constants::Grid::grid_x_shift + cell_dim / 10, MPCD::Constants::Grid::grid_y_shift + cell_dim / 10);
	std::map<int, Eigen::Vector2d> totalCellVelocities;
	std::map<int, int> cellParticles;

	Xoshiro rg_angle(0.0, 2 * 3.141); // not important in this test

	/*
	* This block should place 1 particle in each cell. The mean of every cell should then be the velocity of the particle itself.
	*/
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			Vector2i index(i, j);
			int linearIndex = Grid::convertToLinearIndex(index);
			totalCellVelocities.insert(std::make_pair(linearIndex, zero));
		}
	}
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			Vector2d offset(j * cell_dim, i * cell_dim);
			Vector2d pos = startPos + offset;
			Particle p(pos, startVel);
			particles.push_back(p);
			
			int linearIndex = MPCD::Grid::convertToLinearIndex(p.getCellIndex());

			Vector2d vel = p.getVelocity();
			
			totalCellVelocities[linearIndex] = totalCellVelocities[linearIndex] + vel;
			cellParticles[linearIndex] += 1;
		}
	}

	//updateVelocity(particles, meanCellVelocities, cellRotationAngles);
	Simulation sim(particles);
	std::tuple<std::map<int, Eigen::Vector2d>, std::map<int, double>> cellMeanVelocityAndRotationAngle = sim.calculateCellQuantities(totalCellVelocities, cellParticles);
	std::map<int, Eigen::Vector2d> meanCellVelocities = std::get<0>(cellMeanVelocityAndRotationAngle);
	std::map<int, double> cellRotationAngles = std::get<1>(cellMeanVelocityAndRotationAngle);
	
	for (auto const& [key, val] : meanCellVelocities) {
		double epsilon = 0.0001;
		ASSERT_LT(std::abs(val[0] - startVel[0]) / startVel[0], epsilon);
		ASSERT_LT(std::abs(val[1] - startVel[1]) / startVel[1], epsilon);
	}

}

TEST(Grid, MeanVelocity_AllParticlesInOneCell) {
	std::vector<Particle> particles;
	const int num = MPCD::Constants::number;
	particles.reserve(num);
	const double cell_dim = MPCD::Constants::Grid::cell_dim;

	const int max_rows = MPCD::Constants::Grid::max_rows;
	const int max_cols = MPCD::Constants::Grid::max_cols;

	Xoshiro rg_angle(0, 2 * 3.141); // not importatnt here

	Vector2d pos(3.1 * cell_dim, 3.1 * cell_dim);
	double vx = 1.1111;
	double vy = -12;
	Vector2d vel(vx, vy);
	Vector2d zero(0, 0);

	std::map<int, Eigen::Vector2d> totalCellVelocities;
	std::map<int, int> cellParticles;


	for (int i = 0; i < max_rows; i++) {
		for (int j = 0; j < max_cols; j++) {
			Vector2i index(i, j);
			int linearIndex = Grid::convertToLinearIndex(index);
			totalCellVelocities.insert(std::make_pair(linearIndex, zero));
		}
	}

	for (int i = 0; i < num; i++) {
		Particle p(pos, vel);
		particles.push_back(p);
		int linearIndex = MPCD::Grid::convertToLinearIndex(p.getCellIndex());
		totalCellVelocities[linearIndex] += vel;
		cellParticles[linearIndex] += 1;
	}

	Simulation sim(particles);
	std::tuple<std::map<int, Eigen::Vector2d>, std::map<int, double>> cellMeanVelocityAndRotationAngle = sim.calculateCellQuantities(totalCellVelocities, cellParticles);
	std::map<int, Eigen::Vector2d> meanCellVelocities = std::get<0>(cellMeanVelocityAndRotationAngle);
	std::map<int, double> cellRotationAngles = std::get<1>(cellMeanVelocityAndRotationAngle);

	for (auto const& [key, val] : meanCellVelocities) {
		if (val != zero) {
			double epsilon = 0.0001;
			ASSERT_LT(std::abs((val[0] - vx) / vx), epsilon);
			ASSERT_LT(std::abs((val[1] - vy) / vy), epsilon);
		}
	}
}

TEST(Grid, ConvertToIndex_3x3) {
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

	Eigen::Vector2i cindex_00 = Grid::convertToIndex(A(0, 0), cols);
	Eigen::Vector2i cindex_01 = Grid::convertToIndex(A(0, 1), cols);
	Eigen::Vector2i cindex_02 = Grid::convertToIndex(A(0, 2), cols);
	Eigen::Vector2i cindex_10 = Grid::convertToIndex(A(1, 0), cols);
	Eigen::Vector2i cindex_11 = Grid::convertToIndex(A(1, 1), cols);
	Eigen::Vector2i cindex_12 = Grid::convertToIndex(A(1, 2), cols);
	Eigen::Vector2i cindex_20 = Grid::convertToIndex(A(2, 0), cols);
	Eigen::Vector2i cindex_21 = Grid::convertToIndex(A(2, 1), cols);
	Eigen::Vector2i cindex_22 = Grid::convertToIndex(A(2, 2), cols);

	ASSERT_EQ(cindex_00, index_00);
	ASSERT_EQ(cindex_01, index_01);
	ASSERT_EQ(cindex_02, index_02);
	ASSERT_EQ(cindex_10, index_10);
	ASSERT_EQ(cindex_11, index_11);
	ASSERT_EQ(cindex_12, index_12);
	ASSERT_EQ(cindex_20, index_20);
	ASSERT_EQ(cindex_21, index_21);
	ASSERT_EQ(cindex_22, index_22);

}

TEST(Grid, ConvertToIndex_WholeGrid) {
	const double xvel = 1;
	const double yvel = 1;
	Vector2d vel(xvel, yvel);

	std::vector<Particle> particles;
	const int num = MPCD::Constants::number;
	particles.reserve(num);

	const double cell_dim = MPCD::Constants::Grid::cell_dim;
	Vector2d startPos(MPCD::Constants::Grid::grid_x_shift + cell_dim / 10, MPCD::Constants::Grid::grid_y_shift + cell_dim / 10);

	const int rows = MPCD::Constants::Grid::max_rows;
	const int cols = MPCD::Constants::Grid::max_cols;

	Xoshiro rg_angle(0.0, 2 * 3.141); // not important in this test

	/*
	* This block should place 1 particle in each cell. The mean of every cell should then be the velocity of the particle itself.
	*/
	Particle init(startPos, vel);
	Eigen::Vector2i startIndex = init.getCellIndex();
	for (int i = startIndex[0]; i < rows; i++) {
		int linearIndex_before = Grid::convertToLinearIndex(startIndex) - 1 + i * cols;
		for (int j = startIndex[1]; j < cols; j++) {
			Vector2d offset(j * cell_dim, i * cell_dim);
			Vector2d pos = startPos + offset;
			Particle p(pos, vel);
			particles.push_back(p);

			Vector2i pindex = p.getCellIndex();
			int linearIndex = MPCD::Grid::convertToLinearIndex(p.getCellIndex());
			Vector2i cindex = MPCD::Grid::convertToIndex(linearIndex);
			ASSERT_EQ(pindex[0], cindex[0]);
			ASSERT_EQ(pindex[1], cindex[1]);
			ASSERT_EQ(p.getCellIndex()[0], i);
			ASSERT_EQ(p.getCellIndex()[1], j);
			Vector2i index(i, j);
			ASSERT_EQ(linearIndex, MPCD::Grid::convertToLinearIndex(index));
			ASSERT_EQ(linearIndex, linearIndex_before + 1);
			linearIndex_before = linearIndex;
		}
	}
}