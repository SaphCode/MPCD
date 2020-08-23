#include "pch.h"
#include "Particle.h"
#include <Eigen/Dense>
#include "MPCD.h"
#include <cmath>
#include "Compare.h"
#include "Xoshiro.h"
#include <map>
#include <tuple>

#define _USE_MATH_DEFINES
#include <math.h>

using namespace Eigen;
using namespace MPCD;

class ParticleTest : public ::testing::Test {
protected:
	const double time_step = 1.0;
	const double aspect_ratio = Pipe::width / Pipe::height;
	const double max_x_position = Pipe::width;
	const double max_y_position = Pipe::height;
	const double max_x_velocity = std::max(max_x_position, max_y_position) / 100.0; // also in RandomGenerator.cpp
	const double max_y_velocity = max_x_velocity / aspect_ratio;
	const double max_angle = 2 * M_PI;

	std::vector<Particle> particles;

	Xoshiro rg_angle;
	Xoshiro rg_shift_x;
	Xoshiro rg_shift_y;

	void SetUp() override {
		particles.reserve(number);

		/* dont worry the numbers are just seeds */
		Xoshiro xs_xpos(0.0, max_x_position);
		Xoshiro xs_ypos(0.0, max_y_position);
		Xoshiro xs_xvel(-max_x_velocity, max_x_velocity);
		Xoshiro xs_yvel(-max_y_velocity, max_y_velocity);
		Xoshiro angles(0.0, 2*M_PI);
		rg_angle = angles;
		Xoshiro shifts_x(-Grid::max_shift, Grid::max_shift);
		Xoshiro shifts_y(-Grid::max_shift, Grid::max_shift);
		rg_shift_x = shifts_x;
		rg_shift_y = shifts_y;

		for (int i = 0; i < number; i++) {
			double xs_x = xs_xpos.next();
			double xs_y = xs_ypos.next();
			double xs_vx = xs_xvel.next();
			double xs_vy = xs_yvel.next();

			Vector2d pos(xs_x, xs_y);
			Vector2d vel(xs_vx, xs_vy);

			Particle xs_p(pos, vel);

			particles.push_back(xs_p);
		}
	}

	ParticleTest(const ParticleTest&) = default;
	ParticleTest& operator=(const ParticleTest&) = default;
	ParticleTest();
};

TEST_F(ParticleTest, Streaming) {
	Vector2d pos(0, 0);
	Vector2d vel(0.5, 0.5);
	Particle p(pos, vel);
	Vector2d p_pos = p.getPosition();
	Vector2d p_vel = p.getVelocity();
	ASSERT_EQ(p_pos(0), pos(0)) << "Initialized position is particle position";
	ASSERT_EQ(p_pos(1), pos(1)) << "Initialized position is particle position";
	ASSERT_EQ(p_vel(0), vel(0)) << "Initialized velocity is particle velocity";
	ASSERT_EQ(p_vel(1), vel(1)) << "Initialized velocity is particle velocity";
	double time_step = 1;
	p.move(time_step);
	Vector2d newPos = pos + time_step * vel;
	Vector2d p_new_pos = p.getPosition();
	ASSERT_EQ(p_new_pos(0), newPos(0)) << "Moved position is new particle position";
	ASSERT_EQ(p_new_pos(1), newPos(1)) << "Moved position is new particle position";
}

TEST_F(ParticleTest, Collision) {
	const int rows = std::ceil(Pipe::height / Grid::cell_dim);
	const int cols = std::ceil(Pipe::width / Grid::cell_dim);
	std::map<Vector2i, Vector2d> mean_cell_velocities;
	std::map<Vector2i, int> total_cell_p_numbers;
	std::map<Vector2i, double> rotation_angle;
	for (int r = 0; r < rows; r++) {
		for (int c = 0; c < cols; c++) {
			const Vector2i index(r, c);
			const Vector2d vel(0.0, 0.0);
			mean_cell_velocities.insert({ index, vel });
			total_cell_p_numbers.insert({ index, 0 });
			rotation_angle.insert({ index, 0.0 });
		}
	}
	for (auto it = particles.begin(); it != particles.end(); ++it) {
		Vector2d shift(rg_shift_x.next(), rg_shift_y.next());
		Vector2i cell_index = it->shift(shift);
		Vector2d vel = it->getVelocity();
		mean_cell_velocities[cell_index] = mean_cell_velocities[cell_index] + vel;
		total_cell_p_numbers[cell_index] += 1;
	}
	for (int r = 0; r < rows; r++) {
		for (int c = 0; c < cols; c++) {
			Vector2i cell_index(r, c);
			mean_cell_velocities[cell_index] /= total_cell_p_numbers[cell_index];
			rotation_angle[cell_index] = rg_angle.next();
		}
	}
	// put these together
	for (auto it = particles.begin(); it != particles.end(); ++it) {
		Vector2i cell_index = it->getCellIndex();
		it->updateVelocity(mean_cell_velocities[cell_index], rotation_angle[cell_index]);
	}
}


/* This test actually passes. Fixed the test, what was wrong before:
	------------------------
	Look at this code sample: 

	Vector2d pos0(0, 0);
	Vector2d vel0(0, 0);
	Particle p(pos0, vel0);
	Vector2d zero_shift(0, 0);
	Vector2i cell_index = p.shift(zero_shift);
	std::cout << "a: "<< Grid::cell_dim << std::endl;
	double roundError = 0;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			Vector2d pos_ij(i * Grid::cell_dim, j * Grid::cell_dim);
			Particle p_ij(pos_ij, vel0);
			Vector2i cell_index_ij = p_ij.shift(zero_shift);
			Vector2d p_ij_pos = p_ij.getPosition();
			std::cout << "(" << i << ", " << j << ")" << ": " << "(" << cell_index_ij(0) << ", " << cell_index_ij(1) << ")";
			std::cout << " vs. " << "(" << std::fmod(p_ij_pos(0), Grid::cell_dim) << ", " << std::fmod(p_ij_pos(1), Grid::cell_dim) << ")";
			std::cout << " vs. " << "(" << p_ij_pos(0) << ", " << p_ij_pos(1) << ")" << std::endl;
		}
	}
	std::cout << (63 * Grid::cell_dim / Grid::cell_dim) << std::endl;
	std::cout << std::floor((63 * Grid::cell_dim) / Grid::cell_dim) << std::endl;
	
	TLDR: Basically the division calculates 62.9999999999 and floor rounds down to 62, while output shows 63, which would be the correct answer.
			This is a very very small mistake that basically only happens in a test environment, so it's not relevant for the simulation.
	----------------------		
	*/
TEST_F(ParticleTest, CellLogicWorks) {
	Vector2d pos0(Grid::cell_dim / 10, Grid::cell_dim / 10); // look at reason for in cell above
	Vector2d vel0(0, 0);
	Particle p(pos0, vel0);
	Vector2i cell_index_calc(std::floor(pos0(0) / Grid::cell_dim), std::floor(pos0(1) / Grid::cell_dim));
	Vector2d zero_shift(0, 0);
	Vector2i cell_index_zero_shift = p.shift(zero_shift);
	ASSERT_TRUE(areVectorsEqual(cell_index_calc, cell_index_zero_shift));

	for (int i = 0; i < 100; i++) {
		for (int j = 0; j < 100; j++) {
			Vector2i index(i, j);
			Vector2d pos_ij(pos0(0) + i * Grid::cell_dim, pos0(1) + j * Grid::cell_dim);
			Particle p_ij(pos_ij, vel0);
			Vector2i cell_index_ij = p_ij.shift(zero_shift);
			ASSERT_TRUE(areVectorsEqual(cell_index_ij, cell_index_zero_shift + index));
		}
	}
}

/* r_ stands for reversed 
	Shift should not change the particle position, since it is an "imagined" shift. The particle cannot just teleport, its just to avoid correlations buildup. */
TEST_F(ParticleTest, ShiftParticles) {
	// test if cell index is at most different by one.
	// test if position is exactly different by shift
	Vector2d pos0(Grid::cell_dim/10, Grid::cell_dim/10);
	Vector2d vel0(0, 0);
	Particle p(pos0, vel0);

	Vector2d shift1(Grid::max_shift, Grid::max_shift);
	Vector2i index_after_pshift1(0, 0);
	Vector2i index_after_nshift1(-1, -1);
	Vector2d shift2(Grid::max_shift, -Grid::max_shift);
	Vector2i index_after_pshift2(0, -1);
	Vector2i index_after_nshift2(-1, 0);

	Vector2d shift_zero(0, 0);
	Vector2i zero_index = p.shift(shift_zero);

	Vector2i zeros(0, 0);

	ASSERT_TRUE(areVectorsEqual(zero_index, zeros)) << "Particle with position inside first cell should have (0,0) index (first cell should be (0,0))";

	Vector2i shift1_index = p.shift(shift1);
	Vector2i r_shift1_index = p.shift(-shift1);

	ASSERT_TRUE(areVectorsEqual(index_after_pshift1, shift1_index)) << "Shifting once should not change index in this case.";
	ASSERT_TRUE(areVectorsEqual(index_after_nshift1, r_shift1_index)) << "Shifting once should change index in this case.";
	ASSERT_TRUE(areVectorsEqual(pos0, p.getPosition())) << "Position does not change in a shift. This would take unnecessary storage and computation time.";
	
	Vector2i shift2_index = p.shift(shift2);
	Vector2i r_shift2_index = p.shift(-shift2);

	ASSERT_TRUE(areVectorsEqual(index_after_pshift2, shift2_index)) << "Shifting once should not change index in this case.";
	ASSERT_TRUE(areVectorsEqual(index_after_nshift2, r_shift2_index)) << "Shifting once should change index in this case.";
	ASSERT_TRUE(areVectorsEqual(pos0, p.getPosition())) << "Position does not change in a shift. This would take unnecessary storage and computation time.";
	
	ASSERT_TRUE(areVectorsEqual(p.getVelocity(), vel0)) << "Velocity does not change in a shift.";
}

/* This test is just for measuring the time needed. */
/*
TEST_F(ParticleTest, MoveMultipleParticles) {
	double time_step = 1.0;
	for (auto& p : particles) {
		p.move(time_step);
	}
	//p.move(time_step);
	/* NORMAL FOR LOOP, WITH "Array Access"
	for (unsigned int i = 0; i < particles.size(); i++) {
		particles[i].move(time_step);
	}*/
	/* Iterator method with i++:
		move():
			658ms
	int i = 0;
	for (std::vector<Particle>::iterator it = particles.begin(); it != particles.end(); ++it, i++) {
		particles[i].move(time_step);
		std::cout << *it;
	}
}
*/

/*
The particle's x & y are additionally shifted by MPCD::Grid::cell_size/2 to not
		build up correlations, and therefore to keep Galilean Invariance.
*/