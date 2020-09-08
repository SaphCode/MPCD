#include "pch.h"
#include "MersenneTwister.h"
#include "Xoshiro.h"
#include "Particle.h"
#include <filesystem>
#include <fstream>
#include "seeder.h"
#include "Constants.h"
#include "Out.h"
#include "Locations.h"

#define _USE_MATH_DEFINES
#include <math.h>

// file system stuff
//#ifdef WINDOWS
// THIS CODE WILL ONLY WORK ON WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
/*#else
#include <unistd.h>
#define GetCurrentDir getcwd*/

//#endif


using namespace Eigen;
using namespace MPCD;

class JupyterNotebookTests : public ::testing::Test {
protected:
	const double time_step = 1.0;
	const double aspect_ratio = MPCD::Constants::Pipe::width / MPCD::Constants::Pipe::height;
	const double max_x_position = MPCD::Constants::Pipe::width;
	const double max_y_position = MPCD::Constants::Pipe::height;
	const double max_x_velocity = std::max(max_x_position, max_y_position) / 100.0; // also in RandomGenerator.cpp
	const double max_y_velocity = max_x_velocity / aspect_ratio;
	const double max_angle = 2 * M_PI;

	std::vector<Particle> mers_particles;
	std::vector<double> mers_angles;

	std::vector<Particle> xs_particles;
	std::vector<double> xs_angles;

	void SetUp() override {
		mers_particles.reserve(MPCD::Constants::number);
		mers_angles.reserve(MPCD::Constants::number);

		xs_particles.reserve(MPCD::Constants::number);
		xs_angles.reserve(MPCD::Constants::number);

		MersenneTwister rg_xpos(0.0, max_x_position, DistributionType::UNIFORM);
		MersenneTwister rg_ypos(0.0, max_y_position, DistributionType::UNIFORM);
		MersenneTwister rg_xvel(-max_x_velocity, max_x_velocity, DistributionType::UNIFORM);
		MersenneTwister rg_yvel(-max_y_velocity, max_y_velocity, DistributionType::UNIFORM);
		MersenneTwister rg_angle(0.0, max_angle, DistributionType::UNIFORM);

		/* dont worry the numbers are just seeds */
		Xoshiro xs_xpos(0.0, max_x_position);
		Xoshiro xs_ypos(0.0, max_y_position);
		Xoshiro xs_xvel(-max_x_velocity, max_x_velocity);
		Xoshiro xs_yvel(-max_y_velocity, max_y_velocity);
		Xoshiro xs_angle(0.0, max_angle);

		for (int i = 0; i < MPCD::Constants::number; i++) {
			double mers_x = rg_xpos.next();
			double mers_y = rg_ypos.next();
			double mers_vx = rg_xvel.next();
			double mers_vy = rg_yvel.next();
			double mers_alpha = rg_angle.next();

			Vector2d pos(mers_x, mers_y);
			Vector2d vel(mers_vx, mers_vy);

			Particle mers_p(pos, vel);

			mers_particles.push_back(mers_p);
			mers_angles.push_back(mers_alpha);

			double xs_x = xs_xpos.next();
			double xs_y = xs_ypos.next();
			double xs_vx = xs_xvel.next();
			double xs_vy = xs_yvel.next();
			double xs_alpha = xs_angle.next();

			Vector2d pos_xs(xs_x, xs_y);
			Vector2d vel_xs(xs_vx, xs_vy);

			Particle xs_p(pos_xs, vel_xs);

			xs_particles.push_back(xs_p);
			xs_angles.push_back(xs_alpha);
		}
	}
};

/* Draw a histogram somehow, also check by hand. */
TEST_F(JupyterNotebookTests, MersWrite) { // DISABLED_

	std::filesystem::path cwd = std::filesystem::current_path();
	Out out(cwd.string() + l_data + l_rng);
	out.writeToOut(mers_angles, "mersenne_angles.csv", "alpha");
	out.writeToOut(mers_particles, "mersenne_particles.csv");

	for (int i = 0; i < MPCD::Constants::number; i++) { // initial + 1 move
		mers_particles[i].stream(MPCD::Constants::time_lapse);
	}
	out.writeToOut(mers_particles, "mersenne_particles_after_move.csv");

	int timesteps = 10;
	for (int t = 0; t < timesteps; t++) { // timesteps moves
		for (int i = 0; i < MPCD::Constants::number; i++) {
			mers_particles[i].stream(MPCD::Constants::time_lapse);
		}
	}
	out.writeToOut(mers_particles, "mersenne_particles_after_xx_timesteps.csv");
}

/* Draw a histogram somehow, also check by hand. */
TEST_F(JupyterNotebookTests, XSWrite) { // DISABLED_
	std::filesystem::path cwd = std::filesystem::current_path();
	Out out(cwd.string() + l_data + l_rng);
	out.writeToOut(xs_angles, "xoshiro_angles.csv", "alpha");
	out.writeToOut(xs_particles, "xoshiro_particles.csv");

	for (int i = 0; i < MPCD::Constants::number; i++) { // initial + 1 move
		xs_particles[i].stream(MPCD::Constants::time_lapse);
	}
	out.writeToOut(xs_particles, "xoshiro_particles_after_move.csv");

	int timesteps = 10;
	for (int t = 0; t < timesteps; t++) { // timesteps moves
		for (int i = 0; i < MPCD::Constants::number; i++) {
			xs_particles[i].stream(MPCD::Constants::time_lapse);
		}
	}
	out.writeToOut(xs_particles, "xoshiro_particles_after_xx_timesteps.csv");
}