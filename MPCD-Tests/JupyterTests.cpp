#include "pch.h"
#include "MersenneTwister.h"
#include "Xoshiro.h"
#include "Particle.h"
#include <filesystem>
#include <fstream>
#include "seeder.h"
#include "MPCD.h"

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
	const double aspect_ratio = Pipe::width / Pipe::height;
	const double max_x_position = Pipe::width;
	const double max_y_position = Pipe::height;
	const double max_x_velocity = std::max(max_x_position, max_y_position) / 100.0; // also in RandomGenerator.cpp
	const double max_y_velocity = max_x_velocity / aspect_ratio;
	const double max_angle = 2 * M_PI;

	std::vector<Particle> mers_particles;
	std::vector<double> mers_angles;

	std::vector<Particle> xs_particles;
	std::vector<double> xs_angles;

	void SetUp() override {
		mers_particles.reserve(number);
		mers_angles.reserve(number);

		xs_particles.reserve(number);
		xs_angles.reserve(number);

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

		for (int i = 0; i < number; i++) {
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

	std::ofstream outFileMersPos(cwd.string() + "//Data//RNG//mers_initial_xy.csv"); // navigate to MPCD and create /data folder
	outFileMersPos << "x,y" << "\n"; // header columns

	std::ofstream outFileMersVel(cwd.string() + "//Data//RNG//mers_initial_v_xy.csv");
	outFileMersVel << "vx,vy" << "\n"; // header columns

	std::ofstream outFileMersMovedPos(cwd.string() + "//Data//RNG//mers_moved_xy.csv");
	outFileMersMovedPos << "x,y" << "\n"; // header columns

	std::ofstream outFileMersTimestepsMovedPos(cwd.string() + "//Data//RNG//mers_timesteps_moved_xy.csv");
	outFileMersTimestepsMovedPos << "x,y" << "\n"; // header columns

	std::ofstream outFileMersAngle(cwd.string() + "//Data//RNG//mers_random_angle.csv");
	outFileMersAngle << "alpha" << "\n"; // header columns


	for (int i = 0; i < number; i++) { // initial + 1 move
		Vector2d pos_mers = mers_particles[i].getPosition();
		outFileMersPos << pos_mers(0) << "," << pos_mers(1) << "\n";
		Vector2d vel = mers_particles[i].getVelocity();
		outFileMersVel << vel(0) << "," << vel(1) << "\n";
		mers_particles[i].move(time_step);
		Vector2d moved_pos = mers_particles[i].getPosition();
		outFileMersMovedPos << moved_pos(0) << "," << moved_pos(1) << "\n";
		outFileMersAngle << mers_angles[i] << "\n";
	}

	int timesteps = 10;
	for (int t = 0; t < timesteps; t++) { // timesteps moves
		for (int i = 0; i < number; i++) {
			mers_particles[i].move(time_step);
		}
	}

	for (int i = 0; i < number; i++) {
		Vector2d pos = mers_particles[i].getPosition();
		outFileMersTimestepsMovedPos << pos(0) << "," << pos(1) << "\n";
	}

	outFileMersPos.close();
	outFileMersVel.close();
	outFileMersMovedPos.close();
	outFileMersTimestepsMovedPos.close();
	outFileMersAngle.close();
}

/* Draw a histogram somehow, also check by hand. */
TEST_F(JupyterNotebookTests, XSWrite) { // DISABLED_

	std::filesystem::path cwd = std::filesystem::current_path();

	std::ofstream outFileXSPos(cwd.string() + "//Data//RNG//xs_initial_xy.csv"); // navigate to MPCD and create /data folder
	outFileXSPos << "x,y" << "\n"; // header columns

	std::ofstream outFileXSVel(cwd.string() + "//Data//RNG//xs_initial_v_xy.csv");
	outFileXSVel << "vx,vy" << "\n"; // header columns

	std::ofstream outFileXSMovedPos(cwd.string() + "//Data//RNG//xs_moved_xy.csv");
	outFileXSMovedPos << "x,y" << "\n"; // header columns

	std::ofstream outFileXSTimestepsMovedPos(cwd.string() + "//Data//RNG//xs_timesteps_moved_xy.csv");
	outFileXSTimestepsMovedPos << "x,y" << "\n"; // header columns

	std::ofstream outFileXSAngle(cwd.string() + "//Data//RNG//xs_random_angle.csv");
	outFileXSAngle << "alpha" << "\n"; // header columns


	for (int i = 0; i < number; i++) { // initial + 1 move
		Vector2d pos = xs_particles[i].getPosition();
		outFileXSPos << pos(0) << "," << pos(1) << "\n";
		Vector2d vel = xs_particles[i].getVelocity();
		outFileXSVel << vel(0) << "," << vel(1) << "\n";
		xs_particles[i].move(time_step);
		Vector2d moved_pos = xs_particles[i].getPosition();
		outFileXSMovedPos << moved_pos(0) << "," << moved_pos(1) << "\n";
		outFileXSAngle << xs_angles[i] << "\n";
	}

	int timesteps = 10;
	for (int t = 0; t < timesteps; t++) { // timesteps moves
		for (int i = 0; i < number; i++) {
			xs_particles[i].move(time_step);
		}
	}

	for (int i = 0; i < number; i++) {
		Vector2d pos = xs_particles[i].getPosition();
		outFileXSTimestepsMovedPos << pos(0) << "," << pos(1) << "\n";
	}

	outFileXSPos.close();
	outFileXSVel.close();
	outFileXSMovedPos.close();
	outFileXSTimestepsMovedPos.close();
	outFileXSAngle.close();
}
