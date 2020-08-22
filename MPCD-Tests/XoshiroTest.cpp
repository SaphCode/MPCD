#include "pch.h"
#include "Xoshiro.h"
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <Eigen/Dense>
#include "MPCD.h"
#include "seeder.h"

/* This and only this order gives us pi */
#define _USE_MATH_DEFINES
#include <math.h>

using namespace Eigen;
using namespace MPCD;

class XoshiroTest : public ::testing::Test {
protected:
	const double time_step = 1.0;
	const double aspect_ratio = Pipe::width / Pipe::height;
	const double max_x_position = Pipe::width;
	const double max_y_position = Pipe::height;
	const double max_x_velocity = std::max(Pipe::width, Pipe::height) / 100.0;
	const double max_y_velocity = max_x_velocity / aspect_ratio;
	const double max_angle = 2 * M_PI;

	std::vector<Eigen::Vector2d> xs_positions;
	std::vector<Eigen::Vector2d> xs_velocities;
	std::vector<double> xs_angles;

	void SetUp() override {
		xs_positions.reserve(number);
		xs_velocities.reserve(number);
		xs_angles.reserve(number);

		/* dont worry the numbers are just seeds */
		Xoshiro xs_xpos(0.0, max_x_position);
		Xoshiro xs_ypos(0.0, max_y_position);
		Xoshiro xs_xvel(-max_x_velocity, max_x_velocity);
		Xoshiro xs_yvel(-max_y_velocity, max_y_velocity);
		Xoshiro xs_angle(0.0, max_angle);

		for (int i = 0; i < number; i++) {
			double xs_x = xs_xpos.next();
			double xs_y = xs_ypos.next();
			double xs_vx = xs_xvel.next();
			double xs_vy = xs_yvel.next();
			double xs_alpha = xs_angle.next();

			Vector2d pos(xs_x, xs_y);
			Vector2d vel(xs_vx, xs_vy);

			xs_positions.push_back(pos);
			xs_velocities.push_back(vel);
			xs_angles.push_back(xs_alpha);
		}
	}

	bool inBucketB(int buckets, int b, double variable, double min, double max) {
		double range = max - min;
		double interval = range / buckets;
		double bucket_min = min + b * interval;
		double bucket_max = min + (b + 1) * interval;
		bool result = (variable >= bucket_min) && (variable < bucket_max);
		return result;
	}
};


TEST_F(XoshiroTest, RandomBoundsTest) {

	for (int i = 0; i < number; i++) {
		Vector2d xs_pos(xs_positions[i]);
		Vector2d xs_vel(xs_velocities[i]);
		double xs_alpha = xs_angles[i];

		ASSERT_GE(xs_pos(0), 0.0);
		ASSERT_LE(xs_pos(0), Pipe::width);
		ASSERT_GE(xs_pos(1), 0.0);
		ASSERT_LE(xs_pos(1), Pipe::height);
		ASSERT_GE(xs_vel(0), -max_x_velocity);
		ASSERT_LE(xs_vel(0), max_x_velocity);
		ASSERT_GE(xs_vel(1), -max_y_velocity);
		ASSERT_LE(xs_vel(1), max_y_velocity);
		ASSERT_GE(xs_alpha, 0.0);
		ASSERT_LE(xs_alpha, 2 * M_PI);
	}
}

TEST_F(XoshiroTest, ChiSquaredTest) {

	/* Chi Squared Testing */
	std::vector<double> chi_2_alpha05; // at least 2 degrees of freedom, 1 undefined for alpha = 0.99
	chi_2_alpha05.push_back(3.841); // undefined for 1 and alpha .01
	chi_2_alpha05.push_back(5.991); // 2
	chi_2_alpha05.push_back(7.815);
	chi_2_alpha05.push_back(9.488);
	chi_2_alpha05.push_back(11.070);
	chi_2_alpha05.push_back(12.592);
	chi_2_alpha05.push_back(14.067);
	chi_2_alpha05.push_back(15.507);
	chi_2_alpha05.push_back(16.919);
	chi_2_alpha05.push_back(18.307); // 10
	chi_2_alpha05.push_back(19.675);
	chi_2_alpha05.push_back(21.026);
	chi_2_alpha05.push_back(22.362);
	chi_2_alpha05.push_back(23.685);
	chi_2_alpha05.push_back(24.996);
	chi_2_alpha05.push_back(26.296);
	chi_2_alpha05.push_back(27.587);
	chi_2_alpha05.push_back(28.869);
	chi_2_alpha05.push_back(30.144);

	std::filesystem::path cwd = std::filesystem::current_path();
	std::ofstream outFile;
	outFile.open(cwd.string() + "//Data//RNG//xoshiro_chi2.csv"); // navigate to MPCD and create /data folder
	outFile << "k,observed XS256++,chi^2 probability" << std::endl;
	outFile.close();

	for (int buckets = 2; buckets <= chi_2_alpha05.size(); buckets++) {
		const double expected_bucket_size = (double)number / buckets;

		std::vector<std::vector<double>> xs_buckets_xpos;
		std::vector<std::vector<double>> xs_buckets_ypos;
		std::vector<std::vector<double>> xs_buckets_xvel;
		std::vector<std::vector<double>> xs_buckets_yvel;
		std::vector<std::vector<double>> xs_buckets_angle;
		xs_buckets_xpos.reserve(buckets);
		xs_buckets_ypos.reserve(buckets);
		xs_buckets_xvel.reserve(buckets);
		xs_buckets_yvel.reserve(buckets);
		xs_buckets_angle.reserve(buckets);

		double xs_diffb_xpos_expected_2 = 0;
		double xs_diffb_ypos_expected_2 = 0;
		double xs_diffb_xvel_expected_2 = 0;
		double xs_diffb_yvel_expected_2 = 0;
		double xs_diffb_angle_expected_2 = 0;

		for (int b = 0; b < buckets; b++) {

			std::vector<double> xs_b_xpos;
			std::vector<double> xs_b_ypos;
			std::vector<double> xs_b_xvel;
			std::vector<double> xs_b_yvel;
			std::vector<double> xs_b_angle;

			xs_b_xpos.reserve(std::round(expected_bucket_size));
			xs_b_ypos.reserve(std::round(expected_bucket_size));
			xs_b_xvel.reserve(std::round(expected_bucket_size));
			xs_b_yvel.reserve(std::round(expected_bucket_size));
			xs_b_angle.reserve(std::round(expected_bucket_size));

			ASSERT_EQ(xs_positions.size(), number);
			ASSERT_EQ(xs_velocities.size(), number);
			ASSERT_EQ(xs_angles.size(), number);

			auto size = xs_positions.size();
			for (int i = 0; i < size; i++) {
				Vector2d pos = xs_positions[i];
				if (inBucketB(buckets, b, pos(0), Pipe::x_0, Pipe::width)) {
					xs_b_xpos.push_back(pos(0));
				}
				if (inBucketB(buckets, b, pos(1), Pipe::y_0, Pipe::height)) {
					xs_b_ypos.push_back(pos(1));
				}

				Vector2d vel = xs_velocities[i];
				if (inBucketB(buckets, b, vel(0), -max_x_velocity, max_x_velocity)) {
					xs_b_xvel.push_back(vel(0));
				}
				if (inBucketB(buckets, b, vel(1), -max_y_velocity, max_y_velocity)) {
					xs_b_yvel.push_back(vel(1));
				}

				double alpha = xs_angles[i];
				if (inBucketB(buckets, b, alpha, 0.0, max_angle)) {
					xs_b_angle.push_back(alpha);
				}
			}

			xs_buckets_xpos.push_back(xs_b_xpos);
			xs_buckets_ypos.push_back(xs_b_ypos);
			xs_buckets_xvel.push_back(xs_b_xvel);
			xs_buckets_yvel.push_back(xs_b_yvel);
			xs_buckets_angle.push_back(xs_b_angle);

			xs_diffb_xpos_expected_2 += std::pow(xs_buckets_xpos[b].size() - expected_bucket_size, 2);
			xs_diffb_ypos_expected_2 += std::pow(xs_buckets_ypos[b].size() - expected_bucket_size, 2);
			xs_diffb_xvel_expected_2 += std::pow(xs_buckets_xvel[b].size() - expected_bucket_size, 2);
			xs_diffb_yvel_expected_2 += std::pow(xs_buckets_yvel[b].size() - expected_bucket_size, 2);
			xs_diffb_angle_expected_2 += std::pow(xs_buckets_angle[b].size() - expected_bucket_size, 2);
		}

		double xs_obs_chi_2_xpos = (double)xs_diffb_xpos_expected_2 / expected_bucket_size;
		double xs_obs_chi_2_ypos = (double)xs_diffb_ypos_expected_2 / expected_bucket_size;
		double xs_obs_chi_2_xvel = (double)xs_diffb_xvel_expected_2 / expected_bucket_size;
		double xs_obs_chi_2_yvel = (double)xs_diffb_yvel_expected_2 / expected_bucket_size;
		double xs_obs_chi_2_angle = (double)xs_diffb_angle_expected_2 / expected_bucket_size;

		outFile.open(cwd.string() + "//Data//RNG//xoshiro_chi2.csv", std::ios_base::app); // navigate to MPCD and create /data folder
		outFile << std::fixed << std::setprecision(3);
		outFile << buckets - 1 << "," << xs_obs_chi_2_angle << "," << chi_2_alpha05[buckets - 1] << std::endl;
		outFile.close();

		EXPECT_LE(xs_obs_chi_2_xpos, chi_2_alpha05[buckets - 1]);
		EXPECT_LE(xs_obs_chi_2_ypos, chi_2_alpha05[buckets - 1]);
		EXPECT_LE(xs_obs_chi_2_xvel, chi_2_alpha05[buckets - 1]);
		EXPECT_LE(xs_obs_chi_2_yvel, chi_2_alpha05[buckets - 1]);
		EXPECT_LE(xs_obs_chi_2_angle, chi_2_alpha05[buckets - 1]);
	}
}
