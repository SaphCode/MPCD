#include "pch.h"
#include "MersenneTwister.h"
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <Eigen/Dense>
#include "Constants.h"

/* This and only this order gives us pi */
#define _USE_MATH_DEFINES
#include <math.h>

using namespace Eigen;
using namespace MPCD;

class MersenneTwisterTest : public ::testing::Test {
protected:
	const double time_step = 1.0;
	const double aspect_ratio = MPCD::Constants::Pipe::width / MPCD::Constants::Pipe::height;
	const double min_x_position = MPCD::Constants::Pipe::x_0;
	const double min_y_position = MPCD::Constants::Pipe::y_0;
	const double max_x_position = MPCD::Constants::Pipe::x_max;
	const double max_y_position = MPCD::Constants::Pipe::y_max;
	const double max_x_velocity = std::max(MPCD::Constants::Pipe::width, MPCD::Constants::Pipe::height) / 100.0; // also in RandomGenerator.cpp
	const double max_y_velocity = max_x_velocity / aspect_ratio;
	const double max_angle = 2 * M_PI;

	std::vector<Eigen::Vector2d> mers_positions;
	std::vector<Eigen::Vector2d> mers_velocities;
	std::vector<double> mers_angles;

	void SetUp() override {
		mers_positions.reserve(MPCD::Constants::number);
		mers_velocities.reserve(MPCD::Constants::number);
		mers_angles.reserve(MPCD::Constants::number);

		MersenneTwister rg_xpos(min_x_position, max_x_position, DistributionType::UNIFORM);
		MersenneTwister rg_ypos(min_y_position, max_y_position, DistributionType::UNIFORM);
		MersenneTwister rg_xvel(-max_x_velocity, max_x_velocity, DistributionType::UNIFORM);
		MersenneTwister rg_yvel(-max_y_velocity, max_y_velocity, DistributionType::UNIFORM);
		MersenneTwister rg_angle(0.0, max_angle, DistributionType::UNIFORM);

		for (int i = 0; i < MPCD::Constants::number; i++) {
			double mers_x = rg_xpos.next();
			double mers_y = rg_ypos.next();
			double mers_vx = rg_xvel.next();
			double mers_vy = rg_yvel.next();
			double mers_alpha = rg_angle.next();

			Vector2d pos(mers_x, mers_y);
			Vector2d vel(mers_vx, mers_vy);

			mers_positions.push_back(pos);
			mers_velocities.push_back(vel);
			mers_angles.push_back(mers_alpha);
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

TEST_F(MersenneTwisterTest, RandomBoundsTest) {

	for (int i = 0; i < MPCD::Constants::number; i++) {
		Vector2d mers_pos(mers_positions[i]);
		Vector2d mers_vel(mers_velocities[i]);
		double mers_alpha = mers_angles[i];

		ASSERT_GE(mers_pos(0), min_x_position);
		ASSERT_LE(mers_pos(0), max_x_position);
		ASSERT_GE(mers_pos(1), min_y_position);
		ASSERT_LE(mers_pos(1), max_y_position);
		ASSERT_GE(mers_vel(0), -max_x_velocity);
		ASSERT_LE(mers_vel(0), max_x_velocity);
		ASSERT_GE(mers_vel(1), -max_y_velocity);
		ASSERT_LE(mers_vel(1), max_y_velocity);
		ASSERT_GE(mers_alpha, 0.0);
		ASSERT_LE(mers_alpha, 2 * M_PI);

	}
}

TEST_F(MersenneTwisterTest, ChiSquaredTest) {

	/* Chi Squared Testing */
	std::vector<double> chi_2_alpha05; // at least 2 degrees of freedom, 1 undefined for alpha = 0.99
	chi_2_alpha05.push_back(3.841);
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
	chi_2_alpha05.push_back(30.144); // 19

	std::filesystem::path cwd = std::filesystem::current_path();
	std::ofstream outFile;
	outFile.open(cwd.string() + "//Data//RNG//mersenne_twister_chi2.csv"); // navigate to MPCD and create /data folder
	outFile << "k,observed MT,chi^2 probability" << std::endl;
	outFile.close();
	for (int buckets = 2; buckets <= chi_2_alpha05.size(); buckets++) {
		const double expected_bucket_size = (double)MPCD::Constants::number / buckets;

		std::vector<std::vector<double>> mers_buckets_xpos;
		std::vector<std::vector<double>> mers_buckets_ypos;
		std::vector<std::vector<double>> mers_buckets_xvel;
		std::vector<std::vector<double>> mers_buckets_yvel;
		std::vector<std::vector<double>> mers_buckets_angle;
		mers_buckets_xpos.reserve(buckets);
		mers_buckets_ypos.reserve(buckets);
		mers_buckets_xvel.reserve(buckets);
		mers_buckets_yvel.reserve(buckets);
		mers_buckets_angle.reserve(buckets);

		double mers_diffb_xpos_expected_2 = 0;
		double mers_diffb_ypos_expected_2 = 0;
		double mers_diffb_xvel_expected_2 = 0;
		double mers_diffb_yvel_expected_2 = 0;
		double mers_diffb_angle_expected_2 = 0;



		for (int b = 0; b < buckets; b++) {
			std::vector<double> mers_b_xpos;
			std::vector<double> mers_b_ypos;
			std::vector<double> mers_b_xvel;
			std::vector<double> mers_b_yvel;
			std::vector<double> mers_b_angle;

			mers_b_xpos.reserve(std::round(expected_bucket_size));
			mers_b_ypos.reserve(std::round(expected_bucket_size));
			mers_b_xvel.reserve(std::round(expected_bucket_size));
			mers_b_yvel.reserve(std::round(expected_bucket_size));
			mers_b_angle.reserve(std::round(expected_bucket_size));

			ASSERT_EQ(mers_positions.size(), MPCD::Constants::number);
			ASSERT_EQ(mers_velocities.size(), MPCD::Constants::number);
			ASSERT_EQ(mers_angles.size(), MPCD::Constants::number);

			auto size = mers_positions.size();
			for (int i = 0; i < size; i++) {

				Vector2d pos = mers_positions[i];
				if (inBucketB(buckets, b, pos(0), min_x_position, max_x_position)) {
					mers_b_xpos.push_back(pos(0));
				}
				if (inBucketB(buckets, b, pos(1), min_y_position, max_y_position)) {
					mers_b_ypos.push_back(pos(1));
				}

				Vector2d vel = mers_velocities[i];
				if (inBucketB(buckets, b, vel(0), -max_x_velocity, max_x_velocity)) {
					mers_b_xvel.push_back(vel(0));
				}
				if (inBucketB(buckets, b, vel(1), -max_y_velocity, max_y_velocity)) {
					mers_b_yvel.push_back(vel(1));
				}

				double alpha = mers_angles[i];
				if (inBucketB(buckets, b, alpha, 0.0, max_angle)) {
					mers_b_angle.push_back(alpha);
				}
			}
			mers_buckets_xpos.push_back(mers_b_xpos);
			mers_buckets_ypos.push_back(mers_b_ypos);
			mers_buckets_xvel.push_back(mers_b_xvel);
			mers_buckets_yvel.push_back(mers_b_yvel);
			mers_buckets_angle.push_back(mers_b_angle);

			mers_diffb_xpos_expected_2 += std::pow(mers_buckets_xpos[b].size() - expected_bucket_size, 2);
			mers_diffb_ypos_expected_2 += std::pow(mers_buckets_yvel[b].size() - expected_bucket_size, 2);
			mers_diffb_xvel_expected_2 += std::pow(mers_buckets_xvel[b].size() - expected_bucket_size, 2);
			mers_diffb_yvel_expected_2 += std::pow(mers_buckets_yvel[b].size() - expected_bucket_size, 2);
			mers_diffb_angle_expected_2 += std::pow(mers_buckets_angle[b].size() - expected_bucket_size, 2);


		}

		double mers_obs_chi_2_xpos = (double)mers_diffb_xpos_expected_2 / expected_bucket_size;
		double mers_obs_chi_2_ypos = (double)mers_diffb_ypos_expected_2 / expected_bucket_size;
		double mers_obs_chi_2_xvel = (double)mers_diffb_xvel_expected_2 / expected_bucket_size;
		double mers_obs_chi_2_yvel = (double)mers_diffb_yvel_expected_2 / expected_bucket_size;
		double mers_obs_chi_2_angle = (double)mers_diffb_angle_expected_2 / expected_bucket_size;



		outFile.open(cwd.string() + "//Data//RNG//mersenne_twister_chi2.csv", std::ios_base::app); // navigate to MPCD and create /data folder
		outFile << std::fixed << std::setprecision(3);
		outFile << buckets - 1 << "," << mers_obs_chi_2_angle << "," << chi_2_alpha05[buckets - 1] << std::endl;
		outFile.close();

		EXPECT_LE(mers_obs_chi_2_xpos, chi_2_alpha05[buckets - 1]);
		EXPECT_LE(mers_obs_chi_2_ypos, chi_2_alpha05[buckets - 1]);
		EXPECT_LE(mers_obs_chi_2_xvel, chi_2_alpha05[buckets - 1]);
		EXPECT_LE(mers_obs_chi_2_yvel, chi_2_alpha05[buckets - 1]);
		EXPECT_LE(mers_obs_chi_2_angle, chi_2_alpha05[buckets - 1]);
	}

}

