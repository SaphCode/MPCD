#include "pch.h"
#include "MPCD.h"
#include "Grid.h"

using namespace MPCD;

TEST(MPCD, Timestep) {
	Grid g(MPCD::Constants::min_particles_per_cell);
	std::vector<Particle> particles;
	particles.reserve(MPCD::Constants::number);

	/* dont worry the numbers are just seeds */
	Xoshiro xs_xpos(0.0, 1.0);
	Xoshiro xs_ypos(0.0, 1.0);
	Xoshiro xs_xvel(-0.01, 0.01);
	Xoshiro xs_yvel(-0.01, 0.01);
	Xoshiro angles(0.0, 2 * M_PI);
	Xoshiro rg_angle = angles;

	for (int i = 0; i < MPCD::Constants::number; i++) {
		double xs_x = xs_xpos.next();
		double xs_y = xs_ypos.next();
		double xs_vx = xs_xvel.next();
		double xs_vy = xs_yvel.next();

		Vector2d pos(xs_x, xs_y);
		Vector2d vel(xs_vx, xs_vy);

		Particle xs_p(pos, vel);

		particles.push_back(xs_p);
	}

	double time_lapse = MPCD::Constants::time_lapse;

	Xoshiro rg_angle(0.0, 2 * M_PI);
	double max_shift = g.getMaxShift();
	Xoshiro rg_shift_x(-max_shift, max_shift);
	Xoshiro rg_shift_y(-max_shift, max_shift);

	/*

	Grid g(min_particles_per_cell);

	std::tuple<std::map<int, Vector2d>, std::map<int, double>> origMaps = g.calculateCellValues(particles);
	std::map<int, Vector2d> origMeanVelocities = std::get<0>(origMaps);
	std::map<int, double> origRotationAngles = std::get<1>(origMaps);

	ASSERT_EQ(origMeanVelocities.size(), origRotationAngles.size());

	int steps = 1;
	//for (int step = 0; step < steps; step++) {
	for (auto it = particles.begin(); it != particles.end(); ++it) {
		it->move(time_step);
	}

	std::tuple<std::map<int, Vector2d>, std::map<int, double>> maps = g.calculateCellValues(particles);
	std::map<int, Vector2d> meanCellVelocities = std::get<0>(maps);
	std::map<int, double> rotationAngles = std::get<1>(maps);

	double cell_dim = g.getCellDim();
	for (auto it = particles.begin(); it != particles.end(); ++it) {
		Vector2i cellIndex = it->getCellIndex(it->getPosition(), cell_dim);
		int linearIndex = g.convertToLinearIndex(cellIndex);
		it->updateVelocity(meanCellVelocities[linearIndex], rotationAngles[linearIndex]);
	}

	ASSERT_EQ(true, false); // this is here to remind yourself to fix this test. in mpcd.cpp the new timestep method was built, test that instead.
	//}	

	std::tuple<std::map<int, Vector2d>, std::map<int, double>> maps_new = g.calculateCellValues(particles);
	std::map<int, Vector2d> meanCellVelocities_new = std::get<0>(maps_new);
	std::map<int, double> rotationAngles_new = std::get<1>(maps_new);

	ASSERT_EQ(meanCellVelocities_new.size(), rotationAngles_new.size());

	std::filesystem::path cwd = std::filesystem::current_path();

	Out out(cwd.string() + l_data);
	out.writeToOut(origMeanVelocities, "before_collision.csv", "cell_vx_b,cell_vy_b");
	out.writeToOut(meanCellVelocities_new, "after_collision.csv", "cell_vx_a,cell_vy_a");
	*/
}
