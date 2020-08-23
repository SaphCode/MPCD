#include "Grid.h"

#define _USE_MATH_DEFINES
#include <math.h>

using namespace MPCD;
using namespace Eigen;

Grid::Grid(int minNumberPerCell) {
	min_num_cells = number / minNumberPerCell;
	wanted_num_cells = min_num_cells * 2;
	cell_dim = std::sqrt(std::pow(MPCD::Pipe::width, 2) / (wanted_num_cells * MPCD::Pipe::aspect_ratio));
	max_shift = cell_dim / 2;
	rows = std::ceil(MPCD::Pipe::height / cell_dim);
	cols = std::ceil(MPCD::Pipe::width / cell_dim);
	rg_angle.setMin(0.0);
	rg_angle.setMax(2 * M_PI);
	rg_shift_x.setMin(-max_shift);
	rg_shift_y.setMax(max_shift);
	rg_shift_y.setMin(-max_shift);
	rg_shift_y.setMax(max_shift);
}

Grid::~Grid() {}

std::tuple<std::map<Vector2i, Vector2d>, std::map<Vector2i, double>> Grid::calculateCellValues(std::vector<Particle> particles) {
	std::map<Vector2i, Vector2d> meanCellVelocities;
	std::map<Vector2i, int> numParticles;
	std::map<Vector2i, double> rotationAngles;
	std::map<Vector2i, bool> calculationDone;

	for (auto it = particles.begin(); it != particles.end(); ++it) {
		Vector2d shift(rg_shift_x.next(), rg_shift_y.next());
		Vector2i cell_index = it->shift(shift);
		Vector2d vel = it->getVelocity();
		meanCellVelocities[cell_index] += vel;
		numParticles[cell_index] += 1;
		if (!calculationDone[cell_index]) {
			rotationAngles[cell_index] = rg_angle.next();
		}
	}

	return std::make_pair(meanCellVelocities, rotationAngles);
}

//static int convertToLinearIndex(Eigen::Vector2i index);