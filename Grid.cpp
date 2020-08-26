#include "Grid.h"

#define _USE_MATH_DEFINES
#include <math.h>

using namespace MPCD;
using namespace Eigen;

Grid::Grid(int minNumberPerCell) {
	min_particles_per_cell = minNumberPerCell;
	min_num_cells = MPCD::Constants::number / min_particles_per_cell;
	wanted_num_cells = min_num_cells * 2;
	cell_dim = std::sqrt(std::pow(MPCD::Constants::Pipe::width, 2) / (wanted_num_cells * MPCD::Constants::Pipe::aspect_ratio));
	max_shift = cell_dim / 2;
	rows = std::ceil(MPCD::Constants::Pipe::height / cell_dim);
	cols = std::ceil(MPCD::Constants::Pipe::width / cell_dim);
	rg_angle.setMin(0.0);
	rg_angle.setMax(2 * M_PI);
	rg_shift_x.setMin(-max_shift);
	rg_shift_y.setMax(max_shift);
	rg_shift_y.setMin(-max_shift);
	rg_shift_y.setMax(max_shift);
}

Grid::~Grid() {}

std::tuple<std::map<int, Vector2d>, std::map<int, double>> Grid::calculateCellValues(std::vector<Particle> particles) {
	std::map<int, Vector2d> meanCellVelocities;
	std::map<int, int> numParticles;
	std::map<int, double> rotationAngles;
	std::map<int, bool> calculationDone;

	for (auto it = particles.begin(); it != particles.end(); ++it) {
		Vector2d shift(rg_shift_x.next(), rg_shift_y.next());
		Vector2i cell_index = it->shift(shift, cell_dim);
		int linear_index = convertToLinearIndex(cell_index);
		Vector2d vel = it->getVelocity();
		meanCellVelocities[linear_index] += vel;
		numParticles[linear_index] += 1;
		if (!calculationDone[linear_index]) {
			rotationAngles[linear_index] = rg_angle.next();
		}
	}

	return std::make_pair(meanCellVelocities, rotationAngles);
}

int Grid::convertToLinearIndex(Eigen::Vector2i index) {
	return index[0] * cols + index[1];
}

double Grid::getCellDim() {
	return cell_dim;
}

double Grid::getMaxShift() {
	return max_shift;
}

int Grid::getMinParticlesPerCell() {
	return min_particles_per_cell;
}

int Grid::getMinNumCells() {
	return min_num_cells;
}

int Grid::getWantedNumCells() {
	return wanted_num_cells;
}

int Grid::getRows() {
	return rows;
}

int Grid::getCols() {
	return cols;
}