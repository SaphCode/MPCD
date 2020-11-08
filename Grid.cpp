#include "Grid.h"
#include "Constants.h"
#include <execution>
#include "MaxwellBoltzmann.h"
#include "Cell.h"

using namespace Eigen;

MPCD::Grid::Grid() :
	_shiftGen{ std::random_device()() },
	_unif(-MPCD::Constants::cell_dim / 2, MPCD::Constants::cell_dim / 2),
	_numRows(std::floor((MPCD::Constants::y_max - MPCD::Constants::y_0) / MPCD::Constants::cell_dim)),
	_numCols(std::floor((MPCD::Constants::x_max - MPCD::Constants::x_0) / MPCD::Constants::cell_dim)),
	_a(MPCD::Constants::cell_dim),
	_maxShift(MPCD::Constants::cell_dim / 2),
	_average_particles_per_cell(MPCD::Constants::average_particles_per_cell)
{
}

void MPCD::Grid::updateCoordinates(std::vector<Particle>& particles)
{
	for (auto& [key, cell] : _cells) {
		cell.clear();
	}
	for (auto& p : particles) {
		Eigen::Vector2d particlePos = p.getPosition();
		std::pair<int, int> coordinates = getCoordinates(particlePos);
		assert(coordinates.first >= 0 && coordinates.first <= _numRows);
		assert(coordinates.second >= 0 && coordinates.second <= _numCols);
		_cells.at(coordinates).add(p);
	}
}

std::pair<int, int> MPCD::Grid::getCoordinates(Eigen::Vector2d position) const {
	Eigen::Vector2d shiftedPos = position - _shift;
	int i = std::floor(shiftedPos[1] / _a);
	assert(i >= -1 && i <= _numRows);
	if (i == -1) {
		i = _numRows - 1;
	}
	else if (i == _numRows) {
		i = 0;
	}
	int j = std::floor(shiftedPos[0] / _a);
	assert(j >= -1 && j <= _numCols);
	if (j == -1) {
		j = _numCols - 1;
	}
	else if (j == _numCols) {
		j = 0;
	}
	return std::make_pair(i, j);
}

void MPCD::Grid::collision(bool draw, std::ofstream& outFile)
{
	const double cell_dim = MPCD::Constants::cell_dim;
	const int lastRow = int(std::round(MPCD::Constants::y_max / cell_dim)) - 1;
	const int firstRow = 0;
	#pragma omp parallel for
	for (int i = 0; i <= lastRow; i++) {
		for (int j = 0; j < MPCD::Constants::x_max / cell_dim; j++) {
			std::pair<int, int> key(i, j);
			Cell& cell = _cells.at(key);
			if ((key.first == firstRow) || (key.first == lastRow)) {
				createVirtualParticles(key, cell, firstRow, lastRow, cell_dim);
			}
			// CHECK HERE FOR OBSTACLE CELLS
			double scaling = cell.thermostatScaling();
			#pragma omp critical
			{
				cell.collide(scaling);
				if (draw) {
					cell.draw(key, outFile);
				}
			}
		}
	}
		
}

void MPCD::Grid::createVirtualParticles(const std::pair<int, int>& key, Cell& cell, const int firstRow, const int lastRow, const double cell_dim) {
	const int numParticles = cell.number();
	const int particlesToCreate = MPCD::Constants::average_particles_per_cell - numParticles;
	const double mass = MPCD::Constants::particle_mass; // h2o kg mass
	const double mean = 0;
	const double temperature = MPCD::Constants::temperature;
	MaxwellBoltzmann mb_vel(mean, temperature, mass);
	for (int i = 0; i < particlesToCreate; i++) {
		Vector2d mockPos(key.second * cell_dim, key.first * cell_dim);
		Vector2d vel = mb_vel.next();
		Particle p(mass, mockPos, vel);
		cell.addVirtual(p);
	}
}

void MPCD::Grid::shift() {
	Eigen::Vector2d shift(_unif(_shiftGen), _unif(_shiftGen));
	_shift = shift;
}

void MPCD::Grid::undoShift() {
	Eigen::Vector2d zero(0, 0);
	_shift = zero;
}

int MPCD::Grid::getAverageParticlesPerCell() const
{
	return _average_particles_per_cell;
}

double MPCD::Grid::getA() const
{
	return _a;
}

int MPCD::Grid::getNumRows() const
{
	return _numRows;
}

int MPCD::Grid::getNumCols() const
{
	return _numCols;
}

double MPCD::Grid::getMaxShift() const
{
	return _maxShift;
}
