#include "Grid.h"
#include "Constants.h"
#include <execution>
#include "MaxwellBoltzmann.h"

using namespace Eigen;

MPCD::Grid::Grid() :
	_shiftGen(-MPCD::Constants::cell_dim / 2, MPCD::Constants::cell_dim / 2),
	_numRows(std::floor((MPCD::Constants::y_max - MPCD::Constants::y_0) / MPCD::Constants::cell_dim)),
	_numCols(std::floor((MPCD::Constants::x_max - MPCD::Constants::x_0) / MPCD::Constants::cell_dim)),
	_a(MPCD::Constants::cell_dim),
	_maxShift(MPCD::Constants::cell_dim / 2),
	_average_particles_per_cell(MPCD::Constants::average_particles_per_cell)
{
	for (int i = 0; i < _numRows; i++) {
		for (int j = 0; j < _numCols; j++) {
			Cell cell;
			std::pair<int, int> coords(i, j);
			_cells.insert(std::make_pair(coords, cell));
		}
	}
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
		_cells[coordinates].add(p);
	}
}

/*
void MPCD::Grid::updateCell(Particle p, Eigen::Vector2d positionBeforeMove) {
	Eigen::Vector2d zero(0, 0);
	if (_shift == zero) {
		shift();
	}

	Eigen::Vector2d position = p.getPosition();
	const std::pair<int, int> coordsBeforeMove = getCoordinates(positionBeforeMove);
	const std::pair<int, int> coords = getCoordinates(position);
	const std::pair<int, int> coordsAfterShift = getCoordinates(position + _shift);
	if (coords != coordsBeforeMove) {
		_cells[coordsBeforeMove].remove(p);

		bool cellExists = _cells.count(coords) > 0;
		if (!cellExists) {
			MPCD::Cell cell;
			_cells.insert(std::make_pair(coords, cell));
		}
		_cells[coords].add(p);
	}
	if (coords != coordsAfterShift) {
		_cellsShifted[coords].remove(p); // this is sure to exits, since it is a copy of _cells

		bool cellExists = _cellsShifted.count(coordsAfterShift) > 0;
		if (!cellExists) {
			MPCD::Cell cell;
			_cellsShifted.insert(std::make_pair(coordsAfterShift, cell));
		}
		_cellsShifted[coordsAfterShift].add(p);
	}
	
}
*/
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
	std::mutex m;
	for (auto& [key, cell] : _cells) {
		const double cell_dim = MPCD::Constants::cell_dim;
		const int lastRow = int(std::round(MPCD::Constants::y_max / cell_dim)) - 1;
		const int firstRow = 0;
		if ((key.first == firstRow) || (key.first == lastRow)) {
			createVirtualParticles(key, cell, firstRow, lastRow, cell_dim);
		}
		cell.collide();
		if (draw) {
			cell.draw(m, key, outFile);
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
		cell.add(p);
	}
}

/*
void MPCD::Grid::insert(Particle p) {
	std::pair<int, int> key = getCoordinates(p.getPosition());
	if (_cells.count(key) == 0) {
		Cell c;
		_cells.insert(std::make_pair(key, c));
	}
	_cells[key].add(p);
}
*/

void MPCD::Grid::shift() {
	Eigen::Vector2d shift(_shiftGen.next(), _shiftGen.next());
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


/*
MPCD::Grid MPCD::operator+(const Grid& lhs, const Grid& rhs)
{
	Grid grid;
	grid._cells = lhs._cells;
	for (const auto& [key, val] : rhs._cells) {
		grid._cells[key] = grid._cells[key] + val;
	}
	return grid;
}
*/