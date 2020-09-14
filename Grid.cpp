#include "Grid.h"
#include "Constants.h"
#include <execution>

MPCD::Grid::Grid() {
	for (int i = 0; i < _num_rows; i++) {
		for (int j = 0; j < _num_cols; j++) {
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
		assert(coordinates.first >= 0 && coordinates.first <= _num_rows);
		assert(coordinates.second >= 0 && coordinates.second <= _num_cols);
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
std::pair<int, int> MPCD::Grid::getCoordinates(Eigen::Vector2d position) {
	Eigen::Vector2d shiftedPos = position + _shift;
	int i = std::floor(shiftedPos[1] / _a);
	assert(i >= -1 && i <= _num_rows + 1);
	if (i == -1) {
		i = _num_rows;
	}
	else if (i == _num_rows + 1) {
		i = 0;
	}
	int j = std::floor(shiftedPos[0] / _a);
	assert(j >= -1 && j <= _num_cols + 1);
	if (j == -1) {
		j = _num_cols;
	}
	else if (j == _num_cols + 1) {
		j = 0;
	}
	return std::make_pair(i, j);
}

void MPCD::Grid::collision(bool draw, std::ofstream& outFile)
{
	std::mutex m;
	for (auto& [key, cell] : _cells) {
		cell.collide();
		if (draw) {
			cell.draw(m, key, outFile);
		}
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

int MPCD::Grid::getAverageParticlesPerCell()
{
	return _average_particles_per_cell;
}

double MPCD::Grid::getA()
{
	return _a;
}

int MPCD::Grid::getNumRows()
{
	return _num_rows;
}

int MPCD::Grid::getNumCols()
{
	return _num_cols;
}

double MPCD::Grid::getMaxShift()
{
	return _max_shift;
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