#include "Grid.h"
#include "Constants.h"

MPCD::Grid::Grid() {
	_a = MPCD::Constants::Grid::cell_dim;
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
	div_t divresult_i = div(shiftedPos[1], 
	int i = std::floor(position[1] / _a);
	assert(i >= -1 && i <= MPCD::Constants::Grid::num_cols + 1);
	int j = std::floor(position[0] / _a);
	assert(j >= -1 && j <= MPCD::Constants::Grid::num_rows + 1);
	return std::make_pair(i, j);
}

void MPCD::Grid::insert(Particle p) {
	std::pair<int, int> key = getCoordinates(p.getPosition());
	if (_cells.count(key) == 0) {
		Cell c;
		_cells.insert(std::make_pair(key, c));
	}
	_cells[key].add(p);
}


void MPCD::Grid::shift() {
	Eigen::Vector2d shift(_shiftGen.next(), _shiftGen.next());
	_shift = shift;
}

void MPCD::Grid::undoShift() {
	Eigen::Vector2d zero(0, 0);
	_shift = zero;
}


MPCD::Grid MPCD::operator+(const Grid& lhs, const Grid& rhs)
{
	Grid grid;
	grid._cells = lhs._cells;
	for (const auto& [key, val] : rhs._cells) {
		grid._cells[key] = grid._cells[key] + val;
	}
	return grid;
}
