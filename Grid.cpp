#include "Grid.h"
#include "Constants.h"

int MPCD::Grid::convertToLinearIndex(Eigen::Vector2i index) {
	return index[0] * MPCD::Constants::Grid::cols + index[1];
}
