#include "Grid.h"
#include "Constants.h"

int MPCD::Grid::convertToLinearIndex(Eigen::Vector2i index, int cols) {
	return index[0] * cols + index[1];
}
