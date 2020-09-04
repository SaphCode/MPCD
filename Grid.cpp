#include "Grid.h"
#include "Constants.h"

int MPCD::Grid::convertToLinearIndex(Eigen::Vector2i index) {
	return index[0] * MPCD::Constants::Grid::max_cols + index[1];
}

Eigen::Vector2i MPCD::Grid::convertToIndex(int linearIndex) {
	Eigen::Vector2i index(0, 0);
	int max_cols = MPCD::Constants::Grid::max_cols;
	index[0] = std::floor(linearIndex / max_cols);
	index[1] = linearIndex % max_cols;
	return index;
}

int MPCD::Grid::convertToLinearIndex(Eigen::Vector2i index, int cols) {
	return index[0] * cols + index[1];
}

Eigen::Vector2i MPCD::Grid::convertToIndex(int linearIndex, int cols) {
	Eigen::Vector2i index(0, 0);
	index[0] = std::floor(linearIndex / cols);
	index[1] = linearIndex % cols;
	return index;
}
