#pragma once
#include <Eigen/Dense>

#ifndef GRID_H
#define GRID_H

namespace MPCD {
	class Grid
	{
	public:
		Eigen::Vector2i getCell();
		/* Converts 2d indexes into linear indexes.
			Example:
			A = [[1, 2, 3]
				 [4, 5, 6]
				 [7, 8, 9]]
			A[0,0], index = (0,0)
			returns 0
			A[1,1], index = (1,1)
			returns 4
			@param index needs to be 2d!
			@returns linear index
		*/
		int convertToLinearIndex(Eigen::Vector2i index);

	private:
		Eigen::Matrix2i grid;
	};
}
#endif