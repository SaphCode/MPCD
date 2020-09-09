#pragma once

#ifndef GRID_H
#define GRID_H

#include <Eigen/Dense>
#include "Cell.h"
#include "Particle.h"
#include <map>
#include <mutex>

namespace MPCD {
	class Grid {
	public:
		Grid();
		//void updateCell(Particle p, Eigen::Vector2d positionBefore);
		//void shift();
		//void shiftBack();
		friend Grid operator+(const Grid& lhs, const Grid& rhs);
		std::map<std::pair<int, int>, Cell> _cells;
		std::pair<int, int> getCoordinates(Eigen::Vector2d position);
		void insert(Particle p);
	private:
		double _a;
	};
}
#endif
