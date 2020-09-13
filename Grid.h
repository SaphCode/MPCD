#pragma once

#ifndef GRID_H
#define GRID_H

#include <Eigen/Dense>
#include "Cell.h"
#include "Particle.h"
#include <map>
#include <mutex>
#include "Constants.h"

namespace MPCD {
	class Grid {
		public:
			Grid();
			//void updateCell(Particle p, Eigen::Vector2d positionBefore);
			//void shift();
			//void shiftBack();
			void updateCoordinates(std::vector<Particle>& particles);
			void collision(bool draw, std::ofstream& file);
			//void insert(Particle p);
			void shift();
			void undoShift();
			int getAverageParticlesPerCell();
			double getA();
			int getNumRows();
			int getNumCols();
			double getMaxShift();
		private:
			std::pair<int, int> getCoordinates(Eigen::Vector2d position);
			Xoshiro _shiftGen;
			Eigen::Vector2d _shift;
			const int _average_particles_per_cell = 10;
			const double _a = 1;
			const int _num_rows = std::floor((MPCD::Constants::y_max - MPCD::Constants::y_0) / _a);
			const int _num_cols = std::floor((MPCD::Constants::x_max - MPCD::Constants::x_0) / _a);
			//int _num_cells;
			std::map<std::pair<int, int>, Cell> _cells;
			//constexpr int wanted_num_cells = min_num_cells * 2;
			const double _max_shift = _a / 2;
	};
}
#endif
