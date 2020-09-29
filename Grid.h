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
			int getAverageParticlesPerCell() const;
			double getA() const;
			int getNumRows() const;
			int getNumCols() const;
			double getMaxShift() const;
		private:
			void createVirtualParticles(const std::pair<int, int>& key, Cell& cell, const int firstRow, const int lastRow, const double cell_dim);
			std::pair<int, int> getCoordinates(Eigen::Vector2d position) const;
			const Xoshiro _shiftGen;
			Eigen::Vector2d _shift;
			const int _average_particles_per_cell;
			const double _a;
			const int _numRows;
			const int _numCols;
			//int _num_cells;
			std::map<std::pair<int, int>, Cell> _cells;
			//constexpr int wanted_num_cells = min_num_cells * 2;
			const double _maxShift;
	};
}
#endif
