#pragma once

#ifndef GRID_H
#define GRID_H

#include <Eigen/Dense>
#include "Cell.h"
#include "Particle.h"
#include <map>
#include <mutex>
#include "Constants.h"
#include "CircularObstacle.h"
#include "Wall.h"

namespace MPCD {
	class Grid {
		public:
			Grid();
			void setupCells(std::vector<CircularObstacle> obstacles, std::vector<Wall> walls);
			void updateCoordinates(std::vector<Particle>& particles);
			void collision(bool draw, std::ofstream& file);
			void shift();
			void undoShift();
			int getAverageParticlesPerCell() const;
			double getA() const;
			int getNumRows() const;
			int getNumCols() const;
			double getMaxShift() const;
		private:
			void createVirtualParticles(const std::pair<int, int>& key, Cell& cell, const double cell_dim);
			std::pair<int, int> getCoordinates(Eigen::Vector2d position) const;
			std::mt19937_64 _shiftGen;
			const std::uniform_real_distribution<double> _unif;
			Eigen::Vector2d _shift;
			const int _average_particles_per_cell;
			const double _a;
			const int _numRows;
			const int _numCols;
			std::map<std::pair<int, int>, Cell> _cells;
			const double _maxShift;
	};
}
#endif
