#pragma once

#ifndef MPCD_H
#define MPCD_H

#include "Particle.h"
#include <Eigen/Dense>
#include "Grid.h"
#include "Pipe.h"

#include <map>

namespace MPCD {
	class Simulation {
	private:
		bool _draw = false;
		double _timelapse;

		Grid _grid;
		Pipe _pipe;
		int _w;

		Xoshiro _rgShiftX;
		Xoshiro _rgShiftY;
		double _maxShift;
		int _t;

		void streamingStep(Eigen::Vector2d shift);
		void collisionStep(Eigen::Vector2d shift);
		
	public:
		Simulation(bool draw);
		~Simulation() {}
		/* One timestep */
		void timestep();
	};

	void calculateGrid(Grid& g, std::vector<Particle>& particles, const Eigen::Vector2d shift, const double timelapse);
	
}

#endif // !MPCD_H
