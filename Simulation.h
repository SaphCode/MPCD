#pragma once

#ifndef MPCD_H
#define MPCD_H

#include "Particle.h"
#include <Eigen/Dense>
#include "Grid.h"
#include "Pipe.h"
#include "PhysicalObject.h"
#include "ImmovableObstacle.h"

#include <map>

namespace MPCD {
	class Simulation {
	private:
		bool _draw = false;
		double _timelapse;

		Grid _grid;
		Pipe _pipe;
		int _w;

		int _t;
		std::vector<Particle> _particles;
		std::vector<Obstacle> _obstacles;

		std::vector<ImmovableObstacle> _walls;

		PhysicalObject _imaginaryAttractorWall;

		void streamingStep();
		void collisionStep();

		void writeConstantsToOut(double timelapse, double width, double height, double cell_dim, int averageParticlesPerCell, int timesteps);
		
	public:
		Simulation(bool draw);
		~Simulation() {}
		/* One timestep */
		void timestep();
	};
	
}

#endif // !MPCD_H
