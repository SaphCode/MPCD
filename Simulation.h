#pragma once

#ifndef MPCD_H
#define MPCD_H

#include "Particle.h"
#include <Eigen/Dense>
#include "Grid.h"
#include "Pipe.h"
#include "InteractingBody.h"
#include "IObstacle.h"

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
		std::vector<std::shared_ptr<IObstacle>> _obstacles;
		std::vector<std::shared_ptr<InteractingBody>> _interactors;

		void setup();

		void streamingStep();
		void collisionStep();

		bool isInBoundsOfAnObstacle(Body& b);

		void writeConstantsToOut(double timelapse, double width, double height, double cell_dim, int averageParticlesPerCell, int timesteps);
		
		void setUpParticles(int number, double x_0, double x_max, double y_0, double y_max);

	public:
		Simulation(bool draw);
		~Simulation() {}
		/* One timestep */
		void timestep();
	};
	
}

#endif // !MPCD_H
