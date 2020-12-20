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
#include <fstream>

namespace MPCD {
	class Simulation {
	private:
		const bool _draw;
		double _timelapse;

		Grid _grid;
		Pipe _pipe;
		

		bool _addObstacles = false;

		int _t;
		int _drawInterval;
		int _drawLast;

		std::vector<Particle> _particles;
		std::vector<Monomer> _monomers;
		//std::vector<std::shared_ptr<IObstacle>> _obstacles; // TODO: remove from simulation?
		// idea: make vectors of Wall and Circular obstacles separately,
		// compute interactions and all that stuff separately
		//std::vector<std::shared_ptr<InteractingBody>> _interactors; // TODO: remove from simulation?

		void writeCirclePositionToOut(std::ofstream& outFile, Eigen::Vector2d center_pos, double radius);

		void setup();

		void streamingStep();
		void collisionStep();
		void verlet();

		bool isInBoundsOfAnObstacle(Body& b, std::vector<CircularObstacle> obstacles);

		void writeConstantsToOut(double timelapse, double width, double height, double cell_dim, int averageParticlesPerCell, int timesteps);
		
		void setUpParticles(int number, double x_0, double x_max, double y_0, double y_max, std::vector<CircularObstacle> obstacles);

		void setUpMonomers();

	public:
		Simulation(bool draw, int drawInterval, int drawLast);
		~Simulation() {}
		/* One timestep */
		void timestep();
	};
	
}

#endif // !MPCD_H
