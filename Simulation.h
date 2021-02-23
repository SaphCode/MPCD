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
		const bool _drawParticles;
		double _timestep;
		double _mdTimestep;

		Grid _grid;
		Pipe _pipe;
		

		const bool _addObstacles = true;

		int _t;
		const int _stationaryT;

		std::vector<Particle> _particles;
		std::vector<Monomer> _monomers;
		//std::vector<std::shared_ptr<IObstacle>> _obstacles; // TODO: remove from simulation?
		// idea: make vectors of Wall and Circular obstacles separately,
		// compute interactions and all that stuff separately
		//std::vector<std::shared_ptr<InteractingBody>> _interactors; // TODO: remove from simulation?

		void writeCirclePositionToOut(std::ofstream& outFile, Eigen::Vector2d center_pos, double radius);

		std::vector<CircularObstacle>& setupObstacles(std::ofstream& outFile);

		void streamingStep();
		void collisionStep();
		void verlet();

		bool isInBoundsOfAnObstacle(Body& b, std::vector<CircularObstacle> obstacles);

		void writeConstantsToOut(double timelapse, double width, double height, double cell_dim, int averageParticlesPerCell);
		
		void setUpParticles(int number, double x_0, double x_max, double y_0, double y_max, std::vector<CircularObstacle> obstacles);

		void setUpMonomers();

	public:
		Simulation(bool draw, bool particleDrawing, int stationaryT);
		~Simulation() {}
		/* One timestep */
		void timestep();

		void setup();

		void loadCheckpoint(int timestep, std::string pathParticles);

	};
	
}

#endif // !MPCD_H
