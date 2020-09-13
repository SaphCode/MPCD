#pragma once

#ifndef PIPE_H
#define PIPE_H

#include <Eigen/Dense>
#include "Particle.h"
#include "Obstacle.h"
#include "Constants.h"

namespace MPCD {
	class Pipe
	{
	private:
		std::vector<MPCD::Particle> _particles;
		std::vector<MPCD::Obstacle> _obstacles;
		double _x_0 = MPCD::Constants::x_0;
		double _x_max = MPCD::Constants::x_max;
		double _y_0 = MPCD::Constants::y_0;
		double _y_max = MPCD::Constants::y_max;
		void collide(Particle p);
		void fixOutOfBounds(Particle p);
		bool inBounds(Eigen::Vector2d pos);
	public:
		Pipe();
		~Pipe() {}
		void stream(double lapse, bool draw, std::ofstream& file);
		std::vector<MPCD::Particle>& getParticles();
		void setParticles(std::vector<Particle> particles);
		void setObstacles(std::vector<Obstacle> obstacles);
	};
}

#endif // !PIPE_H
