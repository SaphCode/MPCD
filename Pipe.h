#pragma once

#ifndef PIPE_H
#define PIPE_H

#include <Eigen/Dense>
#include "Particle.h"
#include "IObstacle.h"
#include "Constants.h"

namespace MPCD {
	class Pipe
	{
	private:
		std::vector<std::shared_ptr<IObstacle>> _obstacles;
		double _x_0 = Constants::x_0;
		double _x_max = Constants::x_max;
		double _y_0 = Constants::y_0;
		double _y_max = Constants::y_max;
		void collide(Particle& p);
		void fixOutOfBounds(Particle& p);
		bool inBounds(const Eigen::Vector2d& pos);
		
	public:
		Pipe();
		~Pipe() {}
		void stream(std::vector<Particle>& particles, std::vector<std::shared_ptr<InteractingBody>>& interactors, double lapse, bool draw, std::ofstream& file);
		//const std::vector<MPCD::Particle>& getParticles();
		//void setParticles(std::vector<Particle>& particles);
		void setObstacles(std::vector<std::shared_ptr<IObstacle>>& obstacles);
	};
}

#endif // !PIPE_H
