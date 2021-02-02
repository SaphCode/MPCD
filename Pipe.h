#pragma once

#include <Eigen/Dense>
#include "Particle.h"
#include "Constants.h"
#include "Wall.h"
#include "CircularObstacle.h"
#include "ConstForce.h"
#include "Monomer.h"

namespace MPCD {
	class Pipe
	{
	private:
		std::vector<CircularObstacle> m_obstacles;
		std::vector<Wall> m_walls;
		ConstForce m_constForce;
		double _x_0 = Constants::x_0;
		double _x_max = Constants::x_max;
		double _y_0 = Constants::y_0;
		double _y_max = Constants::y_max;
		void collide(Body& p);
		void fixOutOfBounds(Body& p);
		bool inBounds(const Eigen::Vector2d& pos);
		void calculateInteraction(int currentIndex, Monomer& m, std::vector<Monomer>& monomers);
		void verletVelocity(std::vector<Monomer>& monomers);
		void verletPosition(std::vector<Monomer>& monomers);


		int _w;
	public:
		Pipe(ConstForce force);
		~Pipe() {}
		void stream(std::vector<Particle>& particles, double lapse, bool draw, int t);
		void setObstacles(std::vector<CircularObstacle> obstacles, std::vector<Wall> walls);
		
		void verlet(std::vector<Monomer>& monomers, bool draw, int t);
	};
}
