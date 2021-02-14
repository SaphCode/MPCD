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
		void collide(Body& p);
		void fixOutOfBounds(Body& p);
		bool inBounds(const Eigen::Vector2d& pos);
		void calculateInteraction(int chainIndex, std::vector<Monomer>& monomers);
		void verletVelocity(int chainIndex, std::vector<Monomer>& monomers, double timestep);
		void verletPosition(int chainIndex, std::vector<Monomer>& monomers, double timestep);

		const bool isLegalPosition(const Body& b);

		int _w;
	public:
		Pipe(ConstForce force);
		~Pipe() {}
		void stream(std::vector<Particle>& particles, double lapse, bool draw, int t);
		void setObstacles(std::vector<CircularObstacle> obstacles, std::vector<Wall> walls);
		
		void verlet(std::vector<Monomer>& monomers, bool draw, int t);
	};
}
