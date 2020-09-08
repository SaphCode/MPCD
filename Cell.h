#pragma once

#ifndef CELL_H
#define CELL_H
#include <Eigen/Dense>
#include "Particle.h"
#include "Xoshiro.h"
namespace MPCD {
	class Cell
	{
	public:
		Cell();
		~Cell();
		void add(Particle& p);
		//void remove(Particle& p);
		void collide(Eigen::Vector2d shift);
		friend Cell operator+(const Cell& lhs, const Cell& rhs);
	private:
		std::vector<Particle> _particles;
		Eigen::Vector2d _vel;
		int _num;
		Xoshiro _angleGen;
		Xoshiro _signGen;
	};
}
#endif // !CELL_H
