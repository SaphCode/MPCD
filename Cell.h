#pragma once

#ifndef CELL_H
#define CELL_H
#include <Eigen/Dense>
#include "Particle.h"
#include "Xoshiro.h"
#include <execution>
namespace MPCD {
	class Cell
	{
	public:
		Cell();
		~Cell();
		void add(Particle& p);
		//void remove(Particle& p);
		void collide();
		void clear();
		//friend Cell operator+(const Cell& lhs, const Cell& rhs);
		friend std::ostream& operator<<(std::ostream& os, Cell const& c) {
			os << "Vel: " << std::to_string(c._vel[0]) << ", Num: " << std::to_string(c._num);
			return os;
		}
		void draw(std::mutex& m, std::pair<int, int> index, std::ofstream& ofs);
	private:
		std::vector<Particle> _particles;
		Eigen::Vector2d _vel;
		int _num;
		Xoshiro _angleGen;
		Xoshiro _signGen;
	};
}
#endif // !CELL_H
