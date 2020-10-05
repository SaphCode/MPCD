#pragma once

#ifndef CELL_H
#define CELL_H
#include <Eigen/Dense>
#include "Particle.h"
#include "Xoshiro.h"
#include "Thermostat.h"
#include <execution>
namespace MPCD {
	class Cell
	{
	public:
		Cell(const GammaDistribution& gamma);
		~Cell();
		void add(Particle& p);
		void addVirtual(Particle& p);
		//void remove(Particle& p);
		void collide(double temperatureScalingFactor);
		void clear();
		//friend Cell operator+(const Cell& lhs, const Cell& rhs);
		friend std::ostream& operator<<(std::ostream& os, Cell const& c) {
			os << "Vel: " << std::to_string(c._vel[0]) << ", Num: " << std::to_string(c._particles.size());
			return os;
		}
		void draw(std::mutex& m, std::pair<int, int> index, std::ofstream& ofs) const;
		int number() const;

		double thermostatScaling();

		Eigen::Vector2d getTotalCellVelocity() const;
	private:
		Eigen::Vector2d getMeanVelocity() const;
		std::vector<std::shared_ptr<Particle>> _particles;
		Eigen::Vector2d _vel;
		Eigen::Vector2d _virtualVel;
		int _numVirtual;
		const Xoshiro _angleGen;
		const Xoshiro _signGen;
		Thermostat m_thermostat;
	};
}
#endif // !CELL_H
