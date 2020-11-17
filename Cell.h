#pragma once

#ifndef CELL_H
#define CELL_H
#include <Eigen/Dense>
#include "Particle.h"
#include "Thermostat.h"
#include <execution>
#include <random>
namespace MPCD {
	class Cell
	{
	public:
		Cell();
		~Cell();
		void add(Particle& p);
		void addVirtual(Particle& p);
		void setOccupied(bool occupied) { 
			m_occupied = occupied;
		}
		bool isOccupied() const {
			return m_occupied;
		}
		//void remove(Particle& p);
		void collide(double temperatureScalingFactor);
		void clear();
		//friend Cell operator+(const Cell& lhs, const Cell& rhs);
		friend std::ostream& operator<<(std::ostream& os, Cell const& c) {
			os << "Vel: " << std::to_string(c._vel[0]) << ", Num: " << std::to_string(c._particles.size());
			return os;
		}
		void draw(std::pair<int, int> index, std::ofstream& ofs) const;
		int number() const;

		double thermostatScaling();

		Eigen::Vector2d getTotalCellVelocity() const;
	private:
		bool m_occupied = false;
		Eigen::Vector2d getMeanVelocity() const;
		std::vector<std::shared_ptr<Particle>> _particles;
		Eigen::Vector2d _vel;
		Eigen::Vector2d _virtualVel;
		int _numVirtual;
		std::mt19937_64 _angleGen;
		const std::uniform_real_distribution<double> _unifAngle;
		std::mt19937_64 _signGen;
		const std::uniform_real_distribution<double> _unifSign;
		Thermostat m_thermostat;
	};
}
#endif // !CELL_H
