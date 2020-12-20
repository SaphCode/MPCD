#pragma once

#ifndef CELL_H
#define CELL_H
#include <Eigen/Dense>
#include "Particle.h"
#include "Thermostat.h"
#include "Monomer.h"
#include <execution>
#include <random>
namespace MPCD {
	class Cell
	{
	public:
		Cell();
		~Cell();
		void add(const Particle& p);
		void add(const Monomer& p);
		void addVirtual(const Particle& p);
		void setOccupied(bool occupied) { 
			m_occupied = occupied;
		}
		bool isOccupied() const {
			return m_occupied;
		}
		void clear();
		void draw(std::pair<int, int> index, std::ofstream& ofs) const;

		double thermostatScaling();

		Eigen::Vector2d getMeanVelocity() const {
			return m_cmVelocity;
		}

		double getRotationAngle() const {
			return m_rotationAngle;
		}

		double getScalingFactor() const {
			return m_scalingFactor;
		}

		size_t number() const {
			return _particles.size();
		}

		void calculate();
	private:
		bool m_occupied = false;
		
		
		std::mt19937_64 _angleGen;
		const std::uniform_real_distribution<double> _unifAngle;
		std::mt19937_64 _signGen;
		const std::uniform_real_distribution<double> _unifSign;

		std::vector<Particle> _particles;
		std::vector<Particle> _virtualParticles;
		std::vector<Monomer> _monomers;

		Thermostat m_thermostat;
		Eigen::Vector2d m_cmVelocity;
		double m_scalingFactor;
		double m_rotationAngle;
	};
}
#endif // !CELL_H
