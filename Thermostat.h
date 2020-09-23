#pragma once

#ifndef THERMOSTAT_H
#define THERMOSTAT_H
#include "Xoshiro.h"
#include <Eigen/Dense>
#include "Particle.h"
#include "Constants.h"

namespace MPCD {
	class Thermostat
	{

	public:
		Thermostat();
		double getScalingFactor(const std::vector<Particle>& particles, const int dimension) const;
	private:
		double calculateSum(const double scalingFactor, const std::vector<Particle>& particles) const;
		double calculateA(const double scalingFactor, const std::vector<Particle>& particles, const int dimension) const;
		const double _c = 0.28;
		const double _T0 = MPCD::Constants::temperature;
		const double _particleMass = MPCD::Constants::particle_mass;
		const double _boltzmannConst = MPCD::Constants::k_boltzmann;
		const Xoshiro _scalingFactorGen;
		const Xoshiro _50percentGen;
		const Xoshiro _doWeScaleGen;
		Eigen::Vector2d _flowProfile(Eigen::Vector2d pos) const;
	};

}

#endif // !THERMOSTAT_H