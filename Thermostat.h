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
		const double _c;
		const double _T0;
		const double _particleMass;
		const double _boltzmannConst;
		Eigen::Vector2d _u;
		const Xoshiro _scalingFactorGen;
		const Xoshiro _50percentGen;
		const Xoshiro _doWeScaleGen;
	};

}

#endif // !THERMOSTAT_H