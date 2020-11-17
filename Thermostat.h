#pragma once

#ifndef THERMOSTAT_H
#define THERMOSTAT_H
#include <Eigen/Dense>
#include "Particle.h"
#include "Constants.h"
#include <random>

namespace MPCD {
	class Thermostat
	{

	public:
		Thermostat();
		double getScalingFactor(const std::vector<Particle>& particles, const Eigen::Vector2d cellMeanVelocity);
	private:
		double calculateSum(const std::vector<Particle>& particles, const Eigen::Vector2d cellMeanVelocity) const;
		std::mt19937_64 _gen{ std::random_device()() };
	};

}

#endif // !THERMOSTAT_H