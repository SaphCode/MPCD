#pragma once

#ifndef THERMOSTAT_H
#define THERMOSTAT_H
#include "Xoshiro.h"
#include <Eigen/Dense>
#include "Particle.h"
#include "Constants.h"
#include <random>
#include "GammaDistribution.h"

namespace MPCD {
	class Thermostat
	{

	public:
		Thermostat(const GammaDistribution& gamma, double particleMass, double k_BT);
		double getScalingFactor(const std::vector<std::shared_ptr<Particle>>& particles, const Eigen::Vector2d cellMeanVelocity);
	private:
		double calculateSum(const std::vector<std::shared_ptr<Particle>>& particles, const Eigen::Vector2d cellMeanVelocity) const;
		const double m_k_BT_0;
		const double m_particleMass = MPCD::Constants::particle_mass;
		GammaDistribution m_gammaDist;
	};

}

#endif // !THERMOSTAT_H