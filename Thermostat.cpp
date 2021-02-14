#include "Thermostat.h"
#include <cmath>
#include <algorithm>
#include "Force.h"
#include <iostream>

using namespace Eigen;

MPCD::Thermostat::Thermostat(){}

double MPCD::Thermostat::getScalingFactor(const std::vector<Particle>& particles, const Eigen::Vector2d cellMeanVelocity)
{
	const size_t particleNum = particles.size();
	double alpha = 1;
	if (particleNum > 1) {
		const int df = 2*(particleNum - 1);
		std::gamma_distribution<double> dist(df / 2.0, MPCD::Constants::k_boltzmann * MPCD::Constants::temperature);
		double kineticEnergy = dist(_gen);
		double sum = calculateSum(particles, cellMeanVelocity);
		alpha = std::sqrt((2 * kineticEnergy) / sum);
	}
	return alpha;
} 

double MPCD::Thermostat::calculateSum(const std::vector<Particle>& particles, const Eigen::Vector2d cellMeanVelocity) const
{
	double sum = 0;
	for (const auto& p : particles) {
		Vector2d pVel = p.getVelocity();
		Vector2d diff = (pVel - cellMeanVelocity);
		sum += p.getMass() * diff.dot(diff);
	}
	return sum;
}
