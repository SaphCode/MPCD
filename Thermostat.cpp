#include "Thermostat.h"
#include <cmath>
#include <algorithm>
#include "Force.h"

using namespace Eigen;

MPCD::Thermostat::Thermostat(double particleMass, double k_BT) :
    m_k_BT_0(k_BT),
    m_particleMass(particleMass),
	m_gen{ std::random_device()() }
{

}

double MPCD::Thermostat::getScalingFactor(const std::vector<std::shared_ptr<Particle>>& particles, const Eigen::Vector2d cellMeanVelocity) const
{
	const double degFreedom = 2 * (particles.size() - 1);
	std::gamma_distribution<double> dist(degFreedom / 2, m_k_BT_0);
	double kineticEnergy = dist(m_gen);
	double sum = calculateSum(particles, cellMeanVelocity);
	double alpha = std::sqrt((2 * kineticEnergy) / (m_particleMass * sum));
	return alpha;
}

double MPCD::Thermostat::calculateSum(const std::vector<std::shared_ptr<Particle>>& particles, const Eigen::Vector2d cellMeanVelocity) const
{
	double sum = 0;
	for (const auto& p : particles) {
		Vector2d pVel = p->getVelocity();
		Vector2d diff = (pVel - cellMeanVelocity);
		sum += diff.dot(diff);
	}
	return sum;
}
