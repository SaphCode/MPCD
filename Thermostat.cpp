#include "Thermostat.h"
#include <cmath>
#include <algorithm>
#include "Force.h"

using namespace Eigen;

MPCD::Thermostat::Thermostat(const GammaDistribution& gamma, double particleMass, double k_BT) :
    m_k_BT_0(k_BT),
    m_particleMass(particleMass),
	m_gammaDist(gamma)
{
}

double MPCD::Thermostat::getScalingFactor(const std::vector<std::shared_ptr<Particle>>& particles, const Eigen::Vector2d cellMeanVelocity)
{
	const int particleNum = particles.size();
	double alpha = 1;
	if (particleNum > 1) {
		const int gammadist_alpha = (particleNum - 1);
		double kineticEnergy = m_gammaDist.next(particleNum, gammadist_alpha);
		double sum = calculateSum(particles, cellMeanVelocity);
		alpha = std::sqrt((2 * kineticEnergy) / (m_particleMass * sum));
	}	
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
