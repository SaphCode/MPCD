#include "GammaDistribution.h"
#include <assert.h>

GammaDistribution::GammaDistribution(const double particleMass, const double beta) :
	m_particleMass(particleMass),
	m_beta(beta)
{
	
}

GammaDistribution::~GammaDistribution()
{
}

double GammaDistribution::next(const int numParticles, const double alpha)
{
	auto search = m_distributions.find(numParticles);
	double next = -1;
	if ((search != m_distributions.end())) {
		next = search->second.second(search->first);
	}
	else {
		std::mt19937_64 gen{ std::random_device()() };
		std::gamma_distribution<double> dist(alpha, m_beta);
		next = dist(gen);
		m_distributions.emplace(numParticles, std::make_pair(gen, dist));
	}
	return next;
}
