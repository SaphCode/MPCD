#pragma once

#include <random>
#include <map>

class GammaDistribution
{

public:
	GammaDistribution(double particleMass, double beta);
	~GammaDistribution();

	double next(int numParticles, double alpha);

private:
	const double m_particleMass;
	const double m_beta;
	std::map<int, std::pair<std::mt19937_64, std::gamma_distribution<double>>> m_distributions;

};

