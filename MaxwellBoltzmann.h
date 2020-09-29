#pragma once

#ifndef MAXWELL_BOLTZMANN_H
#define MAXWELL_BOLTZMANN_H
#include <Eigen/Dense>
#include <random>
class MaxwellBoltzmann
{
public:
	MaxwellBoltzmann(double mean, double temperature, double mass);
	Eigen::Vector2d next();
private:
	std::mt19937_64 m_gen_x;
	std::normal_distribution<double> m_dist_x;

	std::mt19937_64 m_gen_y;
	std::normal_distribution<double> m_dist_y;

	const double m_temperature;
	const double m_mass;
	const double m_mean;
	const double m_k_boltzmann;
};
#endif // !MAXWELL_BOLTZMANN_H
