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
	std::mt19937_64 _gen_x;
	std::normal_distribution<double> _dist_x;

	std::mt19937_64 _gen_y;
	std::normal_distribution<double> _dist_y;

	double _temperature;
	double _mass;
	double _mean;
	const double _k_boltzmann = 1.38064852e-23;
};
#endif // !MAXWELL_BOLTZMANN_H
