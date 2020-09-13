#include "MaxwellBoltzmann.h"

MaxwellBoltzmann::MaxwellBoltzmann(double mean, double temperature, double mass)
{
	_temperature = temperature;
	_mass = mass;
	_mean = mean;

	std::random_device rd{};
	std::mt19937_64 gen_x{ rd() };
	_gen_x = gen_x;
	std::mt19937_64 gen_y{ rd() };
	_gen_y = gen_y;

	double stddev = std::sqrt(_k_boltzmann * temperature / mass);

	std::normal_distribution<double> dist_x{ _mean, stddev };
	_dist_x = dist_x;
	std::normal_distribution<double> dist_y{ _mean, stddev };
	_dist_y = dist_y;
}

Eigen::Vector2d MaxwellBoltzmann::next()
{
	double x = _dist_x(_gen_x);
	double y = _dist_y(_gen_y);
	return Eigen::Vector2d(x, y);
}