#include "MersenneTwister.h"
#include <stdexcept>
#include <random>

MersenneTwister::MersenneTwister(const double min, const double max, DistributionType type) {
	_min = min;
	_max = max;
	_type = type;
	std::mt19937 generator(_rand_dev());
	_gen = generator;
	switch (_type) {
	case DistributionType::UNIFORM:
	{
		std::uniform_real_distribution<double> distribution(_min, _max);
		_dist = distribution;
		break;
	}
	case DistributionType::MAXWELL:
	{
		throw std::exception("Not implemented yet.");
		break;
	}
	default:
	{
		throw std::exception("No valid type chosen.");
	}
	}
}

double MersenneTwister::next() {
	double d = _dist(_gen);
	_dist.reset();
	return d;
}
