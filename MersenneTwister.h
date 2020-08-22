#pragma once
#include <random>

#ifndef MERSENNETWISTER_H
#define MERSENNETWISTER_H

enum class DistributionType {
	UNIFORM,
	MAXWELL,
};

class MersenneTwister
{
public:
	/* Creates a MersenneTwister.
	@param min minimum of distribution.
	@param max maximum of distribution.
	*/
	MersenneTwister(double min, double max, DistributionType type);
	~MersenneTwister() {}
	double next();

private:
	DistributionType _type;
	double _min;
	double _max;
	std::uniform_real_distribution<double> _dist;
	std::random_device _rand_dev;
	std::mt19937 _gen;
};


#endif