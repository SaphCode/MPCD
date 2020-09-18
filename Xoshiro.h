#pragma once

#ifndef XOSHIRO_H
#define XOSHIRO_H

class Xoshiro
{
public:
	/* Creates a Xoshiro256plusplus generator. Used only for uniformly real distributed doubles.
	@param min minimum of distribution.
	@param max maximum of distribution. */
	Xoshiro(double min, double max);
	Xoshiro();
	~Xoshiro() {}
	double next() const;

	//void setMax(double max);
	//void setMin(double min);

private:
	const double _min;
	const double _max;
};

#endif