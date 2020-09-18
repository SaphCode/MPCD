#include "Xoshiro.h"
#include "xoshiro256plusplus.h"

Xoshiro::Xoshiro(const double min, const double max) {
	_min = min;
	_max = max;
	x_initialize();
}

Xoshiro::Xoshiro() {
	_max = 1.0;
	_min = 0.0;
	x_initialize();
}

/*
void Xoshiro::setMin(double min) {
	_min = min;
}
*/

/*
void Xoshiro::setMax(double max) {
	_max = max;
}
*/

double Xoshiro::next() const {
	return x_getDouble(_min, _max);
}
