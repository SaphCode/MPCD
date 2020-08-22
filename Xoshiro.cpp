#include "Xoshiro.h"
#include "xoshiro256plusplus.h"

Xoshiro::Xoshiro(double min, double max) {
	_min = min;
	_max = max;
	x_initialize();
}

double Xoshiro::next() {
	return x_getDouble(_min, _max);
}
