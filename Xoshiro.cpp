#include "Xoshiro.h"
#include "xoshiro256plusplus.h"

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
