#include "Cell.h"

#define _USE_MATH_DEFINES
#include <math.h>

MPCD::Cell::Cell() {
	Xoshiro angle(0, 2 * M_PI);
	_angleGen = angle;
	Xoshiro sign(-1, 1);
	_signGen = sign;
	Eigen::Vector2d zero(0, 0);
	_vel = zero;
	_num = 0;
}

MPCD::Cell::~Cell() { }

void MPCD::Cell::add(MPCD::Particle& p) {
	_particles.push_back(p);
	_vel += p.getVelocity();
	_num += 1;
}

/*
void MPCD::Cell::remove(MPCD::Particle& p) {
	_particles.erase(std::remove(_particles.begin(), _particles.end(), p), _particles.end());
	_vel -= p.getVelocity();
	_num -= 1;
}
*/

void MPCD::Cell::collide(Eigen::Vector2d shift) {
	double rotationAngle = _angleGen.next();
	double s = _signGen.next();
	int sign = (s < 0) ? -1 : 1;
	Eigen::Vector2d mean = _vel / _num;
	for (auto p : _particles) {
		p.updateVelocity(mean, sign * rotationAngle);
		p.shift(-shift);
	}
}



MPCD::Cell MPCD::operator+(const Cell& lhs, const Cell& rhs)
{
	Cell cell;
	cell._vel = lhs._vel + rhs._vel;
	cell._num = lhs._num + rhs._num;
	cell._particles.reserve(lhs._particles.size() + rhs._particles.size());
	for (auto const& p : lhs._particles) {
		cell._particles.push_back(p);
	}
	for (auto const& p : rhs._particles) {
		cell._particles.push_back(p);
	}
	return cell;
}
