#include "Cell.h"

#define _USE_MATH_DEFINES
#include <math.h>

#include <fstream>
#include "Constants.h"

MPCD::Cell::Cell() :
	_angleGen{std::random_device()()},
	_signGen{ std::random_device()() },
	_unifAngle(0, 2 * M_PI),
	_unifSign(-1, 1),
	_vel(0, 0),
	_virtualVel(0, 0),
	_numVirtual(0)
{
}

MPCD::Cell::~Cell() { }

void MPCD::Cell::clear() {
	_particles.clear();
	Eigen::Vector2d zero(0, 0);
	_vel = zero;
	_virtualVel = zero;
	_numVirtual = 0;
}

void MPCD::Cell::add(MPCD::Particle& p) {
	_particles.push_back(std::make_shared<Particle>(p)); // TODO: do we need shared here?
	_vel += p.getVelocity();
}

void MPCD::Cell::addVirtual(MPCD::Particle& p) {
	_numVirtual += 1;
	_virtualVel += p.getVelocity();
}

void MPCD::Cell::collide(double temperatureScalingFactor) {
	double rotationAngle = _unifAngle(_angleGen);
	double s = _unifAngle(_angleGen);
	int sign = (s < 0) ? -1 : 1;
	Eigen::Vector2d mean = getMeanVelocity();
	for (auto& p : _particles) {
		p->collide(mean, sign * rotationAngle, temperatureScalingFactor);
	}
}

Eigen::Vector2d MPCD::Cell::getMeanVelocity() const {
	Eigen::Vector2d mean = (_vel + _virtualVel) / (_particles.size() + _numVirtual);
	return mean;
}

Eigen::Vector2d MPCD::Cell::getTotalCellVelocity() const {
	return _vel;
}

void MPCD::Cell::draw(std::pair<int, int> index, std::ofstream& ofs) const {
	ofs << index.first << "," << index.second << "," << _vel[0] << "," << _vel[1] << "," << _particles.size() << "\n";
}

int MPCD::Cell::number() const
{
	return _particles.size();
}

double MPCD::Cell::thermostatScaling() {
	Eigen::Vector2d meanCellVelocity = getMeanVelocity();
	return m_thermostat.getScalingFactor(_particles, meanCellVelocity);
}