#include "Cell.h"

#define _USE_MATH_DEFINES
#include <math.h>

#include <fstream>
#include "Constants.h"
#include <iostream>

MPCD::Cell::Cell() :
	_angleGen{ std::random_device()() },
	_unifAngle(0, 2 * M_PI),
	m_scalingFactor(1),
	m_rotationAngle(0)
{
	_particles.reserve(MPCD::Constants::average_particles_per_cell);
}

MPCD::Cell::~Cell() { }

void MPCD::Cell::clear() {
	_particles.clear();
	_particles.reserve(MPCD::Constants::average_particles_per_cell);
	_virtualParticles.clear();
	_monomers.clear();
	m_cmVelocity = Eigen::Vector2d(0, 0);
	m_scalingFactor = 1;
	m_rotationAngle = 0;
}

void MPCD::Cell::add(const MPCD::Particle& p) {
	_particles.push_back(p); // TODO: do we need shared here?
}

void MPCD::Cell::add(const Monomer& m)
{
	_monomers.push_back(m);
}

void MPCD::Cell::addVirtual(const MPCD::Particle& p) {
	_virtualParticles.push_back(p);
}

void MPCD::Cell::calculate()
{
	if (_particles.size() > 0 || _monomers.size() > 0) {
		Eigen::Vector2d momentum(0, 0);
		double mass = MPCD::Constants::particle_mass;
		for (const auto& p : _particles) {
			momentum += p.getVelocity() * p.getMass();
			assert(mass == p.getMass());
		}
		double virtualMass = MPCD::Constants::particle_mass;
		Eigen::Vector2d virtualMomentum(0, 0);
		for (const auto& vp : _virtualParticles) {
			virtualMomentum += vp.getVelocity() * vp.getMass();
			assert(virtualMass == vp.getMass());
		}

		Eigen::Vector2d monomerMomentum(0, 0);
		double monomerMass = MPCD::Constants::monomer_mass;
		for (const auto& m : _monomers) {
			monomerMomentum += m.getVelocity() * m.getMass();
			assert(monomerMass == m.getMass());
		}

		m_cmVelocity = (momentum + virtualMomentum + monomerMomentum) / (mass * _particles.size() + virtualMass * _virtualParticles.size() + monomerMass * _monomers.size());
		m_scalingFactor = thermostatScaling();
		m_rotationAngle = _unifAngle(_angleGen);
	}
}

void MPCD::Cell::draw(std::pair<int, int> index, std::ofstream& ofs) const {
	// if occupied set true
	ofs << index.first << "," << index.second << "," << m_cmVelocity[0] << "," << m_cmVelocity[1] << "," << _particles.size() << "\n";
}

double MPCD::Cell::thermostatScaling() {
	return m_thermostat.getScalingFactor(_particles, m_cmVelocity);
}