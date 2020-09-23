#include "Thermostat.h"
#include <cmath>
#include <algorithm>
#include "Force.h"

using namespace Eigen;

MPCD::Thermostat::Thermostat() : _scalingFactorGen(1, 1 + _c), _50percentGen(-1, 1), _doWeScaleGen(0, 1) {
}

double MPCD::Thermostat::getScalingFactor(const std::vector<Particle>& particles, const int dimension) const
{
    double S = _scalingFactorGen.next();
    double p = _50percentGen.next();
    if (p < 0) {
        S = 1 / S;
    }
    double ksi = _doWeScaleGen.next();
    double A = calculateA(S, particles, dimension);
    double acceptanceProbabilty = std::min(1.0, A);
    if (ksi < acceptanceProbabilty) {
        return S;
    }
    else {
        return 1;
    }
}

double MPCD::Thermostat::calculateSum(const double scalingFactor, const std::vector<Particle>& particles) const {
    double sum = 0;
    for (const auto& p : particles) {
        Vector2d pVel = p.getVelocity();
        Vector2d u = _flowProfile(p.getPosition());
        Vector2d relVel = pVel - u;
        sum += relVel.dot(relVel) * (scalingFactor * scalingFactor - 1);
    }
    return sum;
}

double MPCD::Thermostat::calculateA(double scalingFactor, const std::vector<Particle>& particles, int dimension) const
{
    double sum = calculateSum(scalingFactor, particles);
    double A = std::pow(scalingFactor, dimension * (particles.size() - 1)) * std::exp(-_particleMass / (2 * _boltzmannConst * _T0) * sum);
    return A;
}

Eigen::Vector2d MPCD::Thermostat::_flowProfile(Eigen::Vector2d pos) const
{
    double R = (MPCD::Constants::y_max - MPCD::Constants::y_0)/2;
    double r = pos[1] - R; // TODO
    double L = MPCD::Constants::x_max - MPCD::Constants::x_0;
    double pressureDifference = Force::const_force;
    double vx = pressureDifference * (R * R - r * r) / (4 * MPCD::Constants::viscosity * L);
    Vector2d flow(vx, 0);
    return flow;
}
