#include "Thermostat.h"
#include <cmath>
#include <algorithm>

using namespace Eigen;

MPCD::Thermostat::Thermostat() {
    Xoshiro scalingFactorGen(1, 1 + _c);
    _scalingFactorGen = scalingFactorGen;
    Xoshiro fiftyPercentGen(-1, 1);
    _50percentGen = fiftyPercentGen;
    Xoshiro doWeScaleGen(0, 1);
    _doWeScaleGen = doWeScaleGen;

    // _u = ; TODO
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
    double acceptanceProbabilty = std::min(1, A);
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
        Vector2d relVel = pVel - _u;
        sum += relVel.dot(relVel) * (scalingFactor * scalingFactor - 1);
    }
    return sum;
}

double MPCD::Thermostat::calculateA(double scalingFactor, const std::vector<Particle>& particles, int dimension) const
{
    Vector2d u = getFlowProfile();
    double sum = calculateSum(scalingFactor, particles);
    double A = std::pow(scalingFactor, dimension * (particles.size() - 1)) * std::exp(-_particleMass / (2 * _boltzmannConst * _T0) * sum);
}
