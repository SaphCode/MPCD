#include "MaxwellBoltzmann.h"
#include "Constants.h"

MaxwellBoltzmann::MaxwellBoltzmann(double mean, double temperature, double mass) :
	m_temperature(temperature),
	m_mass(mass),
	m_mean(mean),
	m_k_boltzmann(MPCD::Constants::k_boltzmann)
{

	std::random_device rd{};
	std::mt19937_64 gen_x{ rd() };
	m_gen_x = gen_x;
	std::mt19937_64 gen_y{ rd() };
	m_gen_y = gen_y;

	double stddev = std::sqrt(m_k_boltzmann * temperature / mass);

	std::normal_distribution<double> dist_x{ m_mean, stddev };
	m_dist_x = dist_x;
	std::normal_distribution<double> dist_y{ m_mean, stddev };
	m_dist_y = dist_y;
}

Eigen::Vector2d MaxwellBoltzmann::next()
{
	double x = m_dist_x(m_gen_x);
	double y = m_dist_y(m_gen_y);
	return Eigen::Vector2d(x, y);
}