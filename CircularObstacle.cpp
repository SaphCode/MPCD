#include "CircularObstacle.h"
#include "Constants.h"

MPCD::CircularObstacle::CircularObstacle(Eigen::Vector2d center, double radius) :
	InteractingBody(std::numeric_limits<double>::infinity(), center, Eigen::Vector2d(0, 0)),
	m_center(center),
	m_radius(radius)
{

}

bool MPCD::CircularObstacle::isInBounds(const Body& o) const
{
	Eigen::Vector2d rel = o.getPosition() - m_center;
	if (rel.stableNorm() <= m_radius) {
		return true;
	}
	return false;
}

Eigen::Vector2d MPCD::CircularObstacle::getOvershoot(const Body& o) const
{
	Eigen::Vector2d oldPos = o.getOldPosition(MPCD::Constants::time_lapse);
	Eigen::Vector2d pos = o.getPosition();
	Eigen::Vector2d relPos = oldPos - pos;

	int sign;
	if ((oldPos - pos)[1] >= 0) {
		sign = 1;
	}
	else {
		sign = -1;
	}

	/*

	double k2 = std::pow(k, 2);
	double sqrt = std::sqrt(
		-std::pow(oldPos[1], 2)
		+ 2 * oldPos[1] * m_center[1]
		+ 2 * oldPos[1] * k * oldPos[0]
		- 2 * oldPos[1] * k * m_center[0]
		- std::pow(m_center[1], 2)
		- 2 * m_center[1] * k * oldPos[0]
		+ 2 * m_center[1] * k * m_center[0]
		+ k2 * std::pow(m_radius, 2)
		- k2 * std::pow(oldPos[0], 2)
		+ 2 * k2 * oldPos[0] * m_center[0]
		- k2 * std::pow(m_center[0], 2)
		+ std::pow(m_radius, 2)
	);
	double b
	double denominator = (1 / (k2 + 1)) *
	*/

	double k = (relPos[1] / relPos[0]);
	double a = 1 + std::pow(k, 2);
	double b = 2 *
		(
			k * (oldPos[1] - m_center[1] - k * oldPos[0])
			- m_center[0]
			);
	double c =
		oldPos[1] * (oldPos[1] - 2 * (k * oldPos[0] + m_center[1]))
		+ oldPos[0] * k * (oldPos[0] * k + 2 * m_center[1])
		+ std::pow(m_center[0], 2) + std::pow(m_center[1], 2) - std::pow(m_radius, 2);

	return Eigen::Vector2d();
}
