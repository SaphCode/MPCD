#include "CircularObstacle.h"

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
	return Eigen::Vector2d();
}
