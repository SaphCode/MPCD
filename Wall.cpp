#include "Wall.h"
#include "Constants.h"

using namespace MPCD;
using namespace Eigen;

bool MPCD::Wall::isInBounds(const Body& o) const
{
	const Eigen::Vector2d pos = o.getPosition();
	if (m_isUpOOB) {
		if (pos[1] >= m_yPos) {
			return true;
		}
	}
	else {
		if (pos[1] <= m_yPos) {
			return true;
		}
	}
	return false;
}

Eigen::Vector2d MPCD::Wall::getOvershoot(const Body& o) const
{
	const Eigen::Vector2d pos = o.getPosition();
	const Eigen::Vector2d oldPos = o.getOldPosition(MPCD::Constants::time_lapse);
	Vector2d rel = pos - oldPos;

	double k = rel[1] / rel[0];

	double yDiff = pos[1] - m_yPos;
	double xDiff = yDiff * 1 / k;
	Eigen::Vector2d overshoot(xDiff, yDiff);

	return overshoot;
}

Eigen::Vector2d MPCD::Wall::interact(InteractingBody& b)
{
	// do nothing for now
	return Vector2d(0, 0);
}

