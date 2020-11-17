#include "Wall.h"
#include "Constants.h"
#include <iostream>

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
	const Eigen::Vector2d oldPos = o.getOldPosition();
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

bool MPCD::Wall::contains(Eigen::Vector2d point) const
{
	if (m_isUpOOB) {
		if (point[1] >= m_yPos) {
			return true;
		}
	}
	else if (!m_isUpOOB) {
		if (point[1] <= m_yPos) {
			return true;
		}
	}
	return false;
}

bool MPCD::Wall::occupies(std::pair<int, int> index, double cell_dim) const
{
	double y0 = index.first * cell_dim;
	double y1 = ((double)index.first + 1) * cell_dim;
	Vector2d point0(0, y0);
	Vector2d point1(0, y1);
	if (contains(point0) || contains(point1)) {
		//std::cout << "Wall contains this cell." << std::endl;
		return true;
	}
	return false;
}

