#include "Wall.h"
#include "Constants.h"
#include <iostream>

using namespace MPCD;
using namespace Eigen;

bool MPCD::Wall::isInBounds(const Body& o) const
{
	const Eigen::Vector2d pos = o.getTorusPosition();
	return contains(pos);
}

Eigen::Vector2d MPCD::Wall::getOvershoot(const Body& o) const
{
	const Eigen::Vector2d pos = o.getTorusPosition();
	const Eigen::Vector2d oldPos = o.getOldTorusPosition();
	
	Vector2d rel = pos - oldPos;

	double k = rel[1] / rel[0];

	double yDiff = pos[1] - m_yPos;
	double xDiff = yDiff * 1.0 / k;
	Eigen::Vector2d overshoot(xDiff, yDiff);

	return overshoot;
}

void MPCD::Wall::interact(InteractingBody& b)
{
	// do nothing for now
}

bool MPCD::Wall::contains(Eigen::Vector2d point) const
{
	bool upOOB = m_yPos >= MPCD::Constants::y_max;
	if (upOOB && point[1] >= m_yPos) {
		return true;
	}
	else if (!upOOB && point[1] <= m_yPos) {
		return true;
	}
	return false;
}

bool MPCD::Wall::occupies(std::pair<int, int> index, Eigen::Vector2d shift, double cell_dim) const
{
	double y0 = index.first * cell_dim + shift[1];
	double y1 = ((double)index.first + 1) * cell_dim + shift[1];
	Vector2d point0(0, y0);
	Vector2d point1(0, y1);
	if (contains(point0) || contains(point1)) {
		return true;
	}
	return false;
}

