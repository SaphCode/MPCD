#include "Wall.h"
#include "Constants.h"

using namespace MPCD;
using namespace Eigen;

Wall::Wall(const double yPos, const ForceType type) : ImmovableObstacle(Eigen::Vector2d(0, yPos), type), _yPos(yPos) {

}

bool MPCD::Wall::isInBounds(const PhysicalObject& o) const
{
	const Eigen::Vector2d pos = o.getPosition();
	if ((pos[1] > MPCD::Constants::y_0) && (pos[1] <= MPCD::Constants::y_max)) {
		return false;
	}
	return true;
}

Eigen::Vector2d MPCD::Wall::getOvershoot(const PhysicalObject& o) const
{
	const Eigen::Vector2d pos = o.getPosition();
	const Eigen::Vector2d oldPos = o.getOldPosition(MPCD::Constants::time_lapse);
	Vector2d rel = pos - oldPos;

	double k = rel[1] / rel[0];

	double yDiff = pos[1] - _yPos;
	double xDiff = yDiff * 1 / k;
	Eigen::Vector2d overshoot(xDiff, yDiff);

	return overshoot;
}

