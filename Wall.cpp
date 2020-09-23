#include "Wall.h"
#include "Constants.h"

using namespace MPCD;

Wall::Wall(const Rectangle rect, const ForceType type) : ImmovableObstacle(Eigen::Vector2d(rect.x0, rect.y0), type), _bounds(rect) {

}

bool MPCD::Wall::isInBounds(const PhysicalObject& o) const
{
	const Eigen::Vector2d pos = o.getPosition();
	if (((pos[0] >= _bounds.x0) && pos[0] <= _bounds.xMax) && ((pos[1] >= _bounds.y0) && (pos[1] <= _bounds.yMax))) {
		return true;
	}
	return false;
}

Eigen::Vector2d MPCD::Wall::getCollisionPoint(const PhysicalObject& o) const
{
	const Eigen::Vector2d pos = o.getPosition();
	const Eigen::Vector2d oldPos = o.getOldPosition(MPCD::Constants::time_lapse);
	Eigen::Vector2d collisionPoint; // TODO
}

