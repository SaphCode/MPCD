#include "InteractingBody.h"

using namespace Eigen;

void InteractingBody::move(const double timelapse)
{
	Body::move(timelapse);
	m_pos += m_effect * timelapse * timelapse / 2;
	Vector2d zero(0, 0);
	m_effect = zero;
}

void InteractingBody::updateVelocity(const double timelapse) {
	m_vel += timelapse * m_effect;
}

void InteractingBody::addEffect(Eigen::Vector2d effect)
{
	m_effect += effect;
}

// TODO: Here might be a mistake, if the velocity is updated before the old Position is needed.
Eigen::Vector2d InteractingBody::getOldPosition(const double timelapse) const
{
	Eigen::Vector2d oldPos = Body::getOldPosition(timelapse);
	oldPos -= m_effect * timelapse * timelapse / 2;
	return oldPos;
}
