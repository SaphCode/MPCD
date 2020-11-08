#include "InteractingBody.h"

using namespace Eigen;

void InteractingBody::move(const double timelapse)
{
	Body::move(timelapse);
	m_pos += m_effect * timelapse * timelapse / 2;
}

void InteractingBody::updateVelocity(const double timelapse) {
	m_vel += timelapse * m_effect;
}

void InteractingBody::addEffect(Eigen::Vector2d effect)
{
	m_effect += effect;
}
