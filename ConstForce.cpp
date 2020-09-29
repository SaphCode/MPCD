#include "ConstForce.h"
#include <limits>

using namespace Eigen;

ConstForce::ConstForce(Eigen::Vector2d acceleration) :
	InteractingBody(std::numeric_limits<double>::infinity(), Vector2d(0, 0), Vector2d(0, 0))
{
	m_acceleration = acceleration;
}

Eigen::Vector2d ConstForce::interact(InteractingBody& b)
{
	Vector2d effect = m_acceleration;
	b.addEffect(effect);
	return effect;
}
