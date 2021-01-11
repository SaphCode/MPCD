#include "ConstForce.h"
#include <limits>

using namespace Eigen;

ConstForce::ConstForce(Eigen::Vector2d force) :
	InteractingBody(std::numeric_limits<double>::infinity(), Vector2d(0, 0), Vector2d(0, 0), BodyType::CONST_FORCE)
{
	m_force = force;
}

void ConstForce::interact(InteractingBody& b)
{
	Vector2d effect = m_force;
	b.addEffect(effect);
}
