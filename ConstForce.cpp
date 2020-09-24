#include "ConstForce.h"
#include <limits>

using namespace Eigen;

ConstForce::ConstForce(double accelerationConstant) :
	InteractingBody(std::numeric_limits<double>::infinity(), Vector2d(0, 0), Vector2d(0, 0))
{
	m_accelerationConstant = accelerationConstant;
}

Eigen::Vector2d ConstForce::interact(InteractingBody& b)
{
	Vector2d effect(m_accelerationConstant, 0);
	b.addEffect(effect);
	return effect;
}
