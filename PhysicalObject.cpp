#include "PhysicalObject.h"
#include "Force.h"

using namespace Eigen;

Vector2d PhysicalObject::getEffect(PhysicalObject o) const
{
    switch (_f) {
    case ForceType::CONST_X:
    {
        Vector2d force(Force::const_force, 0);
        return force;
    }
    break;
    case ForceType::NO_FORCE:
    {
        Vector2d force(0, 0);
        return force;
    }
    break;
    default:
        throw std::exception();
    }
}

Eigen::Vector2d PhysicalObject::getPosition() const
{
    return _pos;
}

Eigen::Vector2d PhysicalObject::getVelocity() const
{
    return _vel;
}

void PhysicalObject::updateVelocity(double timelapse)
{
    Vector2d force = calculateForceOnThis();
    _vel += timelapse * force/_mass;
}

Vector2d PhysicalObject::calculateForceOnThis() const
{
    Vector2d force(0,0);
    for (auto& o : _registeredObjects) {
        force += o->getEffect(*this);
    }
    return force;
}

void PhysicalObject::registerObject(PhysicalObject& o)
{
    _registeredObjects.push_back(&o);
}

void PhysicalObject::unregisterObject(PhysicalObject& o)
{
    _registeredObjects.erase(std::remove(_registeredObjects.begin(), _registeredObjects.end(), &o), _registeredObjects.end());
}
