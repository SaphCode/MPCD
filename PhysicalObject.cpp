#include "PhysicalObject.h"
#include "Force.h"

using namespace Eigen;

Vector2d PhysicalObject::getEffect(const PhysicalObject& o) const
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
        throw std::exception("Not implemented for this force type yet.");
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

void PhysicalObject::updatePosition(double timelapse) {
    Vector2d force = calculateForceOnThis();
    _pos += timelapse * _vel;
    _pos += force / _mass * timelapse * timelapse / 2;
    /*switch (_f) {
    case ForceType::CONST: {
        
    }
    break;
    case ForceType::CONST_X: {
        _pos += force / _mass * timelapse * timelapse / 2;
    }
    break;
    case ForceType::CONST_Y: {
        _pos += force / _mass * timelapse * timelapse / 2;
    }
    break;
    default:
    {
        throw std::exception("Force type not implemented yet.");
    }
    }*/
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
