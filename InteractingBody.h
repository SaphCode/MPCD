#pragma once

#ifndef IINTERACTINGBODY_H
#define IINTERACTINGBODY_H

#include <Eigen/Dense>
#include "Body.h"

enum class BodyType {
	PARTICLE,
	MONOMER,
	OBSTACLE,
	WALL,
	CONST_FORCE
};

class InteractingBody :
	public Body
{
public:
	InteractingBody(double mass, Eigen::Vector2d pos, Eigen::Vector2d vel, BodyType type) :
		Body(mass, pos, vel), m_effect(0, 0),
		m_type(type)
	{

	}

	BodyType getType() const {
		return m_type;
	}

	virtual void interact(InteractingBody& b) {}
	virtual void resetEffect() {
		Eigen::Vector2d zero(0, 0);
		m_effect = zero;
	}
	virtual void move(const double timelapse) override;
	void addEffect(Eigen::Vector2d effect);
	void updateVelocity(const double timelapse);

protected:
	Eigen::Vector2d m_effect;

private:
	BodyType m_type;
};
#endif // !IINTERACTINGBODY_H
