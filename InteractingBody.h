#pragma once

#ifndef IINTERACTINGBODY_H
#define IINTERACTINGBODY_H

#include <Eigen/Dense>
#include "Body.h"

class InteractingBody :
	public Body
{
public:
	InteractingBody(double mass, Eigen::Vector2d pos, Eigen::Vector2d vel) :
		Body(mass, pos, vel), m_effect(0, 0)
	{

	}

	virtual Eigen::Vector2d interact(InteractingBody& b) {

		return Eigen::Vector2d(0, 0); // default implementation
	}
	virtual void resetEffect() {
		Eigen::Vector2d zero(0, 0);
		m_effect = zero;
	}
	virtual void move(const double timelapse) override;
	void addEffect(Eigen::Vector2d effect);
	Eigen::Vector2d getOldPosition(const double timelapse) const override;
	void updateVelocity(const double timelapse);

private:
	Eigen::Vector2d m_effect;
};
#endif // !IINTERACTINGBODY_H
