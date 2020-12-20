#pragma once
#include "InteractingBody.h"

class Monomer :
    public InteractingBody
{
public:

	Monomer(double mass, Eigen::Vector2d pos, Eigen::Vector2d vel, double diameter) :
		InteractingBody(mass, pos, vel, BodyType::MONOMER),
		m_diameter(diameter)
	{}

	void interact(InteractingBody& b) override;

	virtual void move(const double timelapse) override;

	void updateVelocity(const double timelapse, Eigen::Vector2d oldEffect) {
		m_vel += 1 / 2 * timelapse * m_effect / m_mass + 1 / 2 * timelapse * oldEffect / m_mass;
	}

	void collide(Eigen::Vector2d mean_cell_velocity, double rotationAngle, double temperatureScalingFactor);

	std::pair<int, int> getCoordinates() const {
		return m_coords;
	}

	void setCoordinates(std::pair<int, int> coords) {
		m_coords = coords;
	}

	Eigen::Vector2d getEffect() const {
		return m_effect;
	}

	double getDiameter() const {
		return m_diameter;
	}

private:
	std::pair<int, int> m_coords;
	double m_diameter;

	Eigen::Vector2d truncLennardJones(Eigen::Vector2d rel, double tuning, double diameter);
};
