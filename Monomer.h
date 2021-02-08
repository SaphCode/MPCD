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

	void monomerInteraction(Eigen::Vector2d rel, double tuning, double diameter);
	void nonlinearSpring(Eigen::Vector2d rel, double tuning, double diameter);

	Eigen::Vector2d getRelPositionTorus(Eigen::Vector2d otherPos);

private:
	std::pair<int, int> m_coords;
	double m_diameter;

	Eigen::Vector2d truncshiftedLennardJones(Eigen::Vector2d rel, double tuning, double diameter);
	Eigen::Vector2d truncLennardJonesWall(Eigen::Vector2d rel, double tuning, double diameter);

	Eigen::Vector2d lennardJones(Eigen::Vector2d rel, double tuning, double diameter);
};

