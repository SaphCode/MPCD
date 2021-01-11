#pragma once

#include <Eigen/Dense>
#include "InteractingBody.h"

namespace MPCD {

	class Particle : 
		public InteractingBody
	{
	public:
		/* Creates a particle with @param position, @param velocity and unit mass. */
		Particle(double mass, Eigen::Vector2d position, Eigen::Vector2d velocity) :
			InteractingBody(mass, position, velocity, BodyType::PARTICLE)
		{

		}
		
		/* Updates the velocity of the particle using the MPCD cell collision algorithm.
		@param mean_cell_velocity: the mean velocity of the cell the particle belongs to.
				This means that the galilean invariance step, where the particle is shifted, needs to be
				applied first
		@param rotationMatrix: the rotation Matrix for the cell the particle belongs to.
		@param scalingFactor: scaling from cell-level thermostat
		*/
		void collide(Eigen::Vector2d mean_cell_velocity, double rotationAngle, double temperatureScalingFactor);

		void move(double timelapse) override;

		std::pair<int, int> getCoordinates() const {
			return m_coordinates;
		}
		void setCoordinates(std::pair<int, int> coord) {
			m_coordinates = coord;
		}

		void interact(InteractingBody& b) override {
			// no effect on other bodies
		}

		void constForce(Eigen::Vector2d force);
	private:
		std::pair<int, int> m_coordinates;

	};

}
