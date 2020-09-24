#pragma once

#ifndef PARTICLE_H
#define PARTICLE_H
#include <Eigen/Dense>
#include "InteractingBody.h"

namespace MPCD {

	class Particle : 
		public InteractingBody
	{
	public:
		/* Creates a particle with @param position, @param velocity and unit mass. */
		Particle(double mass, Eigen::Vector2d position, Eigen::Vector2d velocity) :
			InteractingBody(mass, position, velocity) {

		}
		
		/* Updates the velocity of the particle using the MPCD cell collision algorithm.
		@param mean_cell_velocity: the mean velocity of the cell the particle belongs to. NOTE:
				This means that the galilean invariance step, where the particle is shifted, needs to be
				applied first
		@param rotationMatrix: the rotation Matrix for the cell the particle belongs to. NOTE: same as above
		*/
		void collide(Eigen::Vector2d mean_cell_velocity, double rotationAngle);

		void move(double timelapse) override;

		Eigen::Vector2d interact(InteractingBody& b) override {
			// no effect on other bodies
			return Eigen::Vector2d(0, 0);
		}

		//friend bool operator==(const Particle& lhs, const Particle& rhs);
	};

}
#endif
