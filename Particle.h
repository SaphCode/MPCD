#pragma once


#ifndef PARTICLE_H
#define PARTICLE_H
#include <Eigen/Dense>
namespace MPCD {

	class Particle
	{
	private:
		Eigen::Vector2d _position;
		Eigen::Vector2d _velocity;
		Eigen::Vector2i _cell_index;
		// Eigen::Vector2i _cell_index; should not need to store this
		// Eigen::Vector2d _shifted_position; should not need to store this

		/* Update the Cell index after movement has occured. Should only be called after shift. (works on the shifted position, not the position)  */
		void updateCellPosition(Eigen::Vector2d shiftedPosition);

		//Vector2d _position;
		//Vector2d _velocity;
		//double mass; // could this be int?, do i need this
	public:
		//Eigen::Vector2d _position;
		//Eigen::Vector2d _velocity;
		/* Creates a particle with @param position, @param velocity and unit mass. */
		Particle(Eigen::Vector2d position, Eigen::Vector2d velocity);
		/* Constructor strictly for std::vector<Particle> */
		Particle(); // default constructor needed for custom class array

		~Particle();

		///* Gives the Cell the Particle is in as a 2d Vector. */
		Eigen::Vector2i getCellIndex();

		/* If you make the variables PRIVATE */

		Eigen::Vector2d getPosition();
		Eigen::Vector2d getVelocity();

		/* Shifts the position of the particle by @param amount and updates the particle's cell position. Keeps the shifted position in a separate variable so 
		 it does not forget where it was.
		 @return the new cell index of the particle*/
		Eigen::Vector2i shift(Eigen::Vector2d amount);
		//void setPosition(Eigen::Vector2d position);
		//void setVelocity(Eigen::Vector2d velocity);
		/* Moves the particle according to its current velocity. */
		void move(double time_step);

		/* Updates the velocity of the particle using the MPCD cell collision algorithm.
		@param mean_cell_velocity: the mean velocity of the cell the particle belongs to. NOTE:
				This means that the galilean invariance step, where the particle is shifted, needs to be
				applied first
		@param rotationMatrix: the rotation Matrix for the cell the particle belongs to. NOTE: same as above
		*/
		void updateVelocity(Eigen::Vector2d mean_cell_velocity, Eigen::Matrix2d rotationMatrix);
		friend std::ostream& operator<<(std::ostream& output, const Particle& p);
	};

	std::ostream& operator<<(std::ostream& output, const Particle& H);
}
#endif