#include "Pipe.h"
#include "Grid.h"
#include <execution>
#include <fstream>
#include <iostream>

MPCD::Pipe::Pipe() {
}

void MPCD::Pipe::setParticles(std::vector<Particle> particles)
{
	_particles = particles;
}

void MPCD::Pipe::setObstacles(std::vector<Obstacle> obstacles)
{
	_obstacles = obstacles;
}

void MPCD::Pipe::stream(double lapse, bool draw, std::ofstream& file) {
	for (auto& p : _particles) {
		//for (auto& p : _particles) {
		p.stream(lapse);
		collide(p);
		/*bool inBounds = false;
		while (!inBounds) {
			inBounds = fixOutOfBounds(p);
		}*/
		fixOutOfBounds(p);
		Eigen::Vector2d particlePos = p.getPosition();
		assert(inBounds(particlePos));

		if (draw) {
			Eigen::Vector2d pos = p.getPosition();
			Eigen::Vector2d vel = p.getVelocity();
			file << pos[0] << "," << pos[1] << "," << vel[0] << "," << vel[1] << "\n";
		}
	}
}

bool MPCD::Pipe::inBounds(Eigen::Vector2d pos) {
	double x = pos[0];
	double y = pos[1];
	bool xInBounds = (x >= _x_0 && x <= _x_max);
	bool yInBounds = (y >= _y_0 && y <= _y_max);
	if (xInBounds && yInBounds) {
		return true;
	}
	return false;
}

void MPCD::Pipe::fixOutOfBounds(Particle& p) {
	Eigen::Vector2d particlePos = p.getPosition();
	Eigen::Vector2d newPos = particlePos;
	double diff_posx = particlePos[0] - _x_max;
	double diff_negx = particlePos[0] - _x_0;
	double diff_posy = particlePos[1] - _y_max;
	double diff_negy = particlePos[1] - _y_0;
	double width = _x_max - _x_0;
	double height = _y_max - _y_0;

	if (diff_posx > 0) {
		double rem = std::fmod(diff_posx, width);
		assert(rem > 0);
		newPos[0] = rem;
	}
	else if (diff_negx < 0) {
		double rem = std::fmod(diff_negx, width);
		assert(rem < 0);
		newPos[0] = _x_max + rem;
	}
	if (diff_posy > 0) {
		double rem = std::fmod(diff_posy, height);
		assert(rem > 0);
		newPos[1] = rem;
	}
	else if (diff_negy < 0) {
		double rem = std::fmod(diff_negy, height);
		assert(rem < 0);
		newPos[1] = _y_max + rem;
	}
	p.correctPosition(newPos);
	/*if ((newPos[0] >= _x_0) && (newPos[0] <= _x_max) && (newPos[1] >= _y_0) && (newPos[1] <= _y_max)) {
		return true;
	}
	else return false;*/
}

void MPCD::Pipe::collide(Particle& p) {

}

std::vector<MPCD::Particle>& MPCD::Pipe::getParticles()
{
	return _particles;
}
