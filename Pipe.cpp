#include "Pipe.h"
#include "Grid.h"
#include <execution>
#include <fstream>
#include <iostream>

MPCD::Pipe::Pipe(ConstForce force) :
	m_constForce(force)
{
}

void MPCD::Pipe::setObstacles(std::vector<CircularObstacle> obstacles, std::vector<Wall> walls)
{
	m_obstacles = obstacles;
	m_walls = walls;
}

void MPCD::Pipe::stream(std::vector<Particle>& particles, double lapse, bool draw, std::ofstream& file) {
	int size = particles.size();
	#pragma omp parallel for
	for (int i = 0; i < size; i++) {
		Particle& p = particles[i];
		for (auto& wall : m_walls) {
			wall.interact(p);
		}
		for (auto& obstacle : m_obstacles) {
			obstacle.interact(p);
		}
		m_constForce.interact(p);

		p.move(lapse);
		collide(p);
		fixOutOfBounds(p);
		p.updateVelocity(lapse);
		p.resetEffect();
		if (draw) {
			Eigen::Vector2d pos = p.getPosition();
			Eigen::Vector2d vel = p.getVelocity();
			#pragma omp critical
			{
				file << pos[0] << "," << pos[1] << "," << vel[0] << "," << vel[1] << "\n";
			}
		}
		
	}
}

bool MPCD::Pipe::inBounds(const Eigen::Vector2d& pos) {
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

	p.correctPosition(newPos);
	if (!inBounds(newPos)) {
		std::cout << "Old pos: (" << particlePos[0] << ", " << particlePos[1] << ")\n";
		std::cout << "New pos: (" << newPos[0] << ", " << newPos[1] << ")\n";
	}
	assert(inBounds(newPos));

}

void MPCD::Pipe::collide(Particle& p) {
	for (auto &o : m_obstacles) {
		if (o.isInBounds(p)) {
			Eigen::Vector2d overshoot = o.getOvershoot(p);
			// no slip
			p.collided(overshoot);
			assert(!(o.isInBounds(p)));
		}
	}
	for (auto& o : m_walls) {
		if (o.isInBounds(p)) {
			Eigen::Vector2d overshoot = o.getOvershoot(p);
			// no slip
			p.collided(overshoot);
			assert(!(o.isInBounds(p)));
		}
	}
	for (const auto& o : m_obstacles) {
		assert(!(o.isInBounds(p))); // double check for now, TODO maybe remove
	}
	for (const auto& o : m_walls) {
		assert(!(o.isInBounds(p))); // double check for now, TODO maybe remove
	}
}

