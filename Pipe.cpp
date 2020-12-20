#include "Pipe.h"
#include "Grid.h"
#include <execution>
#include <iostream>
#include <filesystem>
#include <fstream>

MPCD::Pipe::Pipe(ConstForce force) :
	m_constForce(force)
{
	const int timesteps = Constants::timesteps;
	if (timesteps < 99) {
		_w = 2;
	}
	else if (timesteps < 999) {
		_w = 3;
	}
	else if (timesteps < 9999) {
		_w = 4;
	}
	else if (timesteps < 99999) {
		_w = 5;
	}
	else {
		throw std::exception("timesteps too large or too short: " + timesteps);
	}
}

void MPCD::Pipe::setObstacles(std::vector<CircularObstacle> obstacles, std::vector<Wall> walls)
{
	m_obstacles = obstacles;
	m_walls = walls;
}

void MPCD::Pipe::verlet(std::vector<Particle>& particles, std::vector<Monomer>& monomers, bool draw, int t)
{
	for (int md_t = 0; md_t < MPCD::Constants::num_md_timesteps; md_t++) {
		verletPosition(particles, monomers);
		verletVelocity(particles, monomers);
	}


	if (draw) {
		std::ofstream outFile;
		std::stringstream s;
		std::stringstream av;
		std::string filename;
		s << std::setfill('0') << std::setw(_w) << t;
		av << "av" << Constants::average_particles_per_cell << "_";
		filename = "../../Analysis/" + std::string("Data/") + "monomers_" + av.str() + "timestep" + s.str() + ".csv";//cwd.string()
		outFile = std::ofstream(filename);
		outFile << "x,y,vx,vy,m,sigma\n";

		for (auto& m : monomers) {
			Eigen::Vector2d pos = m.getPosition();
			Eigen::Vector2d vel = m.getVelocity();
			outFile << pos[0] << "," << pos[1] << "," << vel[0] << "," << vel[1] << "," << m.getMass() << "," << m.getDiameter() << "\n";
		}
	}
	
}

void MPCD::Pipe::stream(std::vector<Particle>& particles, double lapse, bool draw, int t) {
	std::ofstream outFile;

	if (draw) {
		std::stringstream s;
		std::stringstream av;
		std::string filename;
		s << std::setfill('0') << std::setw(_w) << t;
		av << "av" << Constants::average_particles_per_cell << "_";
		filename = "../../Analysis/" + std::string("Data/") + "particles_" + av.str() + "timestep" + s.str() + ".csv";//cwd.string()
		outFile = std::ofstream(filename);
		outFile << "x,y,vx,vy\n";
	}
	#pragma omp parallel for
	for (int i = 0; i < particles.size(); i++) {
		Particle& p = particles[i];

		#pragma omp critical 
		{
			for (auto& wall : m_walls) {
				wall.interact(p);
			}
			for (auto& obstacle : m_obstacles) {
				obstacle.interact(p);
			}
			m_constForce.interact(p);
		}
		

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
				outFile << pos[0] << "," << pos[1] << "," << vel[0] << "," << vel[1] << "\n";
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

void MPCD::Pipe::fixOutOfBounds(Body& p) {
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

void MPCD::Pipe::collide(Body& b) {
	for (auto &o : m_obstacles) {
		if (o.isInBounds(b)) {
			Eigen::Vector2d overshoot = o.getOvershoot(b);
			// no slip
			b.collided(overshoot);
			assert(!(o.isInBounds(b)));
		}
	}
	for (auto& o : m_walls) {
		if (o.isInBounds(b)) {
			
			Eigen::Vector2d overshoot = o.getOvershoot(b);
			if (std::abs(overshoot[0]) > 1 || std::abs(overshoot[1]) > 1) {
				std::cout << overshoot << std::endl;
				o.getOvershoot(b);
			}
			// no slip
			b.collided(overshoot);
			assert(!(o.isInBounds(b)));
		}
	}
	for (const auto& o : m_obstacles) {
		assert(!(o.isInBounds(b))); // double check for now, TODO maybe remove
	}
	for (const auto& o : m_walls) {
		if ((o.isInBounds(b))) {
			std::cout << "old pos, new pos, vel" << std::endl;
			std::cout << b.getOldPosition() << std::endl;
			std::cout << b.getPosition() << std::endl;
			std::cout << b.getVelocity() << std::endl;
		}
		assert(!(o.isInBounds(b))); // double check for now, TODO maybe remove
	}
}

void MPCD::Pipe::verletPosition(std::vector<Particle>& particles, std::vector<Monomer>& monomers)
{

	for (int i = 0; i < monomers.size(); i++) {
		Monomer& m = monomers[i];

		m.move(MPCD::Constants::md_timestep); // update position
		collide(m);
		fixOutOfBounds(m);

		#pragma omp parallel for
		for (int j = 0; j < particles.size(); j++) {
			Particle& p = particles[j];

			p.move(MPCD::Constants::md_timestep);
			collide(p);
			fixOutOfBounds(p);
			p.updateVelocity(MPCD::Constants::md_timestep);
			p.resetEffect();
		}

		/*
		std::pair<int, int> coords = m.getCoordinates();

		int box = 10;

		int i_s = coords.first - box < 0 ? _numRows - 1 - (box - coords.first) : coords.first - box;
		int i_t = coords.first + box > _numRows - 1 ? (coords.first + box - (_numRows - 1)) : coords.first + box;
		int j_s = coords.second - box < 0 ? _numCols - 1 - (box - coords.second) : coords.second - box;
		int j_t = coords.second + box > _numCols - 1 ? (coords.second + box - (_numCols - 1)) : coords.second + box;

		int s;
		if (i_s > i_t) {
			for (int i = i_s; i < _numRows; i++) {

			}
		}
		*/
	}
}

void MPCD::Pipe::verletVelocity(std::vector<Particle>& particles, std::vector<Monomer>& monomers)
{
	for (int i = 0; i < monomers.size(); i++) {
		Monomer& m = monomers[i];

		Eigen::Vector2d oldEffect = m.getEffect();
		m.resetEffect();

		calculateInteraction(i, m, particles, monomers);

		m.updateVelocity(MPCD::Constants::md_timestep, oldEffect); // update velocity
	}
}

void MPCD::Pipe::calculateInteraction(int currentIndex, Monomer& m, std::vector<Particle>& particles, std::vector<Monomer>& monomers) {
	#pragma omp parallel for
	for (int j = 0; j < particles.size(); j++) {
		Particle& p = particles[j];

		m.interact(p);
	}

	#pragma omp parallel for
	for (int m_i = 0; m_i < monomers.size(); m_i++) {
		if (currentIndex != m_i) {
			Monomer& m2 = monomers[m_i];

			m.interact(m2);
		}

	}

	#pragma omp parallel for
	for (int k = 0; k < m_obstacles.size(); k++) {
		CircularObstacle& co = m_obstacles[k];
		m.interact(co);
	}

	for (int l = 0; l < m_walls.size(); l++) {
		Wall& wall = m_walls[l];
		m.interact(wall);
	}

	m_constForce.interact(m);
}


