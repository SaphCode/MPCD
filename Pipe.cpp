#include "Pipe.h"
#include "Grid.h"
#include <execution>
#include <iostream>
#include <filesystem>
#include <fstream>

MPCD::Pipe::Pipe(ConstForce force) :
	m_constForce(force)
{
	_w = 5; // width of 0s for filenames
}

void MPCD::Pipe::setObstacles(std::vector<CircularObstacle> obstacles, std::vector<Wall> walls)
{
	m_obstacles = obstacles;
	m_walls = walls;
}

void MPCD::Pipe::verlet(std::vector<Monomer>& monomers, bool draw, int t)
{
	for (int md_t = 0; md_t < MPCD::Constants::num_md_timesteps; md_t++) {
		verletPosition(monomers);
		verletVelocity(monomers);
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

		for (const auto& m : monomers) {
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
	const Eigen::Vector2d constForce(MPCD::Constants::const_force, 0);
	
	#pragma omp parallel for
	for (int i = 0; i < particles.size(); i++) {
		Particle& p = particles[i];

		/* NOT NECESSARY: NO INTERACTION WITH WALL AND OBSTACLES; and const force is simple update.
		{
			for (auto& wall : m_walls) {
				wall.interact(p);
			}
			for (auto& obstacle : m_obstacles) {
				obstacle.interact(p);
			}
			m_constForce.interact(p);
		}
		*/
		
		p.constForce(constForce);

		p.move(lapse);

		collide(p);
		fixOutOfBounds(p);
		p.updateVelocity(lapse);
		p.resetEffect();
		if (draw) {
			const Eigen::Vector2d pos = p.getPosition();
			const Eigen::Vector2d vel = p.getVelocity();
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
	for (const auto &o : m_obstacles) {
		if (o.isInBounds(b)) {
			Eigen::Vector2d overshoot = o.getOvershoot(b);
			// no slip
			b.collided(overshoot);
			assert(!(o.isInBounds(b)));
		}
	}
	for (const auto& o : m_walls) {
		if (o.isInBounds(b)) {
			
			Eigen::Vector2d overshoot = o.getOvershoot(b);
			if (std::abs(overshoot[0]) > 1 || std::abs(overshoot[1]) > 1) {
				//std::cout << overshoot << std::endl;
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

void MPCD::Pipe::verletPosition(std::vector<Monomer>& monomers)
{
	#pragma omp parallel for
	for (int i = 0; i < monomers.size(); i++) {
		Monomer& m = monomers[i];

		m.move(MPCD::Constants::md_timestep); // update position
		collide(m);
		fixOutOfBounds(m);

	}
}

void MPCD::Pipe::verletVelocity(std::vector<Monomer>& monomers)
{
	#pragma omp parallel for
	for (int i = 0; i < monomers.size(); i++) {
		Monomer& m = monomers[i];

		Eigen::Vector2d oldEffect = m.getEffect();
		m.resetEffect();

		calculateInteraction(i, m, monomers);

		m.updateVelocity(MPCD::Constants::md_timestep, oldEffect); // update velocity
	}
}

void MPCD::Pipe::calculateInteraction(int currentIndex, Monomer& m, std::vector<Monomer>& monomers) {

	
	for (int m_i = 0; m_i < monomers.size(); m_i++) {
		if (currentIndex != m_i) {
			const Monomer& m2 = monomers[m_i];
			const Eigen::Vector2d rel = m.getRelPositionTorus(m2.getPosition());
			
			if (currentIndex == m_i - 1 || currentIndex == m_i + 1) {
				m.linearSpring(rel);
			}

			m.monomerInteraction(rel, MPCD::Constants::monomerMonomer_interaction_tuning, m.getDiameter()); // 1/2 b/c of double counting
		}

	}

	for (int k = 0; k < m_obstacles.size(); k++) {
		CircularObstacle& co = m_obstacles[k];
		m.interact(co);
	}

	for (int l = 0; l < m_walls.size(); l++) {
		Wall& wall = m_walls[l];
		m.interact(wall);
	}

	// m_constForce.interact(m);
}


