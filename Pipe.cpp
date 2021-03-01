#include "Pipe.h"
#include "Grid.h"
#include <execution>
#include <iostream>
#include <filesystem>
#include <fstream>

MPCD::Pipe::Pipe(ConstForce force) :
	m_constForce(force)
{
	_w = 6; // width of 0s for filenames
}

void MPCD::Pipe::setObstacles(std::vector<CircularObstacle> obstacles, std::vector<Wall> walls)
{
	m_obstacles = obstacles;
	m_walls = walls;
}

const std::vector<MPCD::CircularObstacle>& MPCD::Pipe::getObstacles() const
{
	return m_obstacles;
}

const std::vector<MPCD::Wall>& MPCD::Pipe::getWalls() const
{
	return m_walls;
}

void MPCD::Pipe::verlet(std::vector<Monomer>& monomers, bool draw, int t)
{
	if (draw) {
		std::ofstream outFile;
		std::stringstream s;
		std::stringstream av;
		std::string filename;
		s << std::setfill('0') << std::setw(_w) << t;
		av << "av" << Constants::average_particles_per_cell << "_";

		std::string external("G:/Bachelor/Data/");
		filename = external + "monomers_" + av.str() + "timestep" + s.str() + ".csv";//cwd.string()
		outFile = std::ofstream(filename);
		outFile << "x,y,vx,vy,m,sigma\n";

		for (const auto& m : monomers) {
			Eigen::Vector2d pos = m.getTorusPosition();
			Eigen::Vector2d vel = m.getVelocity();
			outFile << pos[0] << "," << pos[1] << "," << vel[0] << "," << vel[1] << "," << m.getMass() << "," << m.getDiameter() << "\n";
		}
	}

	for (int md_t = 0; md_t < MPCD::Constants::num_md_timesteps; md_t++) {
		#pragma omp parallel for
		for (int i = 0; i < monomers.size(); i++) {
			verletPosition(i, monomers, MPCD::Constants::md_timestep);
		}
		#pragma omp parallel for
		for (int i = 0; i < monomers.size(); i++) {
			verletVelocity(i, monomers, MPCD::Constants::md_timestep);
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
		std::string external("G:/Bachelor/Data/");
		filename = external + "particles_" + av.str() + "timestep" + s.str() + ".csv";
		outFile = std::ofstream(filename);
		outFile << "x,y,vx,vy\n";
	}
	const Eigen::Vector2d constForce(MPCD::Constants::const_force, 0);
	
	#pragma omp parallel for
	for (int i = 0; i < particles.size(); i++) {
		Particle& p = particles[i];

		p.constForce(constForce);

		p.move(lapse);

		collide(p);
		fixOutOfBounds(p);
		p.updateVelocity(lapse);
		p.resetEffect();
		if (draw) {
			const Eigen::Vector2d pos = p.getTorusPosition();
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
	bool xInBounds = (x >= MPCD::Constants::x_0); // && x <= _x_max
	bool yInBounds = (y >= MPCD::Constants::y_0 && y <= MPCD::Constants::y_max);
	if (xInBounds && yInBounds) {
		return true;
	}
	return false;
}

void MPCD::Pipe::fixOutOfBounds(Body& p) {
	Eigen::Vector2d pos = p.getTorusPosition();

	Eigen::Vector2d newPos = pos;
	double diff_negx = pos[0] - MPCD::Constants::x_0;
	double width = MPCD::Constants::x_max - MPCD::Constants::x_0;
	double height = MPCD::Constants::y_max - MPCD::Constants::y_0;

	/*
	if (diff_posx > 0) {
		double rem = std::fmod(diff_posx, width);
		assert(rem > 0);
		newPos[0] = rem;
	}
	*/
	if (diff_negx < 0) {
		double rem = std::fmod(diff_negx, width);
		assert(rem < 0);
		newPos[0] = MPCD::Constants::x_max + rem;
		p.correctPosition(newPos);
	}

	if (!inBounds(newPos)) {
		//std::cout << "Vel: " << p.getVelocity() << "\n";
		std::cout << "Mass: " << p.getMass() << "\n";
		std::cout << "Vel: " << p.getVelocity() << "\n";
		std::cout << "Old pos: (" << pos[0] << ", " << pos[1] << ")\n";
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
			std::cout << "old pos, new pos, vel" << std::endl; // TODO
			std::cout << b.getOldTorusPosition() << std::endl;
			std::cout << b.getTorusPosition() << std::endl;
			std::cout << b.getVelocity() << std::endl;
		}
		assert(!(o.isInBounds(b))); // double check for now, TODO maybe remove
	}
}

void MPCD::Pipe::verletPosition(int chainIndex, std::vector<Monomer>& monomers, double timestep)
{
	Monomer& m = monomers[chainIndex];
	Eigen::Vector2d oldPos = m.getPosition();
	m.move(timestep); // update position
	//collide(m);
	/*bool isLegalPos = isLegalPosition(m);

	if (isLegalPos) return;
	else {
		
		/*
		std::cout << "Vel: " << m.getVelocity() << "\n";
		std::cout << "F: " << m.getEffect() << "\n";
		std::cout << "Pos: " << m.getPosition() << "\n";
		double count = 0;
		double partition = 0;
		while (!isLegalPos) {
			count++;
			partition = std::pow(2, count);
			m.correctPosition(oldPos);
			m.move(timestep / partition);
			isLegalPos = isLegalPosition(m);
		}
		for (int t = 1; t < count; t++) {
			std::cout << "Correction step " << t + 1 << " of " << count << ".\n";
			verletVelocity(chainIndex, monomers, timestep / partition);
			verletPosition(chainIndex, monomers, timestep / partition);
		}
		std::cout << "Final Vel: " << m.getVelocity() << "\n";
		std::cout << "Final F: " << m.getEffect() << "\n";
		std::cout << "Final Pos: " << m.getPosition() << "\n";
		
	}
	/*
	#pragma omp critical 
	{
		unsigned int count = 0;
		while (!isLegalPos) {
			count += 1;
			timestep /= (double)std::pow(2, count);

			m.correctPosition(oldPos);
			m.move(timestep);
			oldPos = m.getPosition();
			isLegalPos = isLegalPosition(m);
		}

		unsigned int partition = std::ceil(std::pow(2, count));

		for (int p = 0; p < partition; p++) {
			verletVelocity(chainIndex, monomers, timestep);
			verletPosition(chainIndex, monomers, timestep);
		}
	}
	*/
}

const bool MPCD::Pipe::isLegalPosition(const Body& b)
{
	for (const auto& o : m_obstacles) {
		if (o.isInBounds(b)) {
			return false;
		}
	}
	for (const auto& o : m_walls) {
		if (o.isInBounds(b)) {
			return false;
		}
	}
	return true;
}

void MPCD::Pipe::verletVelocity(int chainIndex, std::vector<Monomer>& monomers, double timestep)
{
	Monomer& m = monomers[chainIndex];
	Eigen::Vector2d oldEffect = m.getEffect();
	m.resetEffect();
	calculateInteraction(chainIndex, monomers);
	m.updateVelocity(timestep, oldEffect); // update velocity
}

void MPCD::Pipe::calculateInteraction(int chainIndex, std::vector<Monomer>& monomers) {

	Monomer& m = monomers[chainIndex];
	for (int m_i = 0; m_i < monomers.size(); m_i++) {
		if (chainIndex != m_i) {
			const Monomer& m2 = monomers[m_i];
			const Eigen::Vector2d rel = m.getRelPositionTorus(m2.getPosition());
			
			if ((chainIndex == m_i - 1) || (chainIndex == m_i + 1)) { // right & left neighbor
				m.nonlinearSpring(rel, MPCD::Constants::monomerMonomer_interaction_tuning, m.getDiameter());
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
}


