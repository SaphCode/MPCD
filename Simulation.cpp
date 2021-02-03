#include <iostream>
#include <Eigen/Dense>
#include "Particle.h"
#include "Simulation.h"
#include "Grid.h"

#include <filesystem>
#include <fstream>
#include <cmath>
#include <execution>
#include <chrono>
#include <algorithm>
#include <thread>
#include <future>
#include "MaxwellBoltzmann.h"
#include "Wall.h"
#include "ConstForce.h"
#include "CircularObstacle.h"

#define _USE_MATH_DEFINES
#include <math.h>

using namespace Eigen;
using namespace MPCD;
using namespace std::chrono;

MPCD::Simulation::Simulation(bool draw, bool particleDrawing, int stationaryT) :
	_pipe(ConstForce(Eigen::Vector2d(MPCD::Constants::const_force, 0))),
	_draw(draw),
	_drawParticles(particleDrawing),
	_stationaryT(stationaryT)
{
	_timestep = MPCD::Constants::time_lapse;
	_mdTimestep = MPCD::Constants::md_timestep;
	std::cout << "Streaming timestep: " << _timestep << "\n";
	std::cout << "Monomer timestep: " << _mdTimestep << "\n";
	
	_t = 0;

	setup();
}

void Simulation::setup() {
	int av_particles = MPCD::Constants::average_particles_per_cell;
	double cell_dim = MPCD::Constants::cell_dim;
	double x_0 = MPCD::Constants::x_0;
	double x_max = MPCD::Constants::x_max;
	double width = x_max - x_0;
	double y_0 = MPCD::Constants::y_0;
	double y_max = MPCD::Constants::y_max;
	double height = y_max - y_0;
	double timelapse = MPCD::Constants::time_lapse;

	int number = (int)std::round((double)av_particles * (width / cell_dim) * (height / cell_dim)); // will be a func of Grid

	//std::vector<std::shared_ptr<Obstacle>> obstacles = setUpObstacles(y_0, y_max);
	std::vector<Wall> walls;
	Wall lower(y_0);
	Wall upper(y_max);
	walls.push_back(lower);
	walls.push_back(upper);

	Eigen::Vector2d acceleration(MPCD::Constants::const_force, 0);
	

	const double center_to_center_spacing = MPCD::Obstacles::center_to_center_spacing;
	const int num_per_row = MPCD::Obstacles::num_per_row;
	const double radius = MPCD::Obstacles::radius;
	const double y_upper = MPCD::Obstacles::y_center_upper;
	const double y_lower = MPCD::Obstacles::y_center_lower;

	std::string filename("../../Analysis/Data/constants_");
	std::string obstacles_fp("obstacles");
	std::string csv(".csv");
	std::ofstream outFile(filename + obstacles_fp + csv);
	outFile << "x,y,r" << "\n"; // header columns

	std::vector<CircularObstacle> obstacles;
	if (_addObstacles) {
		for (int i = 0; i < num_per_row; i++) {

			std::cout << "X offset: " << radius * i + center_to_center_spacing * i << std::endl;
			Eigen::Vector2d center_upper(radius * i + center_to_center_spacing * i, y_upper);
			CircularObstacle circ_upper(center_upper, radius);
			writeCirclePositionToOut(outFile, center_upper, radius);
			obstacles.push_back(CircularObstacle(circ_upper));

			Eigen::Vector2d center_lower(radius * i + center_to_center_spacing * i, y_lower);
			CircularObstacle circ_lower(center_lower, radius);
			writeCirclePositionToOut(outFile, center_lower, radius);
			obstacles.push_back(CircularObstacle(circ_lower));
			
		}
	}
	outFile.close();

	_grid.setupCells(obstacles, walls);

	setUpParticles(number, x_0, x_max, y_0, y_max, obstacles);
	_pipe.setObstacles(obstacles, walls);

	if (_draw) {
		writeConstantsToOut(timelapse, width, height, cell_dim, av_particles);
	}
	std::cout << "Setup complete." << std::endl;
}



void MPCD::Simulation::writeCirclePositionToOut(std::ofstream& outFile, Eigen::Vector2d center_pos, double radius) {
	outFile << center_pos[0] << "," << center_pos[1] << "," << radius << std::endl;
}

void Simulation::setUpParticles(int number, double x_0, double x_max, double y_0, double y_max, std::vector<CircularObstacle> obstacles) {
	_particles.reserve(number);

	//dont worry the numbers are just seeds 
	std::mt19937_64 xGen{std::random_device()()};
	std::uniform_real_distribution<double> unifX(x_0, x_max);
	std::mt19937_64 yGen{ std::random_device()() };
	std::uniform_real_distribution<double> unifY(y_0, y_max);


	// MAXWELL BOLTZMANN
	double mass = MPCD::Constants::particle_mass; // h2o kg mass
	double mean = 0;
	double temperature = MPCD::Constants::temperature;

	double max_x_vel = std::sqrt(2 * MPCD::Constants::k_boltzmann * MPCD::Constants::temperature / MPCD::Constants::particle_mass);
	double max_y_vel = std::sqrt(2 * MPCD::Constants::k_boltzmann * MPCD::Constants::temperature / MPCD::Constants::particle_mass);
	MaxwellBoltzmann mb_vel(mean, temperature, mass);

	for (int i = 0; i < number; i++) {
		double xs_x = unifX(xGen);
		double xs_y = unifY(yGen);
		Eigen::Vector2d pos(xs_x, xs_y);
		

		Eigen::Vector2d v = mb_vel.next();

		Particle p(mass, pos, v);
		
		if (isInBoundsOfAnObstacle(p, obstacles)) {
			while (isInBoundsOfAnObstacle(p, obstacles)) {
				p.correctPosition(Eigen::Vector2d(unifX(xGen), unifY(yGen)));
			}
		}

		_particles.push_back(p);
	}
}

void MPCD::Simulation::setUpMonomers()
{
	std::uniform_real_distribution<double> mt_y{MPCD::Constants::y_0 + MPCD::Constants::monomer_diameter, MPCD::Constants::y_max - MPCD::Constants::monomer_diameter};
	std::random_device rd{};
	std::mt19937_64 gen{ rd() };
	
	double start_y = mt_y(gen);
	Eigen::Vector2d start(MPCD::Obstacles::x_end + MPCD::Constants::monomer_diameter, start_y);
		
	double end_y = mt_y(gen);
	Eigen::Vector2d end(MPCD::Obstacles::x_end + 120 * MPCD::Constants::monomer_diameter, end_y);

	Eigen::Vector2d step = (end - start) / MPCD::Constants::num_monomers;
	MaxwellBoltzmann mb(0, MPCD::Constants::temperature, MPCD::Constants::monomer_mass);
	for (int n = 0; n < MPCD::Constants::num_monomers; n++) {
		Eigen::Vector2d monomer_position = start + n * step;
		Eigen::Vector2d monomer_velocity = mb.next();
		_monomers.push_back(Monomer(
			MPCD::Constants::monomer_mass,
			monomer_position,
			monomer_velocity,
			MPCD::Constants::monomer_diameter
		));
	}

}

bool MPCD::Simulation::isInBoundsOfAnObstacle(Body& b, std::vector<CircularObstacle> obstacles) {
	for (auto& o : obstacles) {
		if (o.isInBounds(b)) return true;
	}
	return false;
}


void MPCD::Simulation::writeConstantsToOut(double timelapse, double width, double height, double cell_dim, int averageParticlesPerCell) {
	int num_hypothetical_x_cells = (int)std::round(width / cell_dim);
	int num_hypothetical_y_cells = (int)std::round(height / cell_dim);

	if (width >= height) {
		assert(num_hypothetical_x_cells >= num_hypothetical_y_cells);
	}
	else {
		assert(num_hypothetical_x_cells < num_hypothetical_y_cells);
	}

	double epsilon = 0.05;

	std::filesystem::path cwd = std::filesystem::current_path();

	std::string filename("constants");
	std::string csv(".csv");

	std::ofstream outFile("../../Analysis/Data/" + filename + csv);

	double particle_mass = MPCD::Constants::particle_mass;
	double k_BT = MPCD::Constants::k_boltzmann * MPCD::Constants::temperature;
	outFile << "stationary_t,avg,delta_t,a,w,h,num_particles,m,k_BT,g" << "\n"; // header columns
	outFile << _stationaryT << averageParticlesPerCell << timelapse << "," << cell_dim << "," << width << "," << height << "," << num_hypothetical_x_cells*num_hypothetical_y_cells*averageParticlesPerCell << "," << particle_mass << "," << k_BT << "," << MPCD::Constants::const_force << std::endl;
	outFile.close();
}

void MPCD::Simulation::timestep()
{	
	if (_t == _stationaryT) {
		setUpMonomers();
	}
	if (_t >= _stationaryT) {
		verlet();
	}

	streamingStep();

	_grid.shift();
	_grid.updateCoordinates(_particles, _monomers);

	collisionStep();

	_grid.undoShift();
	
	_t += 1;
}

void MPCD::Simulation::verlet() {
	
	bool draw = true;
	_pipe.verlet(_monomers, draw, _t);
}

void MPCD::Simulation::streamingStep() {

	bool draw = false;
	if (_t == _stationaryT-1) {
		draw = true;
	}
	_pipe.stream(_particles, _timestep, draw, _t);
	
}

void MPCD::Simulation::collisionStep() {

	bool draw = false;
	if (_t > _stationaryT-1) {
		draw = true;
	}
	_grid.calculate(draw, _t);
	_grid.collision(_particles, _monomers);
	
	
}

