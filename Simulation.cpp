#include <iostream>
#include <Eigen/Dense>
#include "Particle.h"
#include "Simulation.h"
#include "Grid.h"

#include <filesystem>
#include <fstream>
#include <cmath>
#include "Locations.h"
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

MPCD::Simulation::Simulation(bool draw) :
	_pipe(ConstForce(Eigen::Vector2d(MPCD::Constants::acceleration_const, 0))),
	_draw(draw)
{
	_timelapse = MPCD::Constants::time_lapse;
	
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

	Eigen::Vector2d acceleration(MPCD::Constants::acceleration_const, 0);
	

	const double x_offset = MPCD::Obstacles::x_offset;
	const int num_circular_obstacles = MPCD::Obstacles::num;
	const double x_dist = MPCD::Obstacles::x_dist;
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
		for (int i = 0; i < num_circular_obstacles / 2; i++) {

			std::cout << "X offset: " << x_offset + radius + 2 * radius * i + x_dist * i << std::endl;
			Eigen::Vector2d center_upper(x_offset + radius + x_dist * i, y_upper);
			CircularObstacle circ_upper(center_upper, radius);
			writeCirclePositionToOut(outFile, center_upper, radius);
			obstacles.push_back(CircularObstacle(circ_upper));

			Eigen::Vector2d center_lower(x_offset + radius + x_dist * i, y_lower);
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
		int timesteps = MPCD::Constants::timesteps;
		writeConstantsToOut(timelapse, width, height, cell_dim, av_particles, timesteps);
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

bool MPCD::Simulation::isInBoundsOfAnObstacle(Body& b, std::vector<CircularObstacle> obstacles) {
	for (auto& o : obstacles) {
		if (o.isInBounds(b)) return true;
	}
	return false;
}


void MPCD::Simulation::writeConstantsToOut(double timelapse, double width, double height, double cell_dim, int averageParticlesPerCell, int timesteps) {
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

	std::string filename("constants_");
	std::string num("av" + std::to_string(averageParticlesPerCell));
	std::string csv(".csv");

	std::ofstream outFile("../../Analysis/" + l_data + filename + num + csv);

	double particle_mass = MPCD::Constants::particle_mass;
	double k_BT = MPCD::Constants::k_boltzmann * MPCD::Constants::temperature;
	outFile << "timesteps,time_lapse,cell_dim,width,height,average_particles_per_cell,total_number_of_particles,particle_mass,k_BT" << "\n"; // header columns
	outFile << timesteps << "," << timelapse << "," << cell_dim << "," << width << "," << height << "," << averageParticlesPerCell << "," << num_hypothetical_x_cells*num_hypothetical_y_cells*averageParticlesPerCell << "," << particle_mass << "," << k_BT << std::endl;
	outFile.close();
}

void MPCD::Simulation::timestep()
{
	_t += 1;
	
	auto t1 = std::chrono::high_resolution_clock::now();
	streamingStep();
	auto t2 = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
	if (_t % 100 == 0) {
		std::cout << "Streaming Step: " << duration << std::endl;
	}
	

	t1 = std::chrono::high_resolution_clock::now();
	collisionStep();
	t2 = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
	if (_t % 100 == 0) {
		std::cout << "Collision Step: " << duration << std::endl;
	}
	
}

void MPCD::Simulation::streamingStep() {

	bool draw;
	if (_t % 100 == 0 || _t + 50 > Constants::timesteps) {
		draw = true;
	}
	else {
		draw = false;
	}

	_pipe.stream(_particles, _timelapse, draw, _t);
	
}

void MPCD::Simulation::collisionStep() {

	_grid.shift();
	_grid.updateCoordinates(_particles);

	bool draw;
	if (_t % 100 == 0 || _t + 100 > Constants::timesteps) {
		draw = true;
	}
	else {
		draw = false;
	}

	_grid.calculate(draw, _t);
	_grid.collision(_particles);
	_grid.undoShift();
	
}

