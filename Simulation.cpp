#include <iostream>
#include <Eigen/Dense>
#include "Particle.h"
#include "Simulation.h"
#include "Grid.h"

#include <filesystem>
#include <fstream>
#include <cmath>
#include "Out.h"
#include "Locations.h"
#include <execution>
#include <chrono>
#include <algorithm>
#include <thread>
#include <future>
#include "MaxwellBoltzmann.h"
#include "Wall.h"
#include "Obstacle.h"
#include "ConstForce.h"

#define _USE_MATH_DEFINES
#include <math.h>

using namespace Eigen;
using namespace MPCD;
using namespace std::chrono;

MPCD::Simulation::Simulation(bool draw)
{
	_timelapse = MPCD::Constants::time_lapse;
	_draw = draw;
	int timesteps = MPCD::Constants::timesteps;

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

	int number = av_particles * (width / cell_dim) * (height / cell_dim); // will be a func of Grid

	//std::vector<std::shared_ptr<Obstacle>> obstacles = setUpObstacles(y_0, y_max);

	Wall lower(y_0, false);
	Wall upper(y_max, true);
	Eigen::Vector2d acceleration(MPCD::Constants::acceleration_const, 0);
	ConstForce constForce(acceleration);

	_obstacles.push_back(std::make_shared<Wall>(lower));
	_obstacles.push_back(std::make_shared<Wall>(upper));

	_interactors.push_back(std::make_shared<Wall>(lower));
	_interactors.push_back(std::make_shared<Wall>(upper));
	_interactors.push_back(std::make_shared<ConstForce>(constForce));

	setUpParticles(number, x_0, x_max, y_0, y_max);
	_pipe.setObstacles(_obstacles);

	if (_draw) {
		int timesteps = MPCD::Constants::timesteps;
		writeConstantsToOut(timelapse, width, height, cell_dim, av_particles, timesteps);
	}

}

void Simulation::setUpParticles(int number, double x_0, double x_max, double y_0, double y_max) {
	_particles.reserve(number);

	//dont worry the numbers are just seeds 
	Xoshiro xs_xpos(x_0, x_max);
	Xoshiro xs_ypos(y_0, y_max);

	// MAXWELL BOLTZMANN
	double mass = MPCD::Constants::particle_mass; // h2o kg mass
	double mean = 0;
	double temperature = MPCD::Constants::temperature;

	double max_x_vel = std::sqrt(2 * MPCD::Constants::k_boltzmann * MPCD::Constants::temperature / MPCD::Constants::particle_mass);
	double max_y_vel = std::sqrt(2 * MPCD::Constants::k_boltzmann * MPCD::Constants::temperature / MPCD::Constants::particle_mass);
	Xoshiro xvel(-max_x_vel, max_x_vel);
	Xoshiro yvel(-max_y_vel, max_y_vel);
	//MaxwellBoltzmann mb_vel(mean, temperature, mass);

	for (int i = 0; i < number; i++) {
		double xs_x = xs_xpos.next();
		double xs_y = xs_ypos.next();
		Eigen::Vector2d pos(xs_x, xs_y);

		double vx = xvel.next();
		double vy = yvel.next();
		Eigen::Vector2d vel(vx, vy);

		Particle p(mass, pos, vel);

		_particles.push_back(p);
	}
}

void MPCD::Simulation::writeConstantsToOut(double timelapse, double width, double height, double cell_dim, int averageParticlesPerCell, int timesteps) {
	int num_hypothetical_x_cells = std::round(width / cell_dim);
	int num_hypothetical_y_cells = std::round(height / cell_dim);

	if (width >= height) {
		assert(num_hypothetical_x_cells >= num_hypothetical_y_cells);
	}
	else {
		assert(num_hypothetical_x_cells < num_hypothetical_y_cells);
	}

	double epsilon = 0.05;

	std::filesystem::path cwd = std::filesystem::current_path();

	std::string filename("//Data//constants_");
	std::string num("av" + std::to_string(averageParticlesPerCell));
	std::string csv(".csv");

	std::ofstream outFile(cwd.string() + filename + num + csv);

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
	std::cout << "Streaming Step: " << duration << std::endl;

	t1 = std::chrono::high_resolution_clock::now();
	collisionStep();
	t2 = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
	std::cout << "Collision Step: " << duration << std::endl;
}

/* O(N) */
void MPCD::Simulation::streamingStep() {
	
	std::filesystem::path cwd;
	std::stringstream s;
	std::stringstream av;
	std::string filename;
	
	cwd = std::filesystem::current_path();
	s << std::setfill('0') << std::setw(_w) << _t;
	av << "av" << _grid.getAverageParticlesPerCell() << "_";
	filename = cwd.string() + l_data + "particles_" + av.str() + "timestep" + s.str() + ".csv";
	std::ofstream outFile(filename);
	

	if (_draw) {
		std::stringstream header("x,y,vx,vy");//,vx,vy"
		outFile << header.str() << '\n';
	}

	_pipe.stream(_particles, _interactors, _timelapse, _draw, outFile);
	


	//reduction(+ : grid)
	//#pragma omp parallel
	//#pragma omp for
	/*const auto processor_count = std::thread::hardware_concurrency();
	int size = _particles.size();
	std::vector<Grid> grids;
	int endBefore = -1;
	std::cout << "Size of particles: " << size << std::endl;
	/*std::vector<std::thread> threads;
	for (int p = 0; p < processor_count; p++) {
		// make new thread
		//, _particles, shift
		Grid g;
		grids.push_back(g);
		double sizeToProcess = size / processor_count;
		int start = std::round(p * sizeToProcess);
		int end = std::round((p + 1) * sizeToProcess);
		if (start == endBefore) {
			start++;
		}
		//std::cout << "Start: " << start << ", End: " << end << std::endl;
		Vector2d shiftT = shift;
		std::vector<Particle> particlesT(_particles.begin() + start, _particles.begin() + end);
		std::thread t(&MPCD::calculateGrid, std::ref(g), std::ref(particlesT), shift, _timelapse);
		threads.push_back(std::move(t));
		//t.detach();
		//calculateGrid(g, particlesT, shift, _timelapse);
		endBefore = end;
	}
	for (int i = 0; i < threads.size(); i++) {
		threads[i].join();
	}
	threads.clear();
	Grid realGrid;
	for (const auto& g : grids) {
		realGrid = realGrid + g;
	}
	std::pair<int, int> coords({ 0,0 });
	std::cout << realGrid._cells[coords];
	/*std::for_each(std::execution::par, _particles.begin(), _particles.end(), [_timelapse, shift, grid](Particle p) {
		p.stream(_timelapse);
		p.shift(shift);
		grid.insert(p);
		});
	_grid = grid;*/
}

/* O(N) */
void MPCD::Simulation::collisionStep() {
	std::filesystem::path cwd;
	std::stringstream s;
	std::stringstream av;
	std::string filename;

	cwd = std::filesystem::current_path();
	s << std::setfill('0') << std::setw(_w) << _t;
	av << "av" << _grid.getAverageParticlesPerCell() << "_";
	filename = cwd.string() + l_data + "cells_" + av.str() + "timestep" + s.str() + ".csv";
	std::ofstream outFile(filename);

	if (_draw) {
		std::stringstream header("i,j,meanX,meanY,num");//,vx,vy"
		outFile << header.str() << '\n';
	}


	_grid.shift();
	_grid.updateCoordinates(_particles);
	_grid.collision(_draw, outFile);
	_grid.undoShift();
}

