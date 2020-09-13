// MPCD.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <Eigen/Dense>
#include "Particle.h"
#include "Simulation.h"
#include "Grid.h"
#include "Constants.h"
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


#define _USE_MATH_DEFINES
#include <math.h>

using namespace Eigen;
using namespace MPCD;
using namespace std::chrono;

MPCD::Simulation::Simulation(bool draw) {
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

	int av_particles = MPCD::Constants::average_particles_per_cell;
	double cell_dim = MPCD::Constants::cell_dim;
	double x_0 = MPCD::Constants::x_0;
	double x_max = MPCD::Constants::x_max;
	double y_0 = MPCD::Constants::y_0;
	double y_max = MPCD::Constants::y_max;

	int number = av_particles * ((x_max - x_0) / cell_dim) * ((y_max - y_0) / cell_dim); // will be a func of Grid
	
	std::vector<Obstacle> obstacles;

	std::vector<Particle> particles;
	particles.reserve(number);

	//dont worry the numbers are just seeds 
	Xoshiro xs_xpos(x_0, x_max);
	Xoshiro xs_ypos(y_0, y_max);

	// MAXWELL BOLTZMANN
	double mass = 2.988e-26; // h2o kg mass
	double mean = 0;
	double temperature = 309.15; // = 36Celsius
	MaxwellBoltzmann mb_vel(mean, temperature, mass);

	for (int i = 0; i < number; i++) {
		double xs_x = xs_xpos.next();
		double xs_y = xs_ypos.next();
		Eigen::Vector2d pos(xs_x, xs_y);

		Eigen::Vector2d vel = mb_vel.next();

		Particle p(pos, vel);

		particles.push_back(p);
	}
	
	_pipe.setParticles(particles);
	_pipe.setObstacles(obstacles);
}

void MPCD::Simulation::timestep()
{
	_t += 1;
	Vector2d shift(_rgShiftX.next(), _rgShiftY.next());
	
	auto t1 = std::chrono::high_resolution_clock::now();
	streamingStep(shift);
	auto t2 = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
	std::cout << "Streaming Step: " << duration << std::endl;

	t1 = std::chrono::high_resolution_clock::now();
	collisionStep(shift);
	t2 = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
	std::cout << "Collision Step: " << duration << std::endl;
}

/* O(N) */
void MPCD::Simulation::streamingStep(Eigen::Vector2d shift) {
	
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

	_pipe.stream(_timelapse, _draw, outFile);
	_grid.updateCoordinates(_pipe.getParticles());


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
void MPCD::Simulation::collisionStep(Eigen::Vector2d shift) {
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
	_grid.collision(_draw, outFile);
}
//
void MPCD::calculateGrid(Grid& grid, std::vector<Particle>& particles, const Eigen::Vector2d shift, const double timelapse) {
	//for (auto& p : particles) {
	//	p.stream(timelapse);
	//	p.shift(shift);
	//	grid.insert(p);
	//}
}
