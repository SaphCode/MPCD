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


#define _USE_MATH_DEFINES
#include <math.h>

using namespace Eigen;
using namespace MPCD;
using namespace std::chrono;

MPCD::Simulation::Simulation(std::vector<Particle>& particles, bool draw) {
	_timelapse = MPCD::Constants::time_lapse;
	_particles = particles;
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
	
	_maxShift = MPCD::Constants::Grid::max_shift;
	Xoshiro rgShiftX(-_maxShift, _maxShift);
	Xoshiro rgShiftY(-_maxShift, _maxShift);
	_rgShiftX = rgShiftX;
	_rgShiftY = rgShiftY;
	_t = 0;
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
	av << "av" << MPCD::Constants::Grid::average_particles_per_cell << "_";
	filename = cwd.string() + l_data + "particles_" + av.str() + "timestep" + s.str() + ".csv";
	std::ofstream outFile(filename);
	

	if (_draw) {
		std::stringstream header("x,y");//,vx,vy"
		outFile << header.str() << '\n';
	}
	Grid g;
	double lapse = _timelapse;

	bool draw = _draw;
	std::mutex m;
	std::for_each(std::execution::par, _particles.begin(), _particles.end(),  [&outFile, lapse, shift, &g, &m, draw](Particle p) {
		//for (auto& p : _particles) {
		p.stream(lapse);
		p.shift(shift);
		
		m.lock();
		if (draw) {
			Vector2d pos = p.getPosition();
			outFile << pos[0] << "," << pos[1] << "\n";
		}
		g.insert(p);
		m.unlock();
		});
	_grid = g;


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
	av << "av" << MPCD::Constants::Grid::average_particles_per_cell << "_";
	filename = cwd.string() + l_data + "cells_" + av.str() + "timestep" + s.str() + ".csv";
	std::ofstream outFile(filename);

	if (_draw) {
		std::stringstream header("i,j,meanX,meanY,num");//,vx,vy"
		outFile << header.str() << '\n';
	}
	bool draw = _draw;
	std::mutex m;
	std::for_each(std::execution::par, _grid._cells.begin(), _grid._cells.end(), [&m, &outFile, shift, draw](std::pair<std::pair<int, int>, Cell> entry) {
		entry.second.collide(shift);
		if (draw) {
			m.lock();
			entry.second.draw(entry.first, outFile);
			m.unlock();
		}
	});
}
//
void MPCD::calculateGrid(Grid& grid, std::vector<Particle>& particles, const Eigen::Vector2d shift, const double timelapse) {
	//for (auto& p : particles) {
	//	p.stream(timelapse);
	//	p.shift(shift);
	//	grid.insert(p);
	//}
}
