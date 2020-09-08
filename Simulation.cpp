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
		std::stringstream header("x,y,vx,vy");
		outFile << header.str() << '\n';
	}

	Grid grid;
	//#pragma omp parallel for reduction(+ : grid)
	for (auto& p : _particles) {
		p.stream(_timelapse);
		p.shift(shift);
		grid.insert(p);
	}
	_grid = grid;
}

/* O(N) */
void MPCD::Simulation::collisionStep(Eigen::Vector2d shift) {
	std::for_each(std::execution::par, _grid._cells.begin(), _grid._cells.end(), [shift](std::pair<std::pair<int, int>, Cell> entry) {
		entry.second.collide(shift);
	});
}
