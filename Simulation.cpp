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
#include <map>
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
	_time_lapse = MPCD::Constants::time_lapse;
	_cell_dim = MPCD::Constants::Grid::cell_dim;
	_particles = particles;
	Xoshiro rg_shift_x(-MPCD::Constants::Grid::max_shift, MPCD::Constants::Grid::max_shift);
	Xoshiro rg_shift_y(-MPCD::Constants::Grid::max_shift, MPCD::Constants::Grid::max_shift);
	Xoshiro rg_angle(0.0, 2 * M_PI);
	Xoshiro rg_sign(-1, 1);
	_rg_shift_x = rg_shift_x;
	_rg_shift_y = rg_shift_y;
	_rg_angle = rg_angle;
	_rg_sign = rg_sign;
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
}

void MPCD::Simulation::timestep(int t)
{
	Vector2d shift(_rg_shift_x.next(), _rg_shift_y.next());

	auto t1 = std::chrono::high_resolution_clock::now();
	moveShiftPrepare(shift, t);
	auto t2 = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
	std::cout << "Parallelized: " << duration << std::endl;

	t1 = std::chrono::high_resolution_clock::now();
	calculateCellQuantities();
	t2 = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
	std::cout << "CalculateCellQuantities: " << duration << std::endl;
	
	t1 = std::chrono::high_resolution_clock::now();
	updateVelocity(shift);
	t2 = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
	std::cout << "updateVelocity: " << duration << std::endl;
	//calculateCellQuantities(_totalCellVelocities, _numCellParticles);
	//updateVelocity(shift);

	/*std::map<std::pair<int, int>, Eigen::Vector2d> meanCellVelocities;
	if (_draw) {
		draw(t);
	}
	*/
	reset(t);
}

/* O(N) */
void MPCD::Simulation::moveShiftPrepare(Eigen::Vector2d shift, int t) {
	std::filesystem::path cwd;
	std::stringstream s;
	std::stringstream av;
	std::string filename;
	
	cwd = std::filesystem::current_path();
	s << std::setfill('0') << std::setw(_w) << t;
	av << "av" << MPCD::Constants::Grid::average_particles_per_cell << "_";
	filename = cwd.string() + l_data + "particles_" + av.str() + "timestep" + s.str() + ".csv";
	std::ofstream outFile(filename);
	

	if (_draw) {
		std::stringstream header("x,y,vx,vy");
		outFile << header.str() << '\n';
	}

	std::mutex m;
	std::for_each(std::execution::par, _particles.begin(), _particles.end(), [&](auto& particle) {
		if (_draw) {
			Vector2d pos = particle.getPosition();
			Vector2d vel = particle.getVelocity();
			m.lock();
			outFile << pos[0] << "," << pos[1] << "," << vel[0] << "," << vel[1] << "\n";
			m.unlock();
		}

		particle.move(); // move const
		particle.shift(shift); // shift // const
		//for (auto& const particle : _particles) {
		Vector2d velocity = particle.getVelocity();
		Vector2d zero(0, 0);
		Vector2i index = particle.getCellIndex();
		std::pair indexp = std::make_pair(index[0], index[1]);

		bool exists = (_totalCellVelocities.count(indexp) > 0) ? true : false;
		//std::lock_guard lock(m);
		if (exists) {
			_totalCellVelocities[indexp] += velocity;
			_numCellParticles[indexp] += 1;
		}
		else {
			m.lock();
			_totalCellVelocities.insert(std::make_pair(indexp, velocity));
			_meanCellVelocities.insert(std::make_pair(indexp, zero));
			_numCellParticles.insert(std::make_pair(indexp, 1));
			m.unlock();
		}
		});
		
	
	/*
	for (auto it = _particles.begin(); it != _particles.end(); ++it) {
		it->move(); // move const
		it->shift(shift); // shift // const
		Vector2d velocity = it->getVelocity();
		Vector2d zero(0, 0);
		Vector2i index = it->getCellIndex();
		std::pair indexp = std::make_pair(index[0], index[1]);

		bool newlyInserted;
		newlyInserted = _totalCellVelocities.insert(std::make_pair(indexp, velocity)).second;
		if (!newlyInserted) {
			_totalCellVelocities[indexp] += velocity;
		} else {
			_meanCellVelocities.insert(std::make_pair(indexp, zero));
		}
		_numCellParticles[indexp] += 1;
		*/
		//it->getCellIndex();
		/*
		Eigen::Vector2d pos = it->getPosition();
		if (pos[0] > max[0]) {
			max[0] = pos[0];
		}
		else if (pos[0] < min[0]) {
			min[0] = pos[0];
		}
		if (pos[1] > max[1]) {
			max[1] = pos[1];
		}
		else if (pos[1] < min[1]) {
			min[1] = pos[1];
		}
		*/
	//}
	/*
	double cell_dim = MPCD::Constants::Grid::cell_dim;

	int min_col = std::floor(min[0] / cell_dim);
	int max_col = std::floor(max[0] / cell_dim);
	int cols = max_col - min_col;
	assert(cols > 0);
	int min_row = std::floor(min[1] / cell_dim);
	int max_row = std::floor(max[1] / cell_dim);
	int rows = max_row - min_row;
	assert(rows > 0);
	}
	*/
	
}



/* O(N) (because num cells should be proportional to particles) */
void MPCD::Simulation::calculateCellQuantities() {
	Vector2d zero(0, 0);
	assert(_totalCellVelocities.size() == _meanCellVelocities.size());
	assert(_meanCellVelocities.size() == _numCellParticles.size());
	std::for_each(std::execution::par, _totalCellVelocities.begin(), _totalCellVelocities.end(), [&](std::pair<std::pair<int, int>, Vector2d> entry) {
		std::pair<int, int> key = entry.first;
		Vector2d val = entry.second;
		if (val == zero) {
			_meanCellVelocities[key] = zero;
		}
		else {
			int numParticles = _numCellParticles[key];
			assert(numParticles != 0);
			_meanCellVelocities[key] = val / numParticles;
		}
		_cellRotationAngle[key] = _rg_angle.next();
		});
}

/* O(N) */
void MPCD::Simulation::updateVelocity(Eigen::Vector2d shift) {
	//std::cout << "--Updating velocities .." << "\n";
	std::for_each(std::execution::par, _particles.begin(), _particles.end(), [&](auto& particle) {
		Vector2i index = particle.getCellIndex(); // const
		double s = _rg_sign.next();
		int sign;
		if (s < 0) {
			sign = -1;
		}
		else if (s >= 0) {
			sign = 1;
		}
		std::pair indexp = std::make_pair(index[0], index[1]);
		particle.updateVelocity(_meanCellVelocities[indexp], sign * _cellRotationAngle[indexp]); // const
		particle.shift(-shift);
		});

}

void MPCD::Simulation::reset(int t) {
	std::filesystem::path cwd;
	std::stringstream s;
	std::stringstream av;
	std::string filename;

	cwd = std::filesystem::current_path();
	s << std::setfill('0') << std::setw(_w) << t;
	av << "av" << MPCD::Constants::Grid::average_particles_per_cell << "_";
	filename = cwd.string() + l_data + "cells_" + av.str() + "timestep" + s.str() + ".csv";
	std::ofstream outFile(filename);

	if (_draw) {
		std::stringstream header("i,j,vx,vy");
		outFile << header.str() << '\n';
	}

	Vector2d reset_vel(0, 0);
	assert(_totalCellVelocities.size() == _meanCellVelocities.size());
	assert(_meanCellVelocities.size() == _numCellParticles.size());
	std::mutex m;
	std::for_each(std::execution::par, _totalCellVelocities.begin(), _totalCellVelocities.end(), [&](std::pair<std::pair<int, int>, Vector2d> entry) {
		std::pair<int, int> key = entry.first;
		Vector2d val = entry.second;
		if (_draw) {
			m.lock();
			outFile << key.first << "," << key.second << "," << val[0] << "," << val[1] << "\n";
			m.unlock();
		}
		_totalCellVelocities[key] = reset_vel;
		_meanCellVelocities[key] = reset_vel;
		_numCellParticles[key] = 0;
		});
}

/*
void MPCD::Simulation::draw(int t) {
	std::filesystem::path cwd = std::filesystem::current_path();
	std::cout << cwd.string() << "," << l_data << std::endl;
	drawParticles(t);
	drawMeanAndFrequencies(t);
}

void MPCD::Simulation::drawParticles(int t) {
	std::stringstream s;
	std::stringstream av;
	s << std::setfill('0') << std::setw(_w) << t;
	av << "_av" << MPCD::Constants::Grid::average_particles_per_cell;
	std::string filename = "particles_timestep_" + s.str() + av.str() + ".csv";
	std::filesystem::path cwd = std::filesystem::current_path();
	Out out(cwd.string() + l_data);
	out.writeToOut(_particles, filename);
}

void MPCD::Simulation::drawMeanAndFrequencies(int t) {
	std::stringstream s;
	std::stringstream av;
	s << std::setfill('0') << std::setw(_w) << t;
	av << "_av" << MPCD::Constants::Grid::average_particles_per_cell;
	std::string filename_mean = "cellvalues_timestep_" + av.str() + s.str() + ".csv";
	std::string filename_frequ = "cell_frequencies_" + av.str() + s.str() + ".csv";
	std::filesystem::path cwd = std::filesystem::current_path();
	Out out(cwd.string() + l_data);
	out.writeToOut(_meanCellVelocities, filename_mean, "i,j,vx,vy");
	out.writeToOut(_numCellParticles, filename_frequ, "i,j,num");
}
*/