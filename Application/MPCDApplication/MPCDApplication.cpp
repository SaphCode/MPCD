// MPCDApplication.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <Eigen/Dense>

#include <iostream>
#include "Constants.h"
#include "Simulation.h"
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>

using namespace Eigen;
using namespace MPCD;

bool cancel = false;
bool draw = true;
bool particleDrawing = false;
int stationaryT = 3000;

Simulation sim(draw, particleDrawing, stationaryT);

int main()
{

	int t = 2999;
	
	//sim.setup();
	sim.loadCheckpoint(2999, "G:/Bachelor/Data/particles_av10_timestep02999.csv"); // , "G:/Bachelor/Data/cells_av10_timestep02999.csv"

	auto start = std::chrono::high_resolution_clock::now();
	while (!cancel) {
		t += 1;
		if (t % 100 == 0) {
			auto finish = std::chrono::high_resolution_clock::now();
			auto seconds = std::chrono::duration_cast<std::chrono::seconds>(finish - start);
			std::cout << "Timestep: " << t << " (" << seconds.count() << "s)" << "\n";
			start = std::chrono::high_resolution_clock::now();
		}
		sim.timestep();
	}
	
}
