// MPCDApplication.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "Constants.h"
#include "Simulation.h"

using namespace Eigen;
using namespace MPCD;

int main()
{
	bool draw = true;
	bool particleDrawing = false;
	int drawInterval = 50;
	int drawLast = 1000;
	Simulation sim(draw, drawInterval, drawLast, particleDrawing);

	int timesteps = MPCD::Constants::timesteps;

	for (int t = 0; t < timesteps; t++) {
		//if (t % drawInterval == 0) {
		std::cout << "Timestep: " << t << "\n";
		//}
		sim.timestep();
	}	
	
}
