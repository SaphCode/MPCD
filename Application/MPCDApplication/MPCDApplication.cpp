// MPCDApplication.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <Eigen/Dense>
#include "Particle.h"
#include "Constants.h"
#include "Xoshiro.h"
#include "Grid.h"
#include "Simulation.h"

#include <fstream>
#include <filesystem>

#define _USE_MATH_DEFINES
#include <math.h>

using namespace Eigen;
using namespace MPCD;



int main()
{
	bool draw = true;
	Simulation sim(draw);

	int timesteps = MPCD::Constants::timesteps;

	for (int t = 0; t < timesteps; t++) {
		std::cout << "Timestep: " << t << "\n";
		sim.timestep();
	}	
	
}
// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
