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

	int t = 0;
	while (!cancel) {
		t += 1;
		if (t % 100 == 0) {
			std::cout << "Timestep: " << t << "\n";
		}
		sim.timestep();
	}
	
}
