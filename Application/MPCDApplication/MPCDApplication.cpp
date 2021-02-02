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

Simulation sim(draw, particleDrawing);

void interrupt_handler(sig_atomic_t s) {
	int input = -1;
	std::cout << "Caught signal: " << s;
	std::cout << "Do you want to interrupt the process? (1/0) 1 is yes, 0 is no\n";
	while (input != 0 && input != 1) {
		std::cin >> input;
	}

	if (input == 1) {
		exit(1);
	}
	else if (input == 0) {
		cancel = false;
	}
	

	std::cout << "Set drawing interval to an integer (from 1 to 100):\n";
	int drawInterval = -1;
	while (drawInterval < 1 && drawInterval > 100) {
		std::cin >> drawInterval;
	}
	sim.setDrawingInterval(drawInterval);
}

int main()
{
	signal(SIGINT, interrupt_handler);

	int t = 0;
	while (!cancel) {
		t += 1;
		if (t % 100 == 0) {
			std::cout << "Timestep: " << t << "\n";
		}
		sim.timestep();
	}
	
}
