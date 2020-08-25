#pragma once

#ifndef MPCD_H
#define MPCD_H

#include "Particle.h"
#include "Grid.h"
#include "Xoshiro.h"

namespace MPCD {
	/* One timestep */
	void timestep(std::vector<MPCD::Particle>& particles, MPCD::Grid g, double time_lapse, Xoshiro& rg_shift_x, Xoshiro& rg_shift_y, Xoshiro& rg_angle);
}

#endif // !MPCD_H
