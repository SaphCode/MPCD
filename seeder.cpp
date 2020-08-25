#include "seeder.h"
#include "Constants.h"

/* This is a fixed-increment version of Java 8's SplittableRandom generator
   See http://dx.doi.org/10.1145/2714064.2660195 and
   http://docs.oracle.com/javase/8/docs/api/java/util/SplittableRandom.html

   It is a very fast generator passing BigCrush, and it can be useful if
   for some reason you absolutely want 64 bits of state.
   This is for generating seeds. should use 2 different algos. */

static uint64_t seed = MPCD::Constants::seed;

uint64_t nextSeed() {
	uint64_t z = (seed += 0x9e3779b97f4a7c15);
	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
	z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
	uint64_t nextSeed = z ^ (z >> 31);
	seed = nextSeed;
	return nextSeed;
}