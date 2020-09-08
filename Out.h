#pragma once

#ifndef OUT_H
#define OUT_H

#include <string>
#include <vector>
#include "Particle.h"
#include <Eigen/Dense>
#include <map>
#include <boost/unordered_map.hpp>

static const std::string data = "//Data";
static const std::string rng = "//RNG";

class Out
{
private:
	std::string _location;
	//void init();
public:
	/* @param location: the desired saving location */
	Out(std::string location);
	/* Takes cwd*/
	Out();
	~Out();
	void writeToOut(std::vector<MPCD::Particle> particles, std::string filename);
	void writeToOut(std::vector<double> numbers, std::string filename, std::string header);
	void writeToOut(std::vector<Eigen::Vector2d> vectors, std::string filename, std::string header);
	void writeToOut(boost::unordered::unordered_map<std::pair<int, int>, Eigen::Vector2d> map, std::string filename, std::string header);
	void writeToOut(boost::unordered::unordered_map<std::pair<int, int>, int> frequencies, std::string filename, std::string header);
};


#endif // !OUT_H