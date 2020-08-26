#include "pch.h"
#include "Out.h"

#include <filesystem>
#include <fstream>

using namespace Eigen;

Out::Out(std::string location) {
	_location = location;
}

Out::Out() {
	_location = std::filesystem::current_path().string();
}

Out::~Out() {}

void Out::writeToOut(std::vector<MPCD::Particle> particles, std::string filename) {
	std::string header = "x,y,vx,vy";
	std::ofstream outFile(_location + "//" + filename);
	outFile << header << std::endl;
	for (auto& p : particles) {
		Vector2d pos = p.getPosition();
		Vector2d vel = p.getVelocity();
		outFile << pos[0] << "," << pos[1] << "," << vel[0] << "," << vel[1] << std::endl;
	}
	outFile.close();
}

void Out::writeToOut(std::vector<double> numbers, std::string filename, std::string header) {
	std::ofstream outFile(_location + "//" + filename);
	outFile << header << std::endl;
	for (auto& n : numbers) {
		outFile << n << std::endl;
	}
}
void Out::writeToOut(std::vector<Eigen::Vector2d> vectors, std::string filename, std::string header) {
	std::ofstream outFile(_location + "//" + filename);
	outFile << header << std::endl;
	for (auto& v : vectors) {
		outFile << v[0] << "," << v[1];
	}
	outFile.close();
}

void Out::writeToOut(std::map<int, Eigen::Vector2d> map, std::string filename, std::string header) {
	std::ofstream outFile(_location + "//" + filename);
	outFile << header << std::endl;
	for (auto it = map.begin(); it != map.end(); ++it) {
		outFile << it->second[0] << "," << it->second[1] << std::endl;
	}
	outFile.close();
}