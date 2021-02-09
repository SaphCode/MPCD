#include "CircularObstacle.h"
#include "Constants.h"
#include <iostream>

MPCD::CircularObstacle::CircularObstacle(Eigen::Vector2d center, double radius) :
	InteractingBody(std::numeric_limits<double>::infinity(), center, Eigen::Vector2d(0, 0), BodyType::OBSTACLE),
	m_center(center),
	m_radius(radius)
{

}

bool MPCD::CircularObstacle::isInBounds(const Body& o) const
{
	Eigen::Vector2d rel = o.getTorusPosition() - m_center;
	if (rel.stableNorm() <= m_radius) {
		//std::cout << "in Bounds" << std::endl;
		return true;
	}
	return false;
}

Eigen::Vector2d MPCD::CircularObstacle::getOvershoot(const Body& o) const
{
	Eigen::Vector2d oldPos = o.getOldTorusPosition();
	Eigen::Vector2d pos = o.getTorusPosition();
	Eigen::Vector2d relPos = pos - oldPos;

	/*
	int sign;
	if ((oldPos - pos)[1] >= 0) {
		sign = 1;
	}
	else {
		sign = -1;
	}
	*/

	/*

	double k2 = std::pow(k, 2);
	double sqrt = std::sqrt(
		-std::pow(oldPos[1], 2)
		+ 2 * oldPos[1] * m_center[1]
		+ 2 * oldPos[1] * k * oldPos[0]
		- 2 * oldPos[1] * k * m_center[0]
		- std::pow(m_center[1], 2)
		- 2 * m_center[1] * k * oldPos[0]
		+ 2 * m_center[1] * k * m_center[0]
		+ k2 * std::pow(m_radius, 2)
		- k2 * std::pow(oldPos[0], 2)
		+ 2 * k2 * oldPos[0] * m_center[0]
		- k2 * std::pow(m_center[0], 2)
		+ std::pow(m_radius, 2)
	);
	double b
	double denominator = (1 / (k2 + 1)) *
	*/

	double k = (relPos[1] / relPos[0]);
	double a = 1 + std::pow(k, 2);
	/*double b1 = 2 *
		(
			k * (oldPos[1] - m_center[1] - k * oldPos[0])
			- m_center[0]
			);*/
	double b = 2 *
		(
			k * oldPos[1] - k*k * oldPos[0] - m_center[0] - m_center[1] * k
			);
	
	//assert(b == b1);
	double c =
		(m_center.dot(m_center) - m_radius * m_radius
		+ 2* (m_center[1] * oldPos[0] * k - k * oldPos[0] * oldPos[1] - m_center[1] * oldPos[1])
		+ k*k * oldPos[0]*oldPos[0] + oldPos[1] * oldPos[1])
		;
	/*double c1 =
		oldPos[1] * (oldPos[1] - 2 * (k * oldPos[0] + m_center[1]))
		+ oldPos[0] * k * (oldPos[0] * k + 2 * m_center[1])
		+ std::pow(m_center[0], 2) + std::pow(m_center[1], 2) - std::pow(m_radius, 2);*/
	
	//assert(c == c1);
	double x_circle_intersection = 0;
	double y_circle_intersection = 0;
	if (oldPos[0] < pos[0]) {
		x_circle_intersection = (-b - std::sqrt(b * b - 4 * a * c)) / (2 * a);
	}
	else {
		x_circle_intersection = (-b + std::sqrt(b * b - 4 * a * c)) / (2 * a);
	}
	y_circle_intersection = (x_circle_intersection - oldPos[0]) * k + oldPos[1];
	assert((x_circle_intersection != 0) && (x_circle_intersection >= MPCD::Constants::x_0) && (x_circle_intersection <= MPCD::Constants::x_max));
	assert((y_circle_intersection >= MPCD::Constants::y_0) && (y_circle_intersection <= MPCD::Constants::y_max));

	Eigen::Vector2d overshoot = pos - Eigen::Vector2d(x_circle_intersection, y_circle_intersection);
	
	//std::cout << "Circle x: " << m_center[0] << ", y: " << m_center[1] << std::endl;
	return overshoot;
}

void MPCD::CircularObstacle::interact(InteractingBody& b)
{
	// do nothing
}

bool MPCD::CircularObstacle::contains(Eigen::Vector2d point) const
{
	Eigen::Vector2d rel = point - m_center;
	if (rel.stableNorm() <= m_radius) return true;
	return false;
}

bool MPCD::CircularObstacle::occupies(std::pair<int, int> index, double cell_dim) const
{
	int fine_grain = 10; // 10 points per outer line, which means 40 per cell

	Eigen::Vector2d p0(index.second * cell_dim, index.first * cell_dim);
	Eigen::Vector2d p1(((double)index.second + 1) * cell_dim, index.first * cell_dim);
	Eigen::Vector2d p2(((double)index.second + 1) * cell_dim, ((double)index.first + 1) * cell_dim);
	Eigen::Vector2d p3(index.second * cell_dim, ((double)index.first + 1) * cell_dim);

	std::vector<Eigen::Vector2d> keyPoints;

	for (int i = 0; i <= fine_grain; i++) {
		double offset = i * cell_dim / 10 ;
		Eigen::Vector2d p01_offset = p0 + Eigen::Vector2d(offset, 0);
		Eigen::Vector2d p12_offset = p1 + Eigen::Vector2d(0, offset);
		Eigen::Vector2d p23_offset = p2 + Eigen::Vector2d(-offset, 0);
		Eigen::Vector2d p30_offset = p3 + Eigen::Vector2d(0, -offset);

		keyPoints.push_back(p01_offset);
		keyPoints.push_back(p12_offset);
		keyPoints.push_back(p23_offset);
		keyPoints.push_back(p30_offset);
	}

	for (auto& keyPoint : keyPoints) {
		if (contains(keyPoint)) {
			return true;
		}
	}

	return false;
}
