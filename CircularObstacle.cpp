#include "CircularObstacle.h"
#include "Constants.h"
#include <iostream>
#include <float.h>

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

	double x_circle_intersection = 0;
	double y_circle_intersection = 0;
	if (std::abs(relPos[0]) < 5e-6) {
		x_circle_intersection = pos[0];
		if (oldPos[1] > pos[1]) {
			y_circle_intersection = std::sqrt(std::pow(m_radius, 2) - std::pow(x_circle_intersection - m_center[0], 2)) + m_center[1];
		}
		else {
			y_circle_intersection = -std::sqrt(std::pow(m_radius, 2) - std::pow(x_circle_intersection - m_center[0], 2)) + m_center[1];
		}
	}
	else {
		double k = (relPos[1] / relPos[0]);
		double a = 1 + std::pow(k, 2);

		double b = 2 *
			(
				k * oldPos[1] - k * k * oldPos[0] - m_center[0] - m_center[1] * k
				);

		double c =
			(m_center.dot(m_center) - m_radius * m_radius
				+ 2 * (m_center[1] * oldPos[0] * k - k * oldPos[0] * oldPos[1] - m_center[1] * oldPos[1])
				+ k * k * oldPos[0] * oldPos[0] + oldPos[1] * oldPos[1])
			;

		double x_circle_intersection = 0;
		double y_circle_intersection = 0;
		double toRoot = b * b - 4 * a * c;
		if (toRoot < 0) {
			std::cout << "Root is smaller 0 for:\n";
			std::cout << "REl pos[0] * 1000: " << relPos[0] * 1000 << std::endl;
			std::cout << "Mass: " << o.getMass() << "\nPos: " << o.getPosition() << "\n";
			std::cout << "Old Pos: " << o.getOldPosition() << "\n";
			std::cout << "Vel: " << m_vel << "\n";
			exit(-1);
		}
		if (oldPos[0] < pos[0]) {
			x_circle_intersection = (-b - std::sqrt(b * b - 4 * a * c)) / (2 * a);
		}
		else {
			x_circle_intersection = (-b + std::sqrt(b * b - 4 * a * c)) / (2 * a);
		}
		y_circle_intersection = (x_circle_intersection - oldPos[0]) * k + oldPos[1];
		assert((x_circle_intersection != 0) && (x_circle_intersection >= MPCD::Constants::x_0) && (x_circle_intersection <= MPCD::Constants::x_max));
		assert((y_circle_intersection >= MPCD::Constants::y_0) && (y_circle_intersection <= MPCD::Constants::y_max));
	}

	Eigen::Vector2d overshoot = pos - Eigen::Vector2d(x_circle_intersection, y_circle_intersection);

	if (_isnan(overshoot[0]) || _isnan(overshoot[1])) {
		exit(-1);
	}
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

bool MPCD::CircularObstacle::occupies(std::pair<int, int> index, Eigen::Vector2d shift, double cell_dim) const
{
	Eigen::Vector2d pll = Eigen::Vector2d(index.second * cell_dim, index.first * cell_dim) + shift;
	Eigen::Vector2d plr = Eigen::Vector2d(index.second * cell_dim, ((double)index.first + 1) * cell_dim) + shift;
	Eigen::Vector2d pul = Eigen::Vector2d(((double)index.second + 1) * cell_dim, index.first * cell_dim) + shift;
	Eigen::Vector2d pur = Eigen::Vector2d(((double)index.second + 1) * cell_dim, ((double)index.first + 1) * cell_dim) + shift;
	return contains(pll) || contains(plr) || contains(pul) || contains(pur);
}
