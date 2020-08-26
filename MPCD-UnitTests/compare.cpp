#include "pch.h"
#include "Compare.h"
#include <Eigen/Dense>

bool areVectorsEqual(Eigen::Vector2d v1, Eigen::Vector2d v2) {
	if ((v1(0) == (double) v2(0)) && ((double) v1(1) == (double) v2(1))) {
		return true;
	}
	else {
		return false;
	}
}

bool areVectorsEqual(Eigen::Vector2i v1, Eigen::Vector2i v2) {
	if ((v1(0) == (double)v2(0)) && ((double)v1(1) == (double)v2(1))) {
		return true;
	}
	else {
		return false;
	}
}