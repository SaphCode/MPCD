#pragma once

#ifndef COMPARE_H
#define COMPARE_H

#include <Eigen/Dense>

bool areVectorsEqual(Eigen::Vector2d v1, Eigen::Vector2d v2);
bool areVectorsEqual(Eigen::Vector2i v1, Eigen::Vector2i v2);

#endif