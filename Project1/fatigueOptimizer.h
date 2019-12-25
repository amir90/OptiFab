#pragma once
#include "HelperTop.h"
#include "FEMSimulator.h"
#include "StressTuner.h"
#include <vector>
#include <iostream>


std::vector<std::vector<Stress>>   ClacStressTimeHistory(std::vector<int> targetElements, std::vector<designVariable> x, double a, double b, int power, std::set < std::pair<int, Eigen::Vector3f>, comp> forceVertexIdSet, double dx, double E, std::vector<std::vector<double>> nodes, std::set<int> constraintVertexIdSet, Eigen::VectorXd u_dot_0, Eigen::VectorXd u0, double dt, int steps, std::vector<std::vector<double>> &sigma_m, std::vector<std::vector<double>> &sigma_a, Eigen::VectorXd& u_test);




