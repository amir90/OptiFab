#pragma once
#include "NSGA3/problem_base.h"
#include "HelperTop.h"
#include "FEMSimulator.h"
#include <MMASolver.h>
#include <cmath>
#include <vector>
#include <iostream>

class ComplianceProblem : public BProblem
{
public:
	ComplianceProblem(std::size_t M, std::size_t k, const std::string &name, std::vector<designVariable> x, std::vector<std::vector<double>> nodes, std::set < std::pair<int, Eigen::Vector3f>, comp> forceVertexIdSet,
		int dx, int numOfVoxelsX, int numOfVoxelsY, int numOfVoxelsZ, double V_allowed_coeff, std::set<int> constraintVertexIdSet);


	virtual std::size_t num_variables() const { return M_ + k_ - 1; }
	virtual std::size_t num_objectives() const { return M_; }

	 bool Evaluate(CIndividual *indv);

	 //add x
	 //add v
	 double r0;
	 int dx, numOfVoxelsX, numOfVoxelsY, numOfVoxelsZ;
	 double V_allowed_coeff;
	 std::vector<std::vector<double>> nodes;
	 std::vector<designVariable> x;
	 std::set < std::pair<int, Eigen::Vector3f>, comp> forceVertexIdSet;
	 std::set<int> constraintVertexIdSet;

	 double change;
	 

protected:
	std::size_t M_; // number of objectives
	std::size_t k_; // number of variables in g(xM)
};
