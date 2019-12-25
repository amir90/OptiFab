#pragma once
#include "NSGAIIIComp.h"
#include "NSGA3/alg_individual.h"
//#include "ComplianceTuner.h"

extern bool debug_flag; 

ComplianceProblem::ComplianceProblem(size_t M, size_t k, const std::string &name, std::vector<designVariable> x, std::vector<std::vector<double>> nodes, std::set < std::pair<int, Eigen::Vector3f>, comp> forceVertexIdSet,
	int dx, int numOfVoxelsX, int numOfVoxelsY, int numOfVoxelsZ, double V_allowed_coeff, std::set<int> constraintVertexIdSet) :
	BProblem(name),
	M_(M), //number of objectives
	k_(k) // number of variables
{
	lbs_.resize(M_ + k_ - 1, 0.0); // lower bound 0.0
	ubs_.resize(M_ + k_ - 1, 1.0); // upper bound 1.0

	this->dx = dx;
	this->constraintVertexIdSet = constraintVertexIdSet;
	this->nodes = nodes;
	this->forceVertexIdSet = forceVertexIdSet;
	this->numOfVoxelsX = numOfVoxelsX;
	this->numOfVoxelsY = numOfVoxelsY;
	this->numOfVoxelsZ = numOfVoxelsZ;
	this->x = x;
}


bool ComplianceProblem::Evaluate(CIndividual *indv) {

	CIndividual::TDecVec &x = indv->vars(); //vector of variables
	CIndividual::TObjVec &f = indv->objs(); //vector of objectives

	for (int i = 0; i < M_; i++) {

		f[i] = 0;

	}

	for (int i = 0; i < k_; i++) {

		this->x[i].rho = x[i];

	}
	
	if (debug_flag) {
		std::cout << "number of variables: " << x.size() << std::endl;
		std::cout << "number of objectives: " << f.size() << std::endl;

		std::cout << "number of forces: " << forceVertexIdSet.size() << std::endl;
		std::cout << "force: " << forceVertexIdSet.begin()->second << std::endl;

		std::cout << "number of constraints: " << constraintVertexIdSet.size() << std::endl;
	}

	if (x.size() != M_ + k_ - 1) return false; // #variables does not match

	//input vars(=values) to x vector

	//calc_rho(nodes, x_design, -1, dx, numOfVoxelsX, numOfVoxelsY, numOfVoxelsZ, r0);

	//std::cout << "Genetic Algorithm - doing FEM" << std::endl;

	//calculate compliance 

	double sum = 0;
	for (int i = 0; i < k_; i++) {

		sum += x[i];
		this->x[i].value = x[i];
		this->x[i].rho = x[i];
	}

	//for (int i = 0; i < M_; i++) {


		Eigen::VectorXd u = doFEM(nodes, this->x, forceVertexIdSet, constraintVertexIdSet, dx, 1);

		f[0] = 0;
		f[1] = 0;

		for (auto j = forceVertexIdSet.begin(); j != forceVertexIdSet.end(); j++) {
			
			int nodeIdx = j->first;
			f[0] += std::abs( j->second(0)*u[3 * nodeIdx] + j->second(1)*u[3 * nodeIdx + 1] + j->second(2)*u[3 * nodeIdx + 2]);

		}

		if (sum > V_allowed_coeff*numOfVoxelsX*numOfVoxelsY*numOfVoxelsZ) {

			f[1] += (f[0]-1)*((sum-V_allowed_coeff*numOfVoxelsX*numOfVoxelsY*numOfVoxelsZ)/sum);
		}

		if (debug_flag) {
		std::cout << "compliance fitness: " << f[0] << std::endl;
		std::cout << "weight fitness: " << f[1] << std::endl;
			}

//	}

	return true;
	

}

