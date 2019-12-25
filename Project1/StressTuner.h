#pragma once
#include "FEMSimulator.h" //need to generalize this
#include<iostream>
#include<fstream>
#include<set>
#include<vector>
#include<MMASolver.h>

struct Stress {
	int idx = 0; //idx of design variable;
	double vonMises = 0; //von mises stress
	Eigen::VectorXd Stresses = Eigen::VectorXd::Zero(6);
	Eigen::VectorXd displacements = Eigen::VectorXd::Zero(24);


};

double calc_DSigmaPN_DSigmaVM(std::vector<Stress> const &stresses, int k, int beginInd, int endInd, double p);
double calcVonMises(Eigen::VectorXd const  &stresses);
Eigen::VectorXd calc_DsigmaVM_Dsigma_a(Stress const &stresses);
double  calc_DniS_Drho_e(double rho_e);
double calc_Drho_e_Dxb(int idx, int j, std::vector<designVariable> const &x, double r0, std::vector<std::vector<double>> const &nodes);
double ni_s(double rho);
double error(int n, double* x_new, double* x_old);
std::vector<Stress> calcStresses(double ni, double E, std::vector<std::vector<double>> const nodes, std::vector<designVariable>  const x, Eigen::VectorXd  u, float dx);
std::vector<Stress> calcStresses(double ni, double E, std::vector<std::vector<double>> const nodes, std::vector<designVariable>  const x, Eigen::VectorXd  u, float dx, std::vector<int> targets);
double optimizeStress(std::vector<Stress> Stresses, Eigen::VectorXd u, int C, int power, float dx, double change, std::vector<std::vector<double>> nodes, std::vector<designVariable>  x,
	double* xnew, double* dg, double* g, double* df, double* xmin, double *xmax, double density, double E, double ni, MMASolver* optimizer, int neighbourLayers,float r0);
