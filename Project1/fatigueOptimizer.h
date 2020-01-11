#pragma once
#include "HelperTop.h"
#include "FEMSimulator.h"
#include "StressTuner.h"
#include <vector>
#include <iostream>


std::vector<std::vector<Stress>>   ClacStressTimeHistory(std::vector<int> targetElements, std::vector<designVariable> x, double a, double b, int power, std::set < std::pair<int, Eigen::Vector3f>, comp> forceVertexIdSet, double dx, double E, std::vector<std::vector<double>> nodes, std::set<int> constraintVertexIdSet, Eigen::VectorXd u_dot_0, Eigen::VectorXd u0, double dt, int steps, std::vector<std::vector<double>> &sigma_m, std::vector<std::vector<Stress>> &sigma_a, Eigen::MatrixXd& u_test, std::vector<Eigen::VectorXd> &u_max, std::vector<Eigen::VectorXd> &u_min);

double dD_PN_k_dgamma_e(int e_id, std::vector<designVariable> x, std::vector<int> Omega_k, double p, std::vector<double> damage, double r0, std::vector<std::vector<double>> nodes);

double dD_PN_k_dDe(int e_id, double p, std::vector<designVariable> x, std::vector<int> Omega_k, std::vector<double> damage);

double dDe_dsigma_a_i(double Sut, double bf, double Sf, std::vector<double> sigma_a, std::vector<double> sigma_m, int cycle, std::vector<int> n);

//Eigen::RowVectorXd dsigma_a_i_d_sigma_a_i_vec(); exists in stress tuner

Eigen::VectorXd dsigma_a_i_vec_dgamma_e(int e_id, double p, std::vector<designVariable> x, Eigen::VectorXd u_max, Eigen::VectorXd u_min, double E);

double dDe_dsigma_m_i(double Sut, double bf, double Sf, std::vector<double> sigma_a, std::vector<double> sigma_m, int cycle, std::vector<int> n);

inline Eigen::RowVectorXd dsigma_m_i_dsigma_m_i_vec();

Eigen::Vector3d dsigma_m_i_vec_dgamma_e(int e_id, double p, std::vector<designVariable> x, Eigen::VectorXd u_max, Eigen::VectorXd u_min);

Eigen::VectorXd Lambda_max_i( double p, std::vector<designVariable> x, std::vector<int> Omega_k, std::vector<Stress> sigma_a, std::vector<double> damage, double Sut, double bf, double Sf, std::vector<double> sigma_m, double E, int cycle, std::vector<int> n);

//Eigen::SparseMatrixBase<Eigen::MatrixXd> dK_dgamma_e();

Eigen::VectorXd Lambda_min_i(double p, std::vector<designVariable> x, std::vector<int> Omega_k, std::vector<Stress> sigma_a, std::vector<double> damage, double Sut, double bf, double Sf, std::vector<double> sigma_m, double E, int cycle, std::vector<int> n);

std::vector<double> full_dD_PN_k_dgamma_e(std::vector<designVariable> &x, std::vector<std::vector<int>> Omega_k, std::vector<int> n, double p, std::vector<Eigen::VectorXd> u_max, std::vector<Eigen::VectorXd> u_min, std::vector<Stress> &sigma_a, std::vector<double> &damage, double Sut, double bf, double Sf, std::vector<double> &sigma_m, double E, std::vector<std::vector<double>> &nodes, double r0);



