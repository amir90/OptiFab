#pragma once
#include "FEMSimulator.h"
#include<iostream>
#include<fstream>
#include<set>

Eigen::MatrixXd B_mat(float x, float y, float z) {
	//x,y,z - location in local coordinate systems

	Eigen::MatrixXd B1(6, 3), B2(6, 3), B3(6, 3), B4(6, 3), B5(6, 3), B6(6, 3), B7(6, 3), B8(6, 3);
	Eigen::MatrixXd B(6, 24);

	B1 << -(1 - y)*(1 - z), 0, 0,
		0, -(1 - x)*(1 - z), 0,
		0, 0, -(1 - x)*(1 - y),
		0, -(1 - x)*(1 - y), -(1 - x)*(1 - z),
		-(1 - x)*(1 - y), 0, -(1 - y)*(1 - z),
		-(1 - x)*(1 - z), -(1 - y)*(1 - z), 0;

	//std::cout << B1.rows() << " " << B1.cols() << std::endl;

	B2 << (1 - y)*(1 - z), 0, 0,
		0, -(1 + x)*(1 - z), 0,
		0, 0, -(1 + x)*(1 - y),
		0, -(1 + x)*(1 - y), -(1 + x)*(1 - z),
		-(1 + x)*(1 - y), 0, (1 - y)*(1 - z),
		-(1 + x)*(1 - z), (1 - y)*(1 - z), 0;

	B3 << -(1 + y)*(1 - z), 0, 0,
		0, (1 - x)*(1 - z), 0,
		0, 0, -(1 - x)*(1 + y),
		0, -(1 - x)*(1 + y), (1 - x)*(1 - z),
		-(1 - x)*(1 + y), 0, -(1 + y)*(1 - z),
		(1 - x)*(1 - z), -(1 + y)*(1 - z), 0;

	B4 << (1 + y)*(1 - z), 0, 0,
		0, (1 + x)*(1 - z), 0,
		0, 0, -(1 + x)*(1 + y),
		0, -(1 + x)*(1 + y), (1 + x)*(1 - z),
		-(1 + x)*(1 + y), 0, (1 + y)*(1 - z),
		(1 + x)*(1 - z), (1 + y)*(1 - z), 0;

	B5 << -(1 - y)*(1 + z), 0, 0,
		0, -(1 - x)*(1 + z), 0,
		0, 0, (1 - x)*(1 - y),
		0, (1 - x)*(1 - y), -(1 - x)*(1 + z),
		(1 - x)*(1 - y), 0, -(1 - y)*(1 + z),
		-(1 - x)*(1 + z), -(1 - y)*(1 + z), 0;

	B6 << (1 - y)*(1 + z), 0, 0,
		0, -(1 + x)*(1 + z), 0,
		0, 0, (1 + x)*(1 - y),
		0, (1 + x)*(1 - y), -(1 + x)*(1 + z),
		(1 + x)*(1 - y), 0, (1 - y)*(1 + z),
		-(1 + x)*(1 + z), (1 - y)*(1 + z), 0;

	B7 << -(1 + y)*(1 + z), 0, 0,
		0, (1 - x)*(1 + z), 0,
		0, 0, (1 - x)*(1 + y),
		0, (1 - x)*(1 + y), (1 - x)*(1 + z),
		(1 - x)*(1 + y), 0, -(1 + y)*(1 + z),
		(1 - x)*(1 + z), -(1 + y)*(1 + z), 0;

	B8 << (1 + y)*(1 + z), 0, 0,
		0, (1 + x)*(1 + z), 0,
		0, 0, (1 + x)*(1 + y),
		0, (1 + x)*(1 + y), (1 + x)*(1 + z),
		(1 + x)*(1 + y), 0, (1 + y)*(1 + z),
		(1 + x)*(1 + z), (1 + y)*(1 + z), 0;

	B.block<6, 3>(0, 0) = B1;				B.block<6, 3>(0, 12) = B5;
	B.block<6, 3>(0, 3) = B2;				B.block<6, 3>(0, 15) = B6;
	B.block<6, 3>(0, 6) = B3;				B.block<6, 3>(0, 18) = B7;
	B.block<6, 3>(0, 9) = B4;				B.block<6, 3>(0, 21) = B8;

	// std::cout << B.rows() << " " << B.cols() << std::endl;

	B = B / 8;



	return B;

}

Eigen::MatrixXd calc_N(float x, float y, float z) {

	Eigen::MatrixXd N1(3, 3), N2(3, 3), N3(3, 3), N4(3, 3), N5(3, 3), N6(3, 3), N7(3, 3), N8(3, 3);
	Eigen::MatrixXd N(3, 24);

	float n1 = (1 - x)*(1 - y)*(1 - z);
	N1 << n1, 0, 0,
		0, n1, 0,
		0, 0, n1;

	float n2 = (1 + x)*(1 - y)*(1 - z);
	N2 << n2, 0, 0,
		0, n2, 0,
		0, 0, n2;

	float n3 = (1 - x)*(1 - y)*(1 - z);
	N3 << n3, 0, 0,
		0, n3, 0,
		0, 0, n3;

	float n4 = (1 + x)*(1 + y)*(1 - z);
	N4 << n4, 0, 0,
		0, n4, 0,
		0, 0, n4;

	float n5 = (1 - x)*(1 - y)*(1 + z);
	N5 << n5, 0, 0,
		0, n5, 0,
		0, 0, n5;

	float n6 = (1 + x)*(1 - y)*(1 + z);
	N6 << n6, 0, 0,
		0, n6, 0,
		0, 0, n6;

	float n7 = (1 - x)*(1 + y)*(1 + z);
	N7 << n7, 0, 0,
		0, n7, 0,
		0, 0, n7;

	float n8 = (1 + x)*(1 + y)*(1 + z);
	N8 << n8, 0, 0,
		0, n8, 0,
		0, 0, n8;

	N.block<3, 3>(0, 0) = N1;				N.block<3, 3>(0, 12) = N5;
	N.block<3, 3>(0, 3) = N2;				N.block<3, 3>(0, 15) = N6;
	N.block<3, 3>(0, 6) = N3;				N.block<3, 3>(0, 18) = N7;
	N.block<3, 3>(0, 9) = N4;;				N.block<3, 3>(0, 21) = N8;

	// std::cout << B.rows() << " " << B.cols() << std::endl;

	N = N / 8;

	return N;

}

Eigen::MatrixXd K_mat() {
	float nu = 0.3;

	double k1 = -(6 * nu - 4) / 9;  double k2 = 1.0 / 12; double k3 = -1.0 / 9;
	double k4 = -(4 * nu - 1) / 12; double k5 = (4 * nu - 1) / 12; double k6 = 1.0 / 18;
	double k7 = 1.0 / 24; double k8 = -1.0 / 12; double k9 = (6 * nu - 5) / 36;
	double k10 = -(4 * nu - 1) / 24; double k11 = -1.0 / 24; double k12 = (4 * nu - 1) / 24;
	double k13 = (3 * nu - 1) / 18; double k14 = (3 * nu - 2) / 18;

	Eigen::MatrixXd k1_b(6, 6);
	k1_b << k1, k2, k2, k3, k5, k5,
		k2, k1, k2, k4, k6, k7,
		k2, k2, k1, k4, k7, k6,
		k3, k4, k4, k1, k8, k8,
		k5, k6, k7, k8, k1, k2,
		k5, k7, k6, k8, k2, k1;

	Eigen::MatrixXd k2_b(6, 6);
	k2_b << k9, k8, k12, k6, k4, k7,
		k8, k9, k12, k5, k3, k5,
		k10, k10, k13, k7, k4, k6,
		k6, k5, k11, k9, k2, k10,
		k4, k3, k5, k2, k9, k12,
		k11, k4, k6, k12, k10, k13;

	Eigen::MatrixXd k3_b(6, 6);
	k3_b << k6, k7, k4, k9, k12, k8,
		k7, k6, k4, k10, k13, k10,
		k5, k5, k3, k8, k12, k9,
		k9, k10, k2, k6, k11, k5,
		k12, k13, k10, k11, k6, k4,
		k2, k12, k9, k4, k5, k3;

	Eigen::MatrixXd k4_b(6, 6);
	k4_b << k14, k11, k11, k13, k10, k10,
		k11, k14, k11, k12, k9, k8,
		k11, k11, k14, k12, k8, k9,
		k13, k12, k12, k14, k7, k7,
		k10, k9, k8, k7, k14, k11,
		k10, k8, k9, k7, k11, k14;

	Eigen::MatrixXd k5_b(6, 6);
	k5_b << k1, k2, k8, k3, k5, k4,
		k2, k1, k8, k4, k6, k11,
		k8, k8, k1, k5, k11, k6,
		k3, k4, k5, k1, k8, k2,
		k5, k6, k11, k8, k1, k8,
		k4, k11, k6, k2, k8, k1;

	Eigen::MatrixXd k6_b(6, 6);
	k6_b << k14, k11, k7, k13, k10, k12,
		k11, k14, k7, k12, k9, k2,
		k7, k7, k14, k10, k2, k9,
		k13, k12, k10, k14, k7, k11,
		k10, k9, k2, k7, k14, k7,
		k12, k2, k9, k11, k7, k14;

	Eigen::MatrixXd k0(24, 24);


	k0.block<6, 6>(0, 0) = k1_b;				k0.block<6, 6>(0, 6) = k2_b;	k0.block<6, 6>(0, 12) = k3_b;				k0.block<6, 6>(0, 18) = k4_b;
	k0.block<6, 6>(6, 0) = k2_b.transpose();	k0.block<6, 6>(6, 6) = k5_b;	k0.block<6, 6>(6, 12) = k6_b;				k0.block<6, 6>(6, 18) = k3_b.transpose();
	k0.block<6, 6>(12, 0) = k3_b.transpose();	k0.block<6, 6>(12, 6) = k6_b;	k0.block<6, 6>(12, 12) = k5_b.transpose();	k0.block<6, 6>(12, 18) = k2_b.transpose();
	k0.block<6, 6>(18, 0) = k4_b;				k0.block<6, 6>(18, 6) = k3_b;	k0.block<6, 6>(18, 12) = k2_b;				k0.block<6, 6>(18, 18) = k1_b.transpose();


	k0 = k0 / ((nu + 1)*(1 - 2 * nu));

	Eigen::MatrixXd permutation(24, 24);
	permutation <<
		1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0;


	/*
	std::cout << "freeDOF" << std::endl;
	for (int i = 0; i < freeDOF.size() ; i++) {
		std::cout << freeDOF[i] << std::endl;
	}
	tested and works
	*/

	//	std::cout << permutation << std::endl;

	k0 = permutation * k0*permutation;

	return k0;

}

Eigen::MatrixXd M_mat() {

	Eigen::MatrixXd m(24,24);

	m << 8, 0, 0, 4, 0, 0, 4, 0, 0, 2, 0, 0, 4, 0, 0, 2, 0, 0, 2, 0, 0, 1, 0, 0,
		 0, 8, 0, 0, 4, 0, 0, 4, 0, 0, 2, 0, 0, 4, 0, 0, 2, 0, 0, 2, 0, 0, 1, 0,
		 0, 0, 8, 0, 0, 4, 0, 0, 4, 0, 0, 2, 0, 0, 4, 0, 0, 2, 0, 0, 2, 0, 0, 1,
		 4, 0, 0, 8, 0, 0, 2, 0, 0, 4, 0, 0, 2, 0, 0, 4, 0, 0, 1, 0, 0, 2, 0, 0,
		 0, 4, 0, 0, 8, 0, 0, 2, 0, 0, 4, 0, 0, 2, 0, 0, 4, 0, 0, 1, 0, 0, 2, 0,
		 0, 0, 4, 0, 0, 8, 0, 0, 2, 0, 0, 4, 0, 0, 2, 0, 0, 4, 0, 0, 1, 0, 0, 2,
		 4, 0, 0, 2, 0, 0, 8, 0, 0, 4, 0, 0, 2, 0, 0, 1, 0, 0, 4, 0, 0, 2, 0, 0,
		 0, 4, 0, 0, 2, 0, 0, 8, 0, 0, 4, 0, 0, 2, 0, 0, 1, 0, 0, 4, 0, 0, 2, 0,
		 0, 0, 4, 0, 0, 2, 0, 0, 8, 0, 0, 4, 0, 0, 2, 0, 0, 1, 0, 0, 4, 0, 0, 2,
		 2, 0, 0, 4, 0, 0, 4, 0, 0, 8, 0, 0, 1, 0, 0, 2, 0, 0, 2, 0, 0, 4, 0, 0,
		 0, 2, 0, 0, 4, 0, 0, 4, 0, 0, 8, 0, 0, 1, 0, 0, 2, 0, 0, 2, 0, 0, 4, 0,
		 0, 0, 2, 0, 0, 4, 0, 0, 4, 0, 0, 8, 0, 0, 1, 0, 0, 2, 0, 0, 2, 0, 0, 4,
		 4, 0, 0, 2, 0, 0, 2, 0, 0, 1, 0, 0, 8, 0, 0, 4, 0, 0, 4, 0, 0, 2, 0, 0,
		 0, 4, 0, 0, 2, 0, 0, 2, 0, 0, 1, 0, 0, 8, 0, 0, 4, 0, 0, 4, 0, 0, 2, 0,
		 0, 0, 4, 0, 0, 2, 0, 0, 2, 0, 0, 1, 0, 0, 8, 0, 0, 4, 0, 0, 4, 0, 0, 2,
		 2, 0, 0, 4, 0, 0, 1, 0, 0, 2, 0, 0, 4, 0, 0, 8, 0, 0, 2, 0, 0, 4, 0, 0,
		 0, 2, 0, 0, 4, 0, 0, 1, 0, 0, 2, 0, 0, 4, 0, 0, 8, 0, 0, 2, 0, 0, 4, 0,
		 0, 0, 2, 0, 0, 4, 0, 0, 1, 0, 0, 2, 0, 0, 4, 0, 0, 8, 0, 0, 2, 0, 0, 4,
		 2, 0, 0, 1, 0, 0, 4, 0, 0, 2, 0, 0, 4, 0, 0, 2, 0, 0, 8, 0, 0, 4, 0, 0,
		 0, 2, 0, 0, 1, 0, 0, 4, 0, 0, 2, 0, 0, 4, 0, 0, 2, 0, 0, 8, 0, 0, 4, 0,
		 0, 0, 2, 0, 0, 1, 0, 0, 4, 0, 0, 2, 0, 0, 4, 0, 0, 2, 0, 0, 8, 0, 0, 4,
		 1, 0, 0, 2, 0, 0, 2, 0, 0, 4, 0, 0, 2, 0, 0, 4, 0, 0, 4, 0, 0, 8, 0, 0,
		 0, 1, 0, 0, 2, 0, 0, 2, 0, 0, 4, 0, 0, 2, 0, 0, 4, 0, 0, 4, 0, 0, 8, 0,
		 0, 0, 1, 0, 0, 2, 0, 0, 2, 0, 0, 4, 0, 0, 2, 0, 0, 4, 0, 0, 4, 0, 0, 8;

	return m/27.0;

}



Eigen::SparseMatrix<double> build_global_stiffness_matrix(std::vector<designVariable> const &x, int numNodes, Eigen::MatrixXd KE, int power, std::vector<int> freeDOF, int numOfFree) {

	//KE - element stiffness matrix;
	// x - designVariable Array
	// numNodes - number of free nodes
	// freeDOF - map from total node number to free node number
	//return a sparse global stiffness matrix

	Eigen::SparseMatrix<double> spMat;

	int m = (numOfFree) * 3;  // number of unknowns (=Total number of free DOFs)
	// Assembly:
	std::vector<T> coefficients;            // list of non-zeros coefficients
	std::vector<int> order = { 0,1,2,4,3,6,5,7 };
	double factor;
	for (int l = 0; l < x.size(); l++) {
		//std::cout << "element: " << l << std::endl;
		factor = std::pow(x[l].rho, power);
		for (int i = 0; i < 8; i++) {
			for (int j = 0; j < 8; j++) {
				if ((freeDOF[x[l].nodeIdx[i]] == -1) || (freeDOF[x[l].nodeIdx[j]] == -1)) {
					continue;
				}
				//	std::cout << "node: " << x[l].nodeIdx[i] << " , " << x[l].nodeIdx[j] << std::endl;

				int row1 = freeDOF[x[l].nodeIdx[i]] * 3; int row2 = freeDOF[x[l].nodeIdx[i]] * 3 + 1; int row3 = freeDOF[x[l].nodeIdx[i]] * 3 + 2;
				int col1 = freeDOF[x[l].nodeIdx[j]] * 3; int col2 = freeDOF[x[l].nodeIdx[j]] * 3 + 1; int col3 = freeDOF[x[l].nodeIdx[j]] * 3 + 2;
				double val11 = KE(order[i] * 3, order[j] * 3)*factor; double val12 = KE(order[i] * 3, order[j] * 3 + 1)*factor; double val13 = KE(order[i] * 3, order[j] * 3 + 2)*factor;
				double val21 = KE(order[i] * 3 + 1, order[j] * 3)*factor; double val22 = KE(order[i] * 3 + 1, order[j] * 3 + 1)*factor; double val23 = KE(order[i] * 3 + 1, order[j] * 3 + 2)*factor;
				double val31 = KE(order[i] * 3 + 2, order[j] * 3)*factor; double val32 = KE(order[i] * 3 + 2, order[j] * 3 + 1)*factor; double val33 = KE(order[i] * 3 + 2, order[j] * 3 + 2)*factor;
				//	std::cout << "val11: " << order[i] * 3 << order[j] * 3 << std::endl;
				//	std::cout << "val12: " << order[i] *3 << order[j] * 3+1 << std::endl;
				//	std::cout << "val13: " << order[i] * 3 << order[j] * 3+2 << std::endl;
				//	std::cout << "val21: " << order[i] * 3+1 << order[j] * 3 << std::endl;
				//	std::cout << "val22: " << order[i] * 3 +1<< order[j] * 3+1 << std::endl;
				//	std::cout << "val23: " << order[i] * 3 +1<< order[j] * 3+2 << std::endl;
				//	std::cout << "val31: " << order[i] * 3 +2<< order[j] * 3 << std::endl;
				//	std::cout << "val32: " << order[i] * 3 +2<< order[j] * 3+1 << std::endl;
				//	std::cout << "val33: " << order[i] * 3 +2<< order[j] * 3+2 << std::endl;

				coefficients.push_back(T(row1, col1, val11)); coefficients.push_back(T(row1, col2, val12)); coefficients.push_back(T(row1, col3, val13));
				coefficients.push_back(T(row2, col1, val21)); coefficients.push_back(T(row2, col2, val22)); coefficients.push_back(T(row2, col3, val23));
				coefficients.push_back(T(row3, col1, val31)); coefficients.push_back(T(row3, col2, val32)); coefficients.push_back(T(row3, col3, val33));
			}
		}
	}

	SpMat A(m, m);
	A.setFromTriplets(coefficients.begin(), coefficients.end());

	//	std::cout << "made sparse global matrix!" <<A<< std::endl;
	return A;

}

Eigen::SparseMatrix<double> build_global_mass_matrix(std::vector<designVariable> const &x, int numNodes, Eigen::MatrixXd ME, int power, std::vector<int> freeDOF, int numOfFree) {

	//ME - element mass matrix;
	// x - designVariable Array
	// numNodes - number of free nodes
	// freeDOF - map from total node number to free node number
	//return a sparse global stiffness matrix

	Eigen::SparseMatrix<double> spMat;

	std::cout << "got here - build global stiffness" << std::endl;
	int m = (numOfFree) * 3;  // number of unknowns (=Total number of free DOFs)
	// Assembly:
	std::vector<T> coefficients;            // list of non-zeros coefficients
	std::vector<int> order = { 0,1,2,4,3,6,5,7 };
	double factor;
	for (int l = 0; l < x.size(); l++) {
		//std::cout << "element: " << l << std::endl;
		factor = std::pow(x[l].rho, power);
		for (int i = 0; i < 8; i++) {
			for (int j = 0; j < 8; j++) {
				if ((freeDOF[x[l].nodeIdx[i]] == -1) || (freeDOF[x[l].nodeIdx[j]] == -1)) {
					continue;
				}
				//	std::cout << "node: " << x[l].nodeIdx[i] << " , " << x[l].nodeIdx[j] << std::endl;

				int row1 = freeDOF[x[l].nodeIdx[i]] * 3; int row2 = freeDOF[x[l].nodeIdx[i]] * 3 + 1; int row3 = freeDOF[x[l].nodeIdx[i]] * 3 + 2;
				int col1 = freeDOF[x[l].nodeIdx[j]] * 3; int col2 = freeDOF[x[l].nodeIdx[j]] * 3 + 1; int col3 = freeDOF[x[l].nodeIdx[j]] * 3 + 2;
				double val11 = ME(order[i] * 3, order[j] * 3)*factor; double val12 = ME(order[i] * 3, order[j] * 3 + 1)*factor; double val13 = ME(order[i] * 3, order[j] * 3 + 2)*factor;
				double val21 = ME(order[i] * 3 + 1, order[j] * 3)*factor; double val22 = ME(order[i] * 3 + 1, order[j] * 3 + 1)*factor; double val23 = ME(order[i] * 3 + 1, order[j] * 3 + 2)*factor;
				double val31 = ME(order[i] * 3 + 2, order[j] * 3)*factor; double val32 = ME(order[i] * 3 + 2, order[j] * 3 + 1)*factor; double val33 = ME(order[i] * 3 + 2, order[j] * 3 + 2)*factor;
				//	std::cout << "val11: " << order[i] * 3 << order[j] * 3 << std::endl;
				//	std::cout << "val12: " << order[i] *3 << order[j] * 3+1 << std::endl;
				//	std::cout << "val13: " << order[i] * 3 << order[j] * 3+2 << std::endl;
				//	std::cout << "val21: " << order[i] * 3+1 << order[j] * 3 << std::endl;
				//	std::cout << "val22: " << order[i] * 3 +1<< order[j] * 3+1 << std::endl;
				//	std::cout << "val23: " << order[i] * 3 +1<< order[j] * 3+2 << std::endl;
				//	std::cout << "val31: " << order[i] * 3 +2<< order[j] * 3 << std::endl;
				//	std::cout << "val32: " << order[i] * 3 +2<< order[j] * 3+1 << std::endl;
				//	std::cout << "val33: " << order[i] * 3 +2<< order[j] * 3+2 << std::endl;

				coefficients.push_back(T(row1, col1, val11)); coefficients.push_back(T(row1, col2, val12)); coefficients.push_back(T(row1, col3, val13));
				coefficients.push_back(T(row2, col1, val21)); coefficients.push_back(T(row2, col2, val22)); coefficients.push_back(T(row2, col3, val23));
				coefficients.push_back(T(row3, col1, val31)); coefficients.push_back(T(row3, col2, val32)); coefficients.push_back(T(row3, col3, val33));
			}
		}
	}

	SpMat A(m, m);
	A.setFromTriplets(coefficients.begin(), coefficients.end());

	//	std::cout << "made sparse global matrix!" <<A<< std::endl;
	return A;


}

Eigen::SparseMatrix<double> build_global_stiffness_matrix_full(std::vector<designVariable> const &x, Eigen::MatrixXd KE, int power, std::vector<std::vector<double>> &nodes) {

	//KE - element stiffness matrix;
	// x - designVariable Array
	// numNodes - number of free nodes
	// freeDOF - map from total node number to free node number
	//return a sparse global stiffness matrix

	Eigen::SparseMatrix<double> spMat;

	std::cout << "got here - build global stiffness" << std::endl;
	int m = (nodes.size()) * 3;  // number of unknowns (=Total number of free DOFs)
	// Assembly:
	std::vector<T> coefficients;            // list of non-zeros coefficients
	std::vector<int> order = { 0,1,2,4,3,6,5,7 };
	double factor;
	for (int l = 0; l < x.size(); l++) {
		//std::cout << "element: " << l << std::endl;
		factor =  std::pow(x[l].rho, power);
		for (int i = 0; i < 8; i++) {
			for (int j = 0; j < 8; j++) {

				//	std::cout << "node: " << x[l].nodeIdx[i] << " , " << x[l].nodeIdx[j] << std::endl;

				int row1 = x[l].nodeIdx[i] * 3; int row2 = x[l].nodeIdx[i] * 3 + 1; int row3 = x[l].nodeIdx[i] * 3 + 2;
				int col1 = x[l].nodeIdx[j] * 3; int col2 = x[l].nodeIdx[j] * 3 + 1; int col3 = x[l].nodeIdx[j] * 3 + 2;
				double val11 = KE(order[i] * 3, order[j] * 3)*factor; double val12 = KE(order[i] * 3, order[j] * 3 + 1)*factor; double val13 = KE(order[i] * 3, order[j] * 3 + 2)*factor;
				double val21 = KE(order[i] * 3 + 1, order[j] * 3)*factor; double val22 = KE(order[i] * 3 + 1, order[j] * 3 + 1)*factor; double val23 = KE(order[i] * 3 + 1, order[j] * 3 + 2)*factor;
				double val31 = KE(order[i] * 3 + 2, order[j] * 3)*factor; double val32 = KE(order[i] * 3 + 2, order[j] * 3 + 1)*factor; double val33 = KE(order[i] * 3 + 2, order[j] * 3 + 2)*factor;
				//	std::cout << "val11: " << order[i] * 3 << order[j] * 3 << std::endl;
				//	std::cout << "val12: " << order[i] *3 << order[j] * 3+1 << std::endl;
				//	std::cout << "val13: " << order[i] * 3 << order[j] * 3+2 << std::endl;
				//	std::cout << "val21: " << order[i] * 3+1 << order[j] * 3 << std::endl;
				//	std::cout << "val22: " << order[i] * 3 +1<< order[j] * 3+1 << std::endl;
				//	std::cout << "val23: " << order[i] * 3 +1<< order[j] * 3+2 << std::endl;
				//	std::cout << "val31: " << order[i] * 3 +2<< order[j] * 3 << std::endl;
				//	std::cout << "val32: " << order[i] * 3 +2<< order[j] * 3+1 << std::endl;
				//	std::cout << "val33: " << order[i] * 3 +2<< order[j] * 3+2 << std::endl;

				coefficients.push_back(T(row1, col1, val11)); coefficients.push_back(T(row1, col2, val12)); coefficients.push_back(T(row1, col3, val13));
				coefficients.push_back(T(row2, col1, val21)); coefficients.push_back(T(row2, col2, val22)); coefficients.push_back(T(row2, col3, val23));
				coefficients.push_back(T(row3, col1, val31)); coefficients.push_back(T(row3, col2, val32)); coefficients.push_back(T(row3, col3, val33));

			}
		}
	}
	SpMat A(m, m);
	A.setFromTriplets(coefficients.begin(), coefficients.end());

	//	std::cout << "made sparse global matrix!" <<A<< std::endl;
	return A;
}

Eigen::VectorXd do_FEM(SpMat &A, Eigen::VectorXd &f) {
	//recieves a global stiffness matrix, the free DOFs, and the forces, and outputs the displacement for the free DOFs


	//std::cout << A.rows() <<" "<< A.cols() << std::endl;
	//Eigen::BiCGSTAB<SpMat> solver;
	//solver.setMaxIterations(61);
	//solver.compute(A);
	//Eigen::VectorXd u = solver.solve(f);
	//std::cout << "#iterations:     " << solver.iterations() << std::endl;
	//std::cout << "estimated error: " << solver.error() << std::endl;
	/* ... update b ... */
	//u = solver.solve(f); // solve again
//	std::cout << Eigen::nbThreads() << std::endl;
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> cg;
	cg.compute(A);
//	Eigen::SimplicialCholesky<SpMat> chol(A);
	//std::cout << A.rows() << " " << A.cols() << " " << f.rows() << " " << f.cols() << std::endl;
	cg.setTolerance(0.01);
//	cg.setMaxIterations(200);
  //  Eigen::BenchTimer t;
	//t.reset(); t.start();
	//Eigen::initParallel();
	std::cout << "start solving for displacements" << std::endl;
	std::cout << "number of cores: " << Eigen::nbThreads()<< std::endl;
	Eigen::VectorXd u = cg.solve(f);//chol.solve(f); //no performance gained with multithreading.. why?
	//t.stop();
	//std::cout << "Real time: " << t.value(1) << " error: "<< cg.error()<< std::endl; // 0=CPU_TIMER, 1=REAL_TIMER
	std::cout << "done solving for displacements" << std::endl;
	return u;
}

Eigen::VectorXd doFEM(std::vector<std::vector<double>> const &nodes, std::vector<designVariable>  const &x, std::set < std::pair<int, Eigen::Vector3f>, comp> &forceVertexIdSet, std::set<int> &constraintVertexIdSet, float dx, double E) {

	Eigen::MatrixXd KE(24, 24); //assumes x[0] value of node is 1 at the start	(?)
	KE = K_mat();
	int power = 3;
	std::vector<int> freeDOF;
	freeDOF.resize(nodes.size());
	int ind = 0; //number of free nodes
	for (int i = 0; i < freeDOF.size(); i++) {

		if (constraintVertexIdSet.find(i) == constraintVertexIdSet.end()) {
			//if (constraints.find(i) == constraints.end()) {

			freeDOF[i] = ind;
			ind++;
		}
		else {

			freeDOF[i] = -1;
		}

	}

	int numOfFixed = nodes.size() - ind; //number of fixed nodes

	Eigen::VectorXd f(3 * ind), u(nodes.size() * 3);

	for (int j = 0; j < 3 * ind; j++) {
		f(j) = 0;
	}

	for (auto j = forceVertexIdSet.begin(); j != forceVertexIdSet.end(); j++) {
		int tempForceid = freeDOF[j->first] * 3;
		f[tempForceid + 0] = j->second[0];// force[temp_f].second[0];
		f[tempForceid + 1] = j->second[1];//force[temp_f].second[1];
		f[tempForceid + 2] = j->second[2]; //force[temp_f].second[2];//

	}

	std::cout << "start building global stiffness matrix" << std::endl;

	//local Ba 
	Eigen::MatrixXd  local_Ba(6, 24);

	SpMat K_global = build_global_stiffness_matrix(x, nodes.size(), KE, power, freeDOF, ind);
	K_global = K_global*E*dx;

	std::cout << "done building global stiffness matrix" << std::endl;

//	std::cout << f.nonZeros() << std::endl;

	Eigen::VectorXd u_free = do_FEM(K_global, f);

	int u_ind = 0;
	for (int k = 0; k < nodes.size(); k++) {

		if (freeDOF[k] != -1) {
			u[k * 3] = u_free[u_ind++];
			u[k * 3 + 1] = u_free[u_ind++];
			u[k * 3 + 2] = u_free[u_ind++];
		}
		else {
			u[k * 3] = 0;
			u[k * 3 + 1] = 0;
			u[k * 3 + 2] = 0;
		}
	}

	return u;

}


Eigen::VectorXd doFEM(std::vector<std::vector<double>> const &nodes, std::vector<designVariable>  const &x, std::set < std::pair<int, Eigen::Vector3f>, comp> &forceVertexIdSet, std::set<int> &constraintVertexIdSet, float dx, double E, SpMat &K_global,Eigen::VectorXd &f, std::vector<int> &freeDOF) {

	Eigen::MatrixXd KE(24, 24); //assumes x[0] value of node is 1 at the start	(?)
	KE = K_mat();
	int power = 3;


	std::cout << "start building global stiffness matrix" << std::endl;

	K_global = K_global * E*dx;

	std::cout << "done building global stiffness matrix" << std::endl;

	Eigen::VectorXd u(nodes.size() * 3);

	Eigen::VectorXd u_free = do_FEM(K_global, f);

	int u_ind = 0;
	for (int k = 0; k < nodes.size(); k++) {

		if (freeDOF[k] != -1) {
			u[k * 3] = u_free[u_ind++];
			u[k * 3 + 1] = u_free[u_ind++];
			u[k * 3 + 2] = u_free[u_ind++];
		}
		else {
			u[k * 3] = 0;
			u[k * 3 + 1] = 0;
			u[k * 3 + 2] = 0;
		}
	}

	return u;

}

void getBCforVoxelTest(std::vector<designVariable> &x, std::vector<std::vector<double>> &nodes, int numOfVoxelsX, int numOfVoxelsY, int numOfVoxelsZ, std::set < std::pair<int, Eigen::Vector3f>, comp> &forceVertexIdSet, std::set<int> &constraintVertexIdSet ) {

	int planeInd = 1;
	bool planeVariable = true;
	std::pair<int, Eigen::Vector3f> tempForce;
	int numOfVoxelsinPlane = 0;
	int startVoxel = -1;
	int nodecounter = 0;
	std::vector<int> vertexList;
	if (planeVariable) {

		if (planeInd == 0) {
			startVoxel = 0;
			numOfVoxelsinPlane = numOfVoxelsY * numOfVoxelsZ;
			nodecounter = 4;
			vertexList = { 0,2,3,5 };

		}
		if (planeInd == 1) {
			startVoxel = numOfVoxelsX - 1;
			numOfVoxelsinPlane = numOfVoxelsY * numOfVoxelsZ;
			nodecounter = 4;
			vertexList = { 1,4,6,7 };
		}
	}
	else {
		//	std::cout << "testing select: " << fid << " " << selectedV << " " << F(fid, selectedV) << " " << F.rows() << " " << V.rows() << std::endl;
		//	tempForce.first = F(fid, selectedV);
	//		tempForce.second = Eigen::Vector3f(forceX, forceY, forceZ);
		numOfVoxelsinPlane = 1;
		nodecounter = 1;
	//	std::cout << "testing end:" << std::endl;
	}

	int currVoxel = startVoxel;

	while (numOfVoxelsinPlane) { //hack for selecting plane or individual vertices

	//	std::cout << "current voxel is: " << currVoxel << std::endl;

		while (nodecounter > 0) {

			if (planeVariable) {
				tempForce.first = x[currVoxel].nodeIdx[vertexList[nodecounter - 1]];
				//		tempForce.second = Eigen::Vector3f(forceX, forceY, forceZ);
			}
			nodecounter--;

			if (true) { //constraints are on - delete forces
			//	std::cout << "got here" << std::endl;
				if (constraintVertexIdSet.find(tempForce.first) == constraintVertexIdSet.end()) {
					constraintVertexIdSet.insert(tempForce.first);
					forceVertexIdSet.erase(tempForce);
			//		std::cout << "got here1" << std::endl;
					std::pair<int, Eigen::Vector3f> tempForce;
			//		tempForce.first = x[0].nodeIdx[0];
			//		tempForce.second = Eigen::Vector3f(0, -10, 0);
			//		forceVertexIdSet.insert(tempForce);
			//		tempForce.first = x[0].nodeIdx[3];
			//		tempForce.second = Eigen::Vector3f(0, -10, 0);
			//		forceVertexIdSet.insert(tempForce);
				}
				else {
					if (!planeVariable) {
						constraintVertexIdSet.erase(tempForce.first);
					}
				}

			}
			else {
				if (forceVertexIdSet.find(tempForce) == forceVertexIdSet.end()) {
		//			std::cout << "got here2" << std::endl;
					constraintVertexIdSet.erase(tempForce.first);
					forceVertexIdSet.insert(tempForce);
		//			std::cout << "got here3" << std::endl;
				}
				else {
					if (!planeVariable) {
						forceVertexIdSet.erase(tempForce);
					}
				}

			}



		}

		numOfVoxelsinPlane--;
		if (planeVariable) {
			currVoxel = moveY(x, currVoxel, 1);
			if (currVoxel == -1) {
				startVoxel = moveZ(x, startVoxel, 1);
				currVoxel = startVoxel;

			}
		}
		nodecounter = 4;

//		std::cout << "got to end of loop" << std::endl;

	}

	planeInd = 0;
	 planeVariable = true;
	 numOfVoxelsinPlane = 0;
	 startVoxel = -1;
	 nodecounter = 0;
	 vertexList;
	if (planeVariable) {

		if (planeInd == 0) {
			startVoxel = 0;
			numOfVoxelsinPlane = numOfVoxelsY * numOfVoxelsZ;
			nodecounter = 4;
			vertexList = { 0,2,3,5 };

		}
		if (planeInd == 1) {
			startVoxel = numOfVoxelsX - 1;
			numOfVoxelsinPlane = numOfVoxelsY * numOfVoxelsZ;
			nodecounter = 4;
			vertexList = { 1,4,6,7 };
		}
	}
	else {
		//	std::cout << "testing select: " << fid << " " << selectedV << " " << F(fid, selectedV) << " " << F.rows() << " " << V.rows() << std::endl;
		//	tempForce.first = F(fid, selectedV);
	//		tempForce.second = Eigen::Vector3f(forceX, forceY, forceZ);
		numOfVoxelsinPlane = 1;
		nodecounter = 1;
	//	std::cout << "testing end:" << std::endl;
	}

	 currVoxel = startVoxel;

	while (numOfVoxelsinPlane) { //hack for selecting plane or individual vertices

//		std::cout << "current voxel is: " << currVoxel << std::endl;

		while (nodecounter > 0) {

			if (planeVariable) {
				tempForce.first = x[currVoxel].nodeIdx[vertexList[nodecounter - 1]];
				tempForce.second = Eigen::Vector3f(0, -20, 0);
			}
			nodecounter--;

			if (false) { //constraints are on - delete forces
//				std::cout << "got here" << std::endl;
				if (constraintVertexIdSet.find(tempForce.first) == constraintVertexIdSet.end()) {
					constraintVertexIdSet.insert(tempForce.first);
					forceVertexIdSet.erase(tempForce);
//					std::cout << "got here1" << std::endl;
				//	std::pair<int, Eigen::Vector3f> tempForce;
					//		tempForce.first = x[0].nodeIdx[0];
					//		tempForce.second = Eigen::Vector3f(0, -10, 0);
					//		forceVertexIdSet.insert(tempForce);
					//		tempForce.first = x[0].nodeIdx[3];
					//		tempForce.second = Eigen::Vector3f(0, -10, 0);
					//		forceVertexIdSet.insert(tempForce);
				}
				else {
					if (!planeVariable) {
						constraintVertexIdSet.erase(tempForce.first);
					}
				}

			}
			else {
				if (forceVertexIdSet.find(tempForce) == forceVertexIdSet.end()) {
//					std::cout << "got here2" << std::endl;
					constraintVertexIdSet.erase(tempForce.first);
					forceVertexIdSet.insert(tempForce);
//					std::cout << "got here3" << std::endl;
				}
				else {
					if (!planeVariable) {
						forceVertexIdSet.erase(tempForce);
					}
				}

			}

		}

		numOfVoxelsinPlane--;
		if (planeVariable) {
			currVoxel = moveY(x, currVoxel, 1);
			if (currVoxel == -1) {
				startVoxel = moveZ(x, startVoxel, 1);
				currVoxel = startVoxel;

			}
		}
		nodecounter = 4;

//		std::cout << "got to end of loop" << std::endl;

	}


}