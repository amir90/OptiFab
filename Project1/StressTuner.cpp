#pragma once
#include "StressTuner.h"


bool compareByLength(const Stress &a, const Stress & b)
{
	return a.vonMises < b.vonMises;
}


double calc_DSigmaPN_DSigmaVM(std::vector<Stress> const &stresses, int k, int beginInd, int endInd, double p) {
	//stresses - stress information vector
	//k - index of current voxel
	//beginInd - beginning of partition
	//endInd - ending of partition
	//p - power

	float N = endInd - beginInd;
	double sum = 0;
	for (int i = beginInd; i < endInd; i++) {
		sum += std::pow(stresses[i].vonMises, p);
	}
	sum = std::pow(sum / N, (1 - p) / p) * (1.0 / N)*std::pow(stresses[k].vonMises, p - 1);

	return sum;
}

double calcVonMises(Eigen::VectorXd const  &stresses) {
	double sigma_ax = stresses[0]; double sigma_ay = stresses[1]; double sigma_az = stresses[2];
	double tau_xy = stresses[3]; double tau_yz = stresses[4]; double tau_xz = stresses[5];
	return std::sqrt(0.5*((sigma_ax - sigma_ay)*(sigma_ax - sigma_ay) + (sigma_ay - sigma_az)*(sigma_ay - sigma_az) + (sigma_az - sigma_ax)*(sigma_az - sigma_ax)) + 3 * (tau_xy*tau_xy + tau_yz * tau_yz + tau_xz * tau_xz));
}

Eigen::VectorXd calc_DsigmaVM_Dsigma_a(Stress const &stresses) {

	Eigen::VectorXd results(6);
	double sigma_ax = stresses.Stresses[0]; double sigma_ay = stresses.Stresses[1]; double sigma_az = stresses.Stresses[2];
	double tau_xy = stresses.Stresses[3]; double tau_yz = stresses.Stresses[4]; double tau_xz = stresses.Stresses[5];
	double vonMises = stresses.vonMises;
	results[0] = (sigma_ax * 2 - sigma_ay - sigma_az) / (vonMises * 2);
	results[1] = (sigma_ay * 2 - sigma_ax - sigma_az) / (vonMises * 2);
	results[2] = (sigma_az * 2 - sigma_ax - sigma_ay) / (vonMises * 2);
	results[3] = 3 * (tau_xy) / (vonMises);
	results[4] = 3 * (tau_yz) / (vonMises);
	results[5] = 3 * (tau_xz) / (vonMises);

	return results;
}

double  calc_DniS_Drho_e(double rho_e) {

	return 0.5*(1.0 / std::sqrt(rho_e));;

}

double calc_Drho_e_Dxb(int idx, int j, std::vector<designVariable> const &x, double r0, std::vector<std::vector<double>> const &nodes) {

	//idx - idx of design variable in nominator
	//j - idx of design variable in denominator
	//get distance between elements

	//for debugging

	//if (idx == j) {
		//	std::cout << "idx: " << idx << std::endl;
//		return 1;
//	}
	//return 0;

	auto pos_e = nodes[x[idx].varIdx];

	auto pos_b = nodes[x[j].varIdx];

	std::vector<double> temp = { pos_e[0] - pos_b[0], pos_e[1] - pos_b[1], pos_e[2] - pos_b[2] };

	double wj = (r0-length(temp))/r0;

	if (wj<0) {
		return 0;
	}

	return std::abs((wj / x[idx].sum_wj));

}

double ni_s(double rho) {
	return std::sqrt(rho);
}

double error(int n, double* x_new, double* x_old) {
	double sum = 0;
	for (int i = 0; i < n; i++) {
		sum += (x_new[i] - x_old[i])*(x_new[i] - x_old[i]);
	}

	return std::sqrt(sum);
}

std::vector<Stress> calcStresses(double ni, double E, std::vector<std::vector<double>> const nodes, std::vector<designVariable>  const x, Eigen::VectorXd  u, float dx) {

	std::cout << "calculating stress field" << std::endl;
	//define stress-index pair vector
	std::vector<Stress> stressVec;
	stressVec.resize(x.size());

	Eigen::MatrixXd B(6, 24), Elasticity(6, 6);
	Eigen::VectorXd  displacement(24);

	Elasticity << 1 - ni, ni, ni, 0, 0, 0,
		ni, 1 - ni, ni, 0, 0, 0,
		ni, ni, 1 - ni, 0, 0, 0,
		0, 0, 0, (1 - 2 * ni) / 2, 0, 0,
		0, 0, 0, 0, (1 - 2 * ni) / 2, 0,
		0, 0, 0, 0, 0, (1 - 2 * ni) / 2;

	Elasticity = E * Elasticity / ((1 + ni)*(1 - 2 * ni));

	Eigen::VectorXd stresses(nodes.size());
	std::vector<int> numStresses;

	numStresses.resize(nodes.size());

	for (int i = 0; i < x.size(); i++) {

		Eigen::VectorXd stress_sum = Eigen::VectorXd::Zero(6, 1);

		displacement(0) = u[x[i].nodeIdx[0] * 3]; displacement(1) = u[x[i].nodeIdx[0] * 3 + 1]; displacement(2) = u[x[i].nodeIdx[0] * 3 + 2];
		displacement(3) = u[x[i].nodeIdx[1] * 3]; displacement(4) = u[x[i].nodeIdx[1] * 3 + 1]; displacement(5) = u[x[i].nodeIdx[1] * 3 + 2];
		displacement(6) = u[x[i].nodeIdx[2] * 3]; displacement(7) = u[x[i].nodeIdx[2] * 3 + 1]; displacement(8) = u[x[i].nodeIdx[2] * 3 + 2];
		displacement(9) = u[x[i].nodeIdx[4] * 3]; displacement(10) = u[x[i].nodeIdx[4] * 3 + 1]; displacement(11) = u[x[i].nodeIdx[4] * 3 + 2];
		displacement(12) = u[x[i].nodeIdx[3] * 3]; displacement(13) = u[x[i].nodeIdx[3] * 3 + 1]; displacement(14) = u[x[i].nodeIdx[3] * 3 + 2];
		displacement(15) = u[x[i].nodeIdx[6] * 3]; displacement(16) = u[x[i].nodeIdx[6] * 3 + 1]; displacement(17) = u[x[i].nodeIdx[6] * 3 + 2];
		displacement(18) = u[x[i].nodeIdx[5] * 3]; displacement(19) = u[x[i].nodeIdx[5] * 3 + 1]; displacement(20) = u[x[i].nodeIdx[5] * 3 + 2];
		displacement(21) = u[x[i].nodeIdx[7] * 3]; displacement(22) = u[x[i].nodeIdx[7] * 3 + 1]; displacement(23) = u[x[i].nodeIdx[7] * 3 + 2];

		for (int j = 0; j < 8; j++) {
			float px, py, pz;
			switch (j) {
			case 0: px = -1;  py = -1;  pz = -1; break;
			case 1: px = 1;  py = -1;  pz = -1; break;
			case 2: px = -1;  py = 1;  pz = -1; break;
			case 3: px = -1;  py = -1;  pz = 1; break;
			case 4: px = 1;  py = 1;  pz = -1; break;
			case 5: px = 1;  py = -1;  pz = 1; break;
			case 6: px = -1;  py = 1;  pz = 1; break;
			case 7: px = 1;  py = 1;  pz = 1; break;
			}

		//	stress_sum = Elasticity * B_mat(px, py, pz)*displacement;
		//	stresses(x[i].nodeIdx[j]) += calcVonMises(stress_sum);
		//	numStresses[x[i].nodeIdx[j]] += 1;

			stressVec[i].idx = i;
			stressVec[i].Stresses = Elasticity * B_mat(0, 0, 0)*displacement/dx;
			stressVec[i].displacements = displacement;
			stressVec[i].vonMises = calcVonMises(stressVec[i].Stresses);

		}
	}
	/*
	std::cout << "done caculating stresses" << std::endl;

	//calculate stresses for each node and display with color graph using libigl

	for (int i = 0; i < nodes.size(); i++) {

		//		std::cout << "node " << i << ": " << stresses(i) << " " << numStresses[i] << std::endl;

		stresses(i) = stresses(i) / numStresses[i];

	}
	//return stresses;
	*/
	return stressVec;
}


std::vector<Stress> calcStresses(double ni, double E, std::vector<std::vector<double>> const nodes, std::vector<designVariable>  const x, Eigen::VectorXd  u, float dx, std::vector<int> targets) {

//	std::cout << "calculating stress field" << std::endl;
	//define stress-index pair vector
	std::vector<Stress> stressVec;
	stressVec.resize(x.size());

	Eigen::MatrixXd B(6, 24), Elasticity(6, 6);
	Eigen::VectorXd  displacement(24);

	Elasticity << 1 - ni, ni, ni, 0, 0, 0,
		ni, 1 - ni, ni, 0, 0, 0,
		ni, ni, 1 - ni, 0, 0, 0,
		0, 0, 0, (1 - 2 * ni) / 2, 0, 0,
		0, 0, 0, 0, (1 - 2 * ni) / 2, 0,
		0, 0, 0, 0, 0, (1 - 2 * ni) / 2;

	Elasticity = E * Elasticity / ((1 + ni)*(1 - 2 * ni));

	Eigen::VectorXd stresses(targets.size());

	for (int i = 0; i < targets.size(); i++) {

		Eigen::VectorXd stress_sum = Eigen::VectorXd::Zero(6, 1);

		displacement(0) = u[x[targets[i]].nodeIdx[0] * 3]; displacement(1) = u[x[targets[i]].nodeIdx[0] * 3 + 1]; displacement(2) = u[x[targets[i]].nodeIdx[0] * 3 + 2];
		displacement(3) = u[x[targets[i]].nodeIdx[1] * 3]; displacement(4) = u[x[targets[i]].nodeIdx[1] * 3 + 1]; displacement(5) = u[x[targets[i]].nodeIdx[1] * 3 + 2];
		displacement(6) = u[x[targets[i]].nodeIdx[2] * 3]; displacement(7) = u[x[targets[i]].nodeIdx[2] * 3 + 1]; displacement(8) = u[x[targets[i]].nodeIdx[2] * 3 + 2];
		displacement(9) = u[x[targets[i]].nodeIdx[4] * 3]; displacement(10) = u[x[targets[i]].nodeIdx[4] * 3 + 1]; displacement(11) = u[x[targets[i]].nodeIdx[4] * 3 + 2];
		displacement(12) = u[x[targets[i]].nodeIdx[3] * 3]; displacement(13) = u[x[targets[i]].nodeIdx[3] * 3 + 1]; displacement(14) = u[x[targets[i]].nodeIdx[3] * 3 + 2];
		displacement(15) = u[x[targets[i]].nodeIdx[6] * 3]; displacement(16) = u[x[targets[i]].nodeIdx[6] * 3 + 1]; displacement(17) = u[x[targets[i]].nodeIdx[6] * 3 + 2];
		displacement(18) = u[x[targets[i]].nodeIdx[5] * 3]; displacement(19) = u[x[targets[i]].nodeIdx[5] * 3 + 1]; displacement(20) = u[x[targets[i]].nodeIdx[5] * 3 + 2];
		displacement(21) = u[x[targets[i]].nodeIdx[7] * 3]; displacement(22) = u[x[targets[i]].nodeIdx[7] * 3 + 1]; displacement(23) = u[x[targets[i]].nodeIdx[7] * 3 + 2];

		for (int j = 0; j < 8; j++) {
			float px, py, pz;
			switch (j) {
			case 0: px = -1;  py = -1;  pz = -1; break;
			case 1: px = 1;  py = -1;  pz = -1; break;
			case 2: px = -1;  py = 1;  pz = -1; break;
			case 3: px = -1;  py = -1;  pz = 1; break;
			case 4: px = 1;  py = 1;  pz = -1; break;
			case 5: px = 1;  py = -1;  pz = 1; break;
			case 6: px = -1;  py = 1;  pz = 1; break;
			case 7: px = 1;  py = 1;  pz = 1; break;
			}

			//	stress_sum = Elasticity * B_mat(px, py, pz)*displacement;
			//	stresses(x[i].nodeIdx[j]) += calcVonMises(stress_sum);
			//	numStresses[x[i].nodeIdx[j]] += 1;

			stressVec[i].idx = i;
			stressVec[i].Stresses = Elasticity * B_mat(0, 0, 0)*displacement / dx;
			stressVec[i].displacements = displacement;
			stressVec[i].vonMises = calcVonMises(stressVec[i].Stresses);

		}
	}
	/*
	std::cout << "done caculating stresses" << std::endl;

	//calculate stresses for each node and display with color graph using libigl

	for (int i = 0; i < nodes.size(); i++) {

		//		std::cout << "node " << i << ": " << stresses(i) << " " << numStresses[i] << std::endl;

		stresses(i) = stresses(i) / numStresses[i];

	}
	//return stresses;
	*/
	return stressVec;
}

double optimizeStress(std::vector<Stress> stressVec, Eigen::VectorXd u, int C, int power, float dx, double change, std::vector<std::vector<double>> nodes, std::vector<designVariable>  x,
	double* xnew, double* dg, double* g, double* df, double* xmin, double *xmax, double density, double E, double ni, MMASolver* optimizer, int neighbourLayers,float r0) {

	//Elasticity
	Eigen::MatrixXd Elasticity(6, 6);

	Elasticity << 1 - ni, ni, ni, 0, 0, 0,
		ni, 1 - ni, ni, 0, 0, 0,
		ni, ni, 1 - ni, 0, 0, 0,
		0, 0, 0, (1 - 2 * ni) / 2, 0, 0,
		0, 0, 0, 0, (1 - 2 * ni) / 2, 0,
		0, 0, 0, 0, 0, (1 - 2 * ni) / 2;

	Elasticity = E * Elasticity / ((1 + ni)*(1 - 2 * ni));

	//optimization parameters
	double Xmax = 1;
	double Xmin = 0.001;
	double movLim = 0.02;

	//FEM parameters
	Eigen::MatrixXd KE = K_mat();

	//Material Parameters (TODO: need to move to store externally)
	double allowableStress = 320000000; 

	//local Ba 
	Eigen::MatrixXd  local_Ba(6, 24);
	local_Ba = (1.0 / dx)*B_mat(0, 0, 0); //at centeroid

	int deriv_Ind = 0;

	for (int i = 0; i < x.size(); i++) {

		Eigen::MatrixXd Ba = Eigen::MatrixXd::Zero(6, nodes.size() * 3);

		Ba.col(x[i].nodeIdx[0] * 3) = local_Ba.col(0);
		Ba.col(x[i].nodeIdx[0] * 3 + 1) = local_Ba.col(1);
		Ba.col(x[i].nodeIdx[0] * 3 + 2) = local_Ba.col(2);

		Ba.col(x[i].nodeIdx[1] * 3) = local_Ba.col(3);
		Ba.col(x[i].nodeIdx[1] * 3 + 1) = local_Ba.col(4);
		Ba.col(x[i].nodeIdx[1] * 3 + 2) = local_Ba.col(5);

		Ba.col(x[i].nodeIdx[2] * 3) = local_Ba.col(6);
		Ba.col(x[i].nodeIdx[2] * 3 + 1) = local_Ba.col(7);
		Ba.col(x[i].nodeIdx[2] * 3 + 2) = local_Ba.col(8);

		Ba.col(x[i].nodeIdx[4] * 3) = local_Ba.col(9);
		Ba.col(x[i].nodeIdx[4] * 3 + 1) = local_Ba.col(10);
		Ba.col(x[i].nodeIdx[4] * 3 + 2) = local_Ba.col(11);

		Ba.col(x[i].nodeIdx[3] * 3) = local_Ba.col(12);
		Ba.col(x[i].nodeIdx[3] * 3 + 1) = local_Ba.col(13);
		Ba.col(x[i].nodeIdx[3] * 3 + 2) = local_Ba.col(14);

		Ba.col(x[i].nodeIdx[6] * 3) = local_Ba.col(15);
		Ba.col(x[i].nodeIdx[6] * 3 + 1) = local_Ba.col(16);
		Ba.col(x[i].nodeIdx[6] * 3 + 2) = local_Ba.col(17);

		Ba.col(x[i].nodeIdx[5] * 3) = local_Ba.col(18);
		Ba.col(x[i].nodeIdx[5] * 3 + 1) = local_Ba.col(19);
		Ba.col(x[i].nodeIdx[5] * 3 + 2) = local_Ba.col(20);

		Ba.col(x[i].nodeIdx[7] * 3) = local_Ba.col(21);
		Ba.col(x[i].nodeIdx[7] * 3 + 1) = local_Ba.col(22);
		Ba.col(x[i].nodeIdx[7] * 3 + 2) = local_Ba.col(23);

		Eigen::VectorXd displacement(24);
		displacement(0) = u[x[i].nodeIdx[0] * 3]; displacement(1) = u[x[i].nodeIdx[0] * 3 + 1]; displacement(2) = u[x[i].nodeIdx[0] * 3 + 2];
		displacement(3) = u[x[i].nodeIdx[1] * 3]; displacement(4) = u[x[i].nodeIdx[1] * 3 + 1]; displacement(5) = u[x[i].nodeIdx[1] * 3 + 2];
		displacement(6) = u[x[i].nodeIdx[2] * 3]; displacement(7) = u[x[i].nodeIdx[2] * 3 + 1]; displacement(8) = u[x[i].nodeIdx[2] * 3 + 2];
		displacement(9) = u[x[i].nodeIdx[4] * 3]; displacement(10) = u[x[i].nodeIdx[4] * 3 + 1]; displacement(11) = u[x[i].nodeIdx[4] * 3 + 2];
		displacement(12) = u[x[i].nodeIdx[3] * 3]; displacement(13) = u[x[i].nodeIdx[3] * 3 + 1]; displacement(14) = u[x[i].nodeIdx[3] * 3 + 2];
		displacement(15) = u[x[i].nodeIdx[6] * 3]; displacement(16) = u[x[i].nodeIdx[6] * 3 + 1]; displacement(17) = u[x[i].nodeIdx[6] * 3 + 2];
		displacement(18) = u[x[i].nodeIdx[5] * 3]; displacement(19) = u[x[i].nodeIdx[5] * 3 + 1]; displacement(20) = u[x[i].nodeIdx[5] * 3 + 2];
		displacement(21) = u[x[i].nodeIdx[7] * 3]; displacement(22) = u[x[i].nodeIdx[7] * 3 + 1]; displacement(23) = u[x[i].nodeIdx[7] * 3 + 2];

		//calc stresses, add to stress and index pairs array

		stressVec[i].idx = i;
		stressVec[i].vonMises = calcVonMises(stressVec[i].Stresses);
		stressVec[i].displacements = displacement;
		x[i].displacements = displacement;

		if (x[i].tunable == false) {
			continue;
		}

		xmax[deriv_Ind] = Min(Xmax, xnew[deriv_Ind] + movLim);
		xmin[deriv_Ind] = Max(Xmin, xnew[deriv_Ind] - movLim);

		df[deriv_Ind] = dx*dx*dx*density; //mass of element

		deriv_Ind++;

	}

	//sort stresses

	std::cout << "sorting stresses" << std::endl;
	std::sort(stressVec.begin(), stressVec.end(), compareByLength);

	int beginInd = 0;

	int endInd;

	deriv_Ind = 0;

	SpMat K_global_full = build_global_stiffness_matrix_full(x, KE, power, nodes);
	K_global_full = dx *E* K_global_full;

	for (int i = 0; i < C; i++) { //i-th stress cluster

	//full global stiffness matrix

		Eigen::SparseMatrix<double> spMat;

		Eigen::VectorXd lambda_i = Eigen::VectorXd::Zero(nodes.size() * 3);

		Eigen::VectorXd  sum1 = Eigen::VectorXd::Zero(nodes.size() * 3);

		Eigen::VectorXd  sum2 = Eigen::VectorXd::Zero(nodes.size() * 3);

		double  sum3 = 0;

		Eigen::SimplicialCholesky<SpMat> chol(K_global_full);

		g[i] = 0;

		for (int k = beginInd; k < endInd; k++) {
			int idx = stressVec[k].idx;

			Eigen::VectorXd DsigmaVM_Dsigma_a(6);

			Eigen::MatrixXd Ba = Eigen::MatrixXd::Zero(6, nodes.size() * 3);

			Ba.col(x[idx].nodeIdx[0] * 3) = local_Ba.col(0);
			Ba.col(x[idx].nodeIdx[0] * 3 + 1) = local_Ba.col(1);
			Ba.col(x[idx].nodeIdx[0] * 3 + 2) = local_Ba.col(2);

			Ba.col(x[idx].nodeIdx[1] * 3) = local_Ba.col(3);
			Ba.col(x[idx].nodeIdx[1] * 3 + 1) = local_Ba.col(4);
			Ba.col(x[idx].nodeIdx[1] * 3 + 2) = local_Ba.col(5);

			Ba.col(x[idx].nodeIdx[2] * 3) = local_Ba.col(6);
			Ba.col(x[idx].nodeIdx[2] * 3 + 1) = local_Ba.col(7);
			Ba.col(x[idx].nodeIdx[2] * 3 + 2) = local_Ba.col(8);

			Ba.col(x[idx].nodeIdx[4] * 3) = local_Ba.col(9);
			Ba.col(x[idx].nodeIdx[4] * 3 + 1) = local_Ba.col(10);
			Ba.col(x[idx].nodeIdx[4] * 3 + 2) = local_Ba.col(11);

			Ba.col(x[idx].nodeIdx[3] * 3) = local_Ba.col(12);
			Ba.col(x[idx].nodeIdx[3] * 3 + 1) = local_Ba.col(13);
			Ba.col(x[idx].nodeIdx[3] * 3 + 2) = local_Ba.col(14);

			Ba.col(x[idx].nodeIdx[6] * 3) = local_Ba.col(15);
			Ba.col(x[idx].nodeIdx[6] * 3 + 1) = local_Ba.col(16);
			Ba.col(x[idx].nodeIdx[6] * 3 + 2) = local_Ba.col(17);

			Ba.col(x[idx].nodeIdx[5] * 3) = local_Ba.col(18);
			Ba.col(x[idx].nodeIdx[5] * 3 + 1) = local_Ba.col(19);
			Ba.col(x[idx].nodeIdx[5] * 3 + 2) = local_Ba.col(20);

			Ba.col(x[idx].nodeIdx[7] * 3) = local_Ba.col(21);
			Ba.col(x[idx].nodeIdx[7] * 3 + 1) = local_Ba.col(22);
			Ba.col(x[idx].nodeIdx[7] * 3 + 2) = local_Ba.col(23);

			DsigmaVM_Dsigma_a = calc_DsigmaVM_Dsigma_a(stressVec[k]);

			sum2 = sum2 + calc_DSigmaPN_DSigmaVM(stressVec, k, beginInd, endInd, power) * Ba.transpose()*Elasticity.transpose()*DsigmaVM_Dsigma_a;
		}

		//calculate lambda_i (adjoint variable)
		lambda_i = chol.solve(sum2);

		std::cout << "begin calculating derivatives with respect to tunable DV" << std::endl;

		for (int j = 0; j < x.size(); j++) { //derivative of stresses with respect to stress i

			std::cout << "calculated :" << j << "/" << x.size() << std::endl;

			if (x[j].tunable == false) { continue; }

			SpMat K_deriv((nodes.size()) * 3, (nodes.size()) * 3);

			float sum4 = 0;

			//build sparse matrix - the global stiffness derivative matrix
			for (int l = 0; l < x[j].influencesVoxels.size(); l++) {

				int voxelIdx = x[j].influencesVoxels[l];

				// Assembly:
				std::vector<T> coefficients;            // list of non-zeros coefficients
				std::vector<int> order = { 0,1,2,4,3,6,5,7 };
				double factor;
				factor = power * E* std::pow(x[voxelIdx].rho, power - 1);
				for (int i_ind = 0; i_ind < 8; i_ind++) {
					for (int j_ind = 0; j_ind < 8; j_ind++) {

						int row1 = x[voxelIdx].nodeIdx[i_ind] * 3; int row2 = x[voxelIdx].nodeIdx[i_ind] * 3 + 1; int row3 = x[voxelIdx].nodeIdx[i_ind] * 3 + 2;
						int col1 = x[voxelIdx].nodeIdx[j_ind] * 3; int col2 = x[voxelIdx].nodeIdx[j_ind] * 3 + 1; int col3 = x[voxelIdx].nodeIdx[j_ind] * 3 + 2;
						double val11 = KE(order[i_ind] * 3, order[j_ind] * 3)*factor; double val12 = KE(order[i_ind] * 3, order[j_ind] * 3 + 1)*factor; double val13 = KE(order[i_ind] * 3, order[j_ind] * 3 + 2)*factor;
						double val21 = KE(order[i_ind] * 3 + 1, order[j_ind] * 3)*factor; double val22 = KE(order[i_ind] * 3 + 1, order[j_ind] * 3 + 1)*factor; double val23 = KE(order[i_ind] * 3 + 1, order[j_ind] * 3 + 2)*factor;
						double val31 = KE(order[i_ind] * 3 + 2, order[j_ind] * 3)*factor; double val32 = KE(order[i_ind] * 3 + 2, order[j_ind] * 3 + 1)*factor; double val33 = KE(order[i_ind] * 3 + 2, order[j_ind] * 3 + 2)*factor;

						coefficients.push_back(T(row1, col1, val11)); coefficients.push_back(T(row1, col2, val12)); coefficients.push_back(T(row1, col3, val13));
						coefficients.push_back(T(row2, col1, val21)); coefficients.push_back(T(row2, col2, val22)); coefficients.push_back(T(row2, col3, val23));
						coefficients.push_back(T(row3, col1, val31)); coefficients.push_back(T(row3, col2, val32)); coefficients.push_back(T(row3, col3, val33));
					}
				}

				K_deriv.setFromTriplets(coefficients.begin(), coefficients.end());

				sum1 += dx * K_deriv*calc_Drho_e_Dxb(j, voxelIdx, x, r0, nodes)*u;

				coefficients.clear();

			}

			endInd = std::floor(((i + 1)*stressVec.size()) / C); // begining and ending index for stress group in stress array

			//build global Ba

			for (int k = beginInd; k < endInd; k++) {
				int idx = stressVec[k].idx;

				Eigen::VectorXd DsigmaVM_Dsigma_a(6);

				Eigen::MatrixXd Ba = Eigen::MatrixXd::Zero(6, nodes.size() * 3);

				Ba.col(x[idx].nodeIdx[0] * 3) = local_Ba.col(0);
				Ba.col(x[idx].nodeIdx[0] * 3 + 1) = local_Ba.col(1);
				Ba.col(x[idx].nodeIdx[0] * 3 + 2) = local_Ba.col(2);

				Ba.col(x[idx].nodeIdx[1] * 3) = local_Ba.col(3);
				Ba.col(x[idx].nodeIdx[1] * 3 + 1) = local_Ba.col(4);
				Ba.col(x[idx].nodeIdx[1] * 3 + 2) = local_Ba.col(5);

				Ba.col(x[idx].nodeIdx[2] * 3) = local_Ba.col(6);
				Ba.col(x[idx].nodeIdx[2] * 3 + 1) = local_Ba.col(7);
				Ba.col(x[idx].nodeIdx[2] * 3 + 2) = local_Ba.col(8);

				Ba.col(x[idx].nodeIdx[4] * 3) = local_Ba.col(9);
				Ba.col(x[idx].nodeIdx[4] * 3 + 1) = local_Ba.col(10);
				Ba.col(x[idx].nodeIdx[4] * 3 + 2) = local_Ba.col(11);

				Ba.col(x[idx].nodeIdx[3] * 3) = local_Ba.col(12);
				Ba.col(x[idx].nodeIdx[3] * 3 + 1) = local_Ba.col(13);
				Ba.col(x[idx].nodeIdx[3] * 3 + 2) = local_Ba.col(14);

				Ba.col(x[idx].nodeIdx[6] * 3) = local_Ba.col(15);
				Ba.col(x[idx].nodeIdx[6] * 3 + 1) = local_Ba.col(16);
				Ba.col(x[idx].nodeIdx[6] * 3 + 2) = local_Ba.col(17);

				Ba.col(x[idx].nodeIdx[5] * 3) = local_Ba.col(18);
				Ba.col(x[idx].nodeIdx[5] * 3 + 1) = local_Ba.col(19);
				Ba.col(x[idx].nodeIdx[5] * 3 + 2) = local_Ba.col(20);

				Ba.col(x[idx].nodeIdx[7] * 3) = local_Ba.col(21);
				Ba.col(x[idx].nodeIdx[7] * 3 + 1) = local_Ba.col(22);
				Ba.col(x[idx].nodeIdx[7] * 3 + 2) = local_Ba.col(23);

				DsigmaVM_Dsigma_a = calc_DsigmaVM_Dsigma_a(stressVec[k]);

				sum2 = sum2 + std::pow(x[stressVec[k].idx].rho, 0.5) *calc_DSigmaPN_DSigmaVM(stressVec, k, beginInd, endInd, power) * Ba.transpose()*Elasticity.transpose()*DsigmaVM_Dsigma_a;

				sum3 += calc_DSigmaPN_DSigmaVM(stressVec, k, beginInd, endInd, power)*DsigmaVM_Dsigma_a.transpose()*calc_DniS_Drho_e(x[stressVec[k].idx].rho)*calc_Drho_e_Dxb(j, stressVec[k].idx, x, r0, nodes)*Elasticity*Ba * u;

			}

			/*
			if (deriv_Ind == 0) {
				std::cout << "sum3: " << sum3 << std::endl;
				std::cout << "sum1: " << sum1 << std::endl;
				std::cout << "lambda_i: " << lambda_i << std::endl;
				exit(-1);
			}
			*/
			dg[deriv_Ind] = sum3 - std::pow(x[stressVec[j].idx].rho, 0.5)*lambda_i.transpose()*sum1;

			deriv_Ind++;

		}

		for (int k = beginInd; k < endInd; k++) {
			g[i] += std::pow(stressVec[k].vonMises, 3);

		}

		g[i] = std::pow((double)1.0 / (endInd - beginInd)*g[i], (double)1.0 / 3) - allowableStress;


		beginInd = endInd;

		std::cout << "finished calculating tunable derivatives" << std::endl;
	}
			
				std::cout << "done calculating derivatives - starting mma" << std::endl;
				/*
				for (int i = 0; i < x.size(); i++) {
					std::cout << "df[" << i << "]: " << df[i] << std::endl;
				}
				for (int i = 0; i < C; i++) {
					std::cout << "g[" << i << "]: " << g[i] << std::endl;
				}
				for (int i = 0; i < x.size()*C; i++) {
					std::cout << "dg[" << i << "]: " << dg[i] << std::endl;
				}
				*/
				optimizer->Update(xnew, df, g, dg, xmin, xmax);
				int tempInd = 0;

				change = -1;
				for (int i = 0; i < x.size(); i++) {
					if (x[i].tunable == false) {
						continue;
					}
					change = Max(change, (std::abs(x[i].value - xnew[tempInd]))/xnew[tempInd]);
					x[i].value = xnew[tempInd];
			//		std::cout << "x[" << tempInd << "]: " << xnew[tempInd] << std::endl;
					tempInd++;

				}

				return change;
}


