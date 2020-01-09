#include "fatigueOptimizer.h"

extern bool debug_flag;

void updateMaxStresses(Eigen::VectorXd tempStress, Eigen::VectorXd &maxStress, Eigen::VectorXd &minStress) {

	if (tempStress(0) > maxStress(0)) {

		maxStress(0) = tempStress(0);

	}

	if (tempStress(1) > maxStress(1)) {

		maxStress(1) = tempStress(1);

	}

	if (tempStress(2) > maxStress(2)) {

		maxStress(2) = tempStress(2);

	}

	if (tempStress(3) > maxStress(3)) {

		maxStress(3) = tempStress(3);

	}

	if (tempStress(4) > maxStress(4)) {

		maxStress(4) = tempStress(4);

	}

	if (tempStress(5) > maxStress(5)) {

		maxStress(5) = tempStress(5);

	}

	if (tempStress(0) < minStress(0)) {

		minStress(0) = tempStress(0);

	}

	if (tempStress(1) < minStress(1)) {

		minStress(1) = tempStress(1);

	}

	if (tempStress(2) < minStress(2)) {

		minStress(2) = tempStress(2);

	}

	if (tempStress(3) < minStress(3)) {

		minStress(3) = tempStress(3);

	}

	if (tempStress(4) < minStress(4)) {

		minStress(4) = tempStress(4);

	}

	if (tempStress(5) < minStress(5)) {

		minStress(5) = tempStress(5);

	}

}


std::vector<std::vector<Stress>>  ClacStressTimeHistory(std::vector<int> targetElements, std::vector<designVariable> x, double a, double b, int power, std::set < std::pair<int, Eigen::Vector3f>, comp> forceVertexIdSet, double dx, double E, std::vector<std::vector<double>> nodes, std::set<int> constraintVertexIdSet, Eigen::VectorXd u_dot_0, Eigen::VectorXd u0, double dt, int steps, std::vector<std::vector<double>> &sigma_m, std::vector<std::vector<Stress>> &sigma_a,Eigen::MatrixXd& u_test)
{

	u_test = Eigen::MatrixXd(steps, 3*nodes.size());

	std::vector<int> freeDOF;
	freeDOF.resize(nodes.size());
	sigma_m.resize(targetElements.size());
	sigma_a.resize(targetElements.size());
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

	Eigen::VectorXd f(3 * ind);

	for (int j = 0; j < 3 * ind; j++) {
		f(j) = 0;
	}

	for (auto j = forceVertexIdSet.begin(); j != forceVertexIdSet.end(); j++) {
		int tempForceid = freeDOF[j->first] * 3;
		f[tempForceid + 0] = j->second[0];// force[temp_f].second[0];
		f[tempForceid + 1] = j->second[1];//force[temp_f].second[1];
		f[tempForceid + 2] = j->second[2]; //force[temp_f].second[2];//

	}

	int numOfFixed = nodes.size() - ind; //number of fixed nodes

	SpMat K = build_global_stiffness_matrix(x, nodes.size(), K_mat(), power, freeDOF, ind)*E*dx;

	//need to build the mass matrix

//	std::cout << "finished building global stiffness matrix" << std::endl;

	SpMat M = build_global_mass_matrix(x, nodes.size(), M_mat() , power, freeDOF, ind )*dx;

//	std::cout << "finished building global mass matrix" << std::endl;

	//build the damping matrix = linear combiniation of mass and stiffness matrices.

	SpMat C = M * a;

	C += K * b; //better performance-wise according to eigen sparse matrix documentation

	Eigen::VectorXd u_prev(3*ind), u_prev_dot(3 * ind), u_prev_dotdot(3 * ind);
	Eigen::VectorXd u_curr(3 * ind), u_curr_dot(3 * ind), u_curr_dotdot(3 * ind);
	//Mu''(t) + Cu'(t) + Ku(t) = f(t)
	
	//calc u''

	//u = u0; u_dot = u_dot_0; //for time step 0

	//initialize newmarks method

	u_curr = Eigen::VectorXd::Zero(3*ind);

	u_curr_dot = Eigen::VectorXd::Zero(3*ind);

	u_prev = Eigen::VectorXd::Zero(3 * ind);

	u_prev_dot = Eigen::VectorXd::Zero(3 * ind);

	u_prev_dotdot = Eigen::VectorXd::Zero(3 * ind);

	Eigen::VectorXd A = f - K * u_curr - C * u_curr_dot ;

	Eigen::SimplicialCholesky<SpMat> chol(M);

	double beta = 0.25;
	double alpha = 0.5;

	Eigen::SimplicialCholesky<SpMat> chol2(M + alpha*dt*C + beta*dt*dt*K);

	u_prev_dotdot = chol.solve(A);

	//test: calculate u''(t)

	//std::cout << u_dotdot << std::endl; //tested

	//Eigen::VectorXd u_test(steps);

	std::vector<std::vector<Stress>> StressMat; //vector of stress vectors

	StressMat.resize(targetElements.size());

	for (int i = 0; i < steps; i++) {

		//newmarks method
		
		A = -K * (u_prev + dt * u_prev_dot + (1 - 2*beta) / 2.0*dt*dt*u_prev_dotdot) - C * (u_prev_dot + (1 - alpha)*dt*u_prev_dotdot);

		u_curr_dotdot = chol2.solve(A);

		u_curr_dot = u_prev_dot + dt * (1 - alpha) * u_prev_dotdot + dt * alpha*u_curr_dotdot;

		u_curr = u_prev + dt * u_prev_dot + (1 - 2 * beta) / 2.0 * dt*dt*u_prev_dotdot + beta * dt*dt*u_curr_dotdot;

		u_prev_dotdot = u_curr_dotdot;

		u_prev_dot = u_curr_dot;

		u_prev = u_curr;

		/* dynamic FEM tested and works
		Eigen::VectorXd u_check = doFEM(nodes, x, forceVertexIdSet, constraintVertexIdSet, dx, 1);

		std::cout << u_check(3 * ind - 1) << std::endl;
		*/
		
		//calculate von mises stresses for required elements

		Eigen::VectorXd u(3*nodes.size());

		// create global displacement vector
		int u_ind = 0; 
		for (int k = 0; k < nodes.size(); k++) {

			if (freeDOF[k] != -1) {
				u[k * 3] = u_curr[u_ind++];
				u[k * 3 + 1] = u_curr[u_ind++];
				u[k * 3 + 2] = u_curr[u_ind++];
			}
			else {
				u[k * 3] = 0;
				u[k * 3 + 1] = 0;
				u[k * 3 + 2] = 0;
			}
		}

		u_test.row(i) = u.transpose();
		
		auto stresses = calcStresses(0.3, 1, nodes, x, u, dx, targetElements);

	//	std::cout << "done calculating stress" << std::endl;

		for (int j = 0; j < targetElements.size(); j++) { //stresses at time t at target elemetns
			
			StressMat[j].push_back(stresses[j]); 

		}
		
	} //loop for all time steps

	//for debugging: inject stress Mat

	
	std::ofstream debugFile;

	if (debug_flag) {

		//std::vector<double> stressList = {0, 80, 125, 50, 20,140, 275, 240, 125, 115, 170, 140, 160, 60};

		std::vector<double> stressListX = { 0, 70, 120, 51, 30, 75, 100, 20, -60, -50, 40, 45, 0, 20 };

		std::vector<double> stressListY = { 0, 30, 40, 20, 10, 32, 38, 8, -28, -18, 15, 18, -3, 8 };

		std::vector<double> stressListXY = { 0, 50, 75, 40, 20, 52, 60, 8, -50, -40, 20, 25, -3, 5 };


		StressMat[0].clear();
		Stress tempStress;
		for (int i = 0; i < stressListX.size(); i++) {

			tempStress.Stresses = Eigen::VectorXd::Zero(6);
			tempStress.Stresses(0) = stressListX[i];
			tempStress.Stresses(1) = stressListY[i];
			tempStress.Stresses(3) = stressListXY[i];
			tempStress.vonMises = calcVonMises(tempStress.Stresses);
			StressMat[0].push_back(tempStress);

		}

		std::cout << "debug_flag is true" << std::endl;
	}

//	std::cout << "begin multi axial cycle counting" << std::endl;

	for (int i = 0; i < targetElements.size(); i++) { //for every target element

		//MultiAxial cycle counting method

		//find largest von mises stress
		int minS=0;

		for (int j = 1; j < StressMat[i].size(); j++) {

			if (StressMat[i][minS].vonMises < StressMat[i][j].vonMises) {

				minS = j;

			}

		}

	//	std::cout << "maximum is: " << minS << std::endl;

		//create array that begins with largest stress (relative time array).

		std::vector<Stress> relArr(StressMat[i].size());

		for (int j = 0; j < StressMat[i].size(); j++) {
			int curr_id = (j + minS) % StressMat[i].size();
			relArr[j] = StressMat[i][curr_id];

		}

		std::vector<int> cycleArr(relArr.size(),-1); //keeps track which node is in which cycle: cycleArr[node id] = cycle id

		int numOfStressesLeft = relArr.size(); //number of stresses not assigned a cycle

		int numOfCurrCycle = 0;

		//Cycles.push_back(std::pair<int, int>(0, StressMat[i].size())); //initialize with all the stresses

		std::set<int> blockIndices; //keeps tracks of where a block needs to begin

		blockIndices.insert(0); //begin at time t=0

		std::vector<Stress> sigma_a_raw;

		std::vector<Stress> sigma_m_raw;

		//change to relative time vector

		while (blockIndices.size() > 0) { //loop until all elements are in a cycle

			Eigen::VectorXd stressMax(6);
			Eigen::VectorXd stressMin(6);

			//calculate equivalent stresses for each target element

			int t0 = *blockIndices.begin();

			blockIndices.erase(t0);

			int offset = t0; //where the block really begins

			t0 = t0 - 1;

//			std::cout << "offset: " << offset << std::endl;

			int blockSize = 0; //how many stresses there are in the "block"

			if (numOfCurrCycle == 0) {
				t0 = 0;
				blockSize = relArr.size();
			}
			else {

				for (int j = 0; j < relArr.size(); j++) { //find maximum in block and set offset to it

					blockSize++;

					int curr_idx = (j + offset) % relArr.size();

					if (cycleArr[curr_idx] != -1) {
						break;
					}

//					std::cout << "stress in " << curr_idx << ": " << relArr[curr_idx].vonMises << std::endl;


				}

				blockSize--;

			}
//			std::cout << "t0 index: " << t0 << std::endl;

		//	std::cout << "blocksize: " << blockSize << std::endl;

			//t0 is now the start of the new cycle

			cycleArr[t0] = numOfCurrCycle;

			stressMax(0) = relArr[t0].Stresses(0); stressMin(0) = relArr[t0].Stresses(0);

			stressMax(1) = relArr[t0].Stresses(1); stressMin(1) = relArr[t0].Stresses(1);

			stressMax(2) = relArr[t0].Stresses(2); stressMin(2) = relArr[t0].Stresses(2);

			stressMax(3) = relArr[t0].Stresses(3); stressMin(3) = relArr[t0].Stresses(3);

			stressMax(4) = relArr[t0].Stresses(4); stressMin(4) = relArr[t0].Stresses(4);

			stressMax(5) = relArr[t0].Stresses(5); stressMin(5) = relArr[t0].Stresses(5);

			int NumOfStresses = blockSize;//StressMat[i].size();

			double delta_curr;

			double delta_prev=0;

			int lastIdx = 0;

			int currIdx = (t0 - offset + 1) % NumOfStresses+offset;

//			std::cout << "delta: " << delta_curr << std::endl;

			bool blockFlag = false; //are you in a block?

			while (cycleArr[currIdx]==-1   && currIdx!=t0) {

				delta_curr = calcVonMises(relArr[currIdx].Stresses - relArr[t0].Stresses);

//				std::cout << "current index: " << currIdx <<" delta: "<< delta_curr<< " block: "<< blockFlag <<std::endl;

				//Collect all points that are not decreasing - this is the first major reversal

				if ((delta_curr - delta_prev)/delta_curr>= 0.001) {

//					std::cout <<"delta test: "<< (delta_curr - delta_prev) / delta_prev << std::endl;

					cycleArr[currIdx] = numOfCurrCycle; 

					delta_prev = delta_curr;

			//		updateMaxStresses(relArr[currIdx].Stresses, stressMax, stressMin);

//					std::cout <<"biggest delta: "<< delta_curr << std::endl;

					blockFlag = false;

					lastIdx = currIdx;
				}
				else {

					if (!blockFlag) { //save index as start of new block

						blockIndices.insert(currIdx);

					}

					blockFlag = true;

				}

				currIdx = (currIdx-offset + 1) % NumOfStresses + offset;

			}


			numOfCurrCycle++;

			
//			for (int j = 0; j < cycleArr.size(); j++) {

//				std::cout << cycleArr[j] << std::endl;
//			}

//			std::cout << "......" << std::endl;

//			for (auto j = blockIndices.begin(); j != blockIndices.end(); j++) {

//				std::cout << *j << std::endl;
//			}
			

//			std::cout << "calc sigma_a and sigma_m" << std::endl;

//			std::cout << t0 << "\n" << currIdx-1 << std::endl;

			Stress tempStress;

			tempStress.Stresses = (-relArr[t0].Stresses + relArr[lastIdx].Stresses) / 2;

			tempStress.vonMises = calcVonMises(tempStress.Stresses);

			sigma_a[i].push_back(tempStress);

			Eigen::VectorXd meanStress = (relArr[t0].Stresses + relArr[lastIdx].Stresses)/2;

			sigma_m[i].push_back(meanStress(0)+ meanStress(1) + meanStress(2));

		}

		if (debug_flag) {
			break;
		}

	}

	return StressMat;
}



double dD_PN_k_dgamma_e(int e_id, std::vector<designVariable> x, std::vector<int> Omega_k, double p, std::vector<double> damage, double r0, std::vector<std::vector<double>> nodes) {

	// calculate the derivative of damage cluster D_PN_k with respect to the value of the voxel e
	//inputs: 
		// e_id - index of voxel e
		// x - design variable (voxel) vector. includes all topological and geomterical information about voxel grid.
		// Omega_k - vector of voxels indices in damage cluster k
		// damage - vector of damages and ids.

	double sum1 = 0;
	double sum2 = 0;

	for (int i = 0; i < Omega_k.size(); i++) {
		int currId = Omega_k[i];
		sum1 += damage[currId] * x[currId].rho;
			if (x[currId].influencesVoxels.find(i) != x[currId].influencesVoxels.end()) {
				sum2 += calc_Drho_e_Dxb(e_id, i, x, r0, nodes)*std::pow(damage[currId],p);
			}
	}

	return (1.0 / p)*std::pow(sum1, 1.0 / p - 1)*sum2;
	
}

double dD_PN_k_dDe(int e_id, double p, std::vector<designVariable> x,std::vector<int> Omega_k, std::vector<double> damage) {
	// tested and outputs a number

	double sum1 = 0;

	for (int i = 0; i < Omega_k.size(); i++) {
		int currId = Omega_k[i];
		sum1 += damage[currId] * x[currId].rho;
	}

	//note: if sum1 is close to 0, numerically unstable. TODO: need to normalize
	return (p - 1.0 / p)*std::pow(sum1, 1.0 / p - 1)*std::pow(damage[e_id], p - 1)*x[e_id].rho;

}

double dDe_dsigma_a_i(double Sut, double bf, double Sf, std::vector<double> sigma_a, std::vector<double> sigma_m, int cycle, std::vector<int> n) {

	double sum = 0;

	int eps; //small constant for approximation of Psi_(2,max)

	int numOfCycles = sigma_a.size();

	double factor;
	
	for (int i = 0; i < numOfCycles; i++) {

		double Psi = sigma_m[i] / 2 + std::sqrt(sigma_m[i]*sigma_m[i]+ 0.00001) / 2;

		if (i == cycle) {

			factor = 2 * n[i] * std::pow((Sut*sigma_a[i] / (Sut - Psi)), -1.0 / bf);

			sum += factor;

		} else {

			sum += 2 * n[i] * std::pow((Sut*sigma_a[i] / (Sut - Psi)), -1.0 / bf);
		}

	}

	return std::pow(sum, 1.0 / Sf - 1.0)*1.0 / Sf * factor / (-bf * sigma_a[cycle]);

}

Eigen::VectorXd dsigma_a_i_vec_dgamma_e(int e_id, double p, std::vector<designVariable> x, Eigen::VectorXd u_max, Eigen::VectorXd u_min, double E) {

	return 0.5*calc_DniS_Drho_e(x[e_id].rho)*E*std::pow(x[e_id].rho,p)*B_mat(0,0,0)*(u_max-u_min);

}

double dDe_dsigma_m_i (double Sut, double bf, double Sf, std::vector<double> sigma_a, std::vector<double> sigma_m, int cycle, std::vector<int> n) {

	double sum = 0;

	int eps; //small constant for approximation of Psi_(2,max)

	int numOfCycles = sigma_a.size();

	double factor;
	
	double dPsi;

	for (int i = 0; i < numOfCycles; i++) {

		double Psi = sigma_m[i] / 2 + std::sqrt(sigma_m[i] * sigma_m[i] + 0.00001) / 2;

		if (i == cycle) {

			dPsi = 0.5 - 0.5*sigma_m[i] / std::sqrt(sigma_m[i] * sigma_m[i] + eps);

			factor = 2*n[i]*std::pow(Sut*sigma_a[i]/(Sut-Psi),-1.0/bf)/(Sut - Psi)*(1/Sf)*(1/bf) ;

		}

			sum += 2 * n[i] * std::pow((Sut*sigma_a[i] / (Sut - Psi)), -1.0 / bf);
		}

		return sum * factor*dPsi;

	}


inline Eigen::RowVectorXd dsigma_m_i_dsigma_m_i_vec() {

	Eigen::RowVectorXd vec(6);

	vec << 1, 1, 1, 0, 0, 0;
	return  vec;

}


Eigen::Vector3d dsigma_m_i_vec_dgamma_e(int e_id, double p, std::vector<designVariable> x, Eigen::VectorXd u_max, Eigen::VectorXd u_min) {

	return 0.5*calc_DniS_Drho_e(x[e_id].rho)*B_mat(0, 0, 0)*(u_max + u_min);

}


Eigen::VectorXd Lambda_max_i(int e_id, double p, std::vector<designVariable> x, std::vector<int> Omega_k, std::vector<Stress> sigma_a, std::vector<double> damage, double Sut, double bf, double Sf, std::vector<double> sigma_m, double E, int cycle, std::vector<int> n) {

	//tested and produces a number.

	std::vector<double> Sigma_a_vonMises(sigma_a.size());

	for (int i=0; i<sigma_a.size(); i++) {

		Sigma_a_vonMises[i] = calcVonMises(sigma_a[i].Stresses);

	}

	Eigen::VectorXd sum = Eigen::VectorXd::Zero(24);
	for (int i = 0; i < Omega_k.size(); i++) {
		sum += dD_PN_k_dDe(i, p, x, Omega_k, damage)*(dDe_dsigma_a_i(Sut, bf, Sf, Sigma_a_vonMises, sigma_m, cycle, n)*calc_DsigmaVM_Dsigma_a(sigma_a[i]).transpose()*ni_s(x[i].rho)*E*std::pow(x[i].rho, p)*(K_mat()*B_mat(0, 0, 0).transpose()).transpose() +
			dDe_dsigma_m_i(Sut, bf, Sf, Sigma_a_vonMises, sigma_m, cycle, n)*dsigma_m_i_dsigma_m_i_vec()*E*std::pow(x[i].rho, p)*(K_mat()*B_mat(0, 0, 0).transpose()).transpose());

		std::cout << sum << std::endl;
	}


	sum = sum / 2;
	Eigen::MatrixXd K = E * std::pow(x[e_id].rho, p)*K_mat();
	Eigen::VectorXd lambda = K.inverse()*sum;
	return lambda;

}


Eigen::VectorXd Lambda_min_i(int e_id, double p, std::vector<designVariable> x, std::vector<int> Omega_k, std::vector<Stress> sigma_a, std::vector<double> damage, double Sut, double bf, double Sf, std::vector<double> sigma_m, double E, int cycle, std::vector<int> n) {

	std::vector<double> Sigma_a_vonMises(sigma_a.size());

	for (int i = 0; i < sigma_a.size(); i++) {

		Sigma_a_vonMises[i] = calcVonMises(sigma_a[i].Stresses);

	}

	Eigen::VectorXd sum = Eigen::VectorXd::Zero(24);
	for (int i = 0; i < Omega_k.size(); i++) {
		sum += dD_PN_k_dDe(i, p, x, Omega_k, damage)*(-dDe_dsigma_a_i(Sut, bf, Sf, Sigma_a_vonMises, sigma_m, cycle, n)*calc_DsigmaVM_Dsigma_a(sigma_a[i]).transpose()*ni_s(x[i].rho)*E*std::pow(x[i].rho, p)*(K_mat()*B_mat(0, 0, 0).transpose()).transpose() +
			dDe_dsigma_m_i(Sut, bf, Sf, Sigma_a_vonMises, sigma_m, cycle, n)*dsigma_m_i_dsigma_m_i_vec()*E*std::pow(x[i].rho, p)*(K_mat()*B_mat(0, 0, 0).transpose()).transpose());

	}

	Eigen::MatrixXd K = E * std::pow(x[e_id].rho, p)*K_mat();
	sum = K.inverse()*sum / 2;
	Eigen::VectorXd lambda = K.inverse()*sum;
	return lambda;

}
