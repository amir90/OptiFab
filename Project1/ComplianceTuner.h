#pragma once
#include "FEMSimulator.h" //need to generalize this

//function that performs the optimization step

double optimizeCompliance(double* xnew, double* dg, double* g, double* dc, double* xmin, double *xmax, double* dc_filtered, double V_allowed, MMASolver* mma, std::vector<std::vector<double>>& nodes, std::vector<designVariable> &x, double dx,
	double numOfVoxelsX, double numOfVoxelsY, double numOfVoxelsZ, float r0)
{

	//calculate sensitivty filter
//		calc_sensitivty_filter(dc_filtered, dc, nodes, x, 2, dx, numOfVoxelsX, numOfVoxelsY, numOfVoxelsZ);

//	for (int k = 0; k < x.size(); k++) {
//		std::cout << "testing dc: " << dc[k] << std::endl;
//	}

	mma->Update(xnew, dc, g, dg, xmin, xmax);

	double change = -1;
	g[0] = -V_allowed;
	for (int i = 0; i < x.size(); i++) {
		change = Max(change, std::abs(x[i].value - xnew[i]) / x[i].value);
		x[i].value = xnew[i];
	}
	//	std::cout <<"element "<<i<<": "<< x[i].value << std::endl;

	calc_rho(nodes, x, 1, dx, numOfVoxelsX, numOfVoxelsY, numOfVoxelsZ,r0);

	for (int i = 0; i < x.size(); i++) {

		g[0] += x[i].value;
	}

	return change;
	//	makeMesh(x, nodes);
}