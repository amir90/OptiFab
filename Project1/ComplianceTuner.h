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

	double movLim = 0.2;
	int ind = 0;
	for (int i = 0; i < x.size(); i++) {

		if (!x[i].tunable) {
			continue;
		}

		xmax[ind] = std::min(1.0, xnew[ind] + movLim);
		xmin[ind] = std::max(0.001, xnew[ind] - movLim);
		ind++;
	}

	mma->Update(xnew, dc, g, dg, xmin, xmax);

	double change = -1;
	ind = 0;
	g[0] = -V_allowed;
	for (int i = 0; i < x.size(); i++) {
		if (!x[i].tunable) {
			continue;
		}

		change = Max(change, std::abs(x[i].value - xnew[ind]) / x[i].value);
		x[i].value = xnew[ind];

		ind++;
	}
	//	std::cout <<"element "<<i<<": "<< x[i].value << std::endl;

	//calc_rho(nodes, x, 1, dx, numOfVoxelsX, numOfVoxelsY, numOfVoxelsZ,r0);

	for (int i = 0; i < x.size(); i++) {

		if (x[i].tunable) {
			g[0] += x[i].value;
		}
	}

	return change;
	//	makeMesh(x, nodes);
}