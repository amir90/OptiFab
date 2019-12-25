
#include<iostream>
#include<fstream>
#include <MMASolver.h>
#include<igl\cotmatrix.h>
#include<igl\readOFF.h>
#include <igl\readOBJ.h>
#include <igl\readMESH.h>
#include<igl\opengl\glfw\Viewer.h>
#include "HelperTop.h"
//#include <igl/jet.h>




double Min(double d1, double d2) {
	return d1 < d2 ? d1 : d2;
}

double Max(double d1, double d2) {
	return d1 > d2 ? d1 : d2;
}

int Min(int d1, int d2) {
	return d1 < d2 ? d1 : d2;
}

int Max(int d1, int d2) {
	return d1 > d2 ? d1 : d2;
}


int moveX(std::vector<designVariable> const &x, int curridx, int move) {
	//return the voxel id of voxel that is "move" voxels in x direction from voxel with curridx

	if (move == 0) { return curridx; }
	int tempidx = curridx;
	if (move > 0) {
		for (int i = 0; i < move; i++) {
			tempidx = x[tempidx].neighbors[4];
			if (tempidx == -1) {
				break;
			}
		}
	}
	else {
		for (int i = 0; i < -1 * move; i++) {
			tempidx = x[tempidx].neighbors[5];
			if (tempidx == -1) {
				break;
			}
		}
	}

	return tempidx;
}

int moveY(std::vector<designVariable> const &x, int curridx, int move) {
	//return the voxel id of voxel that is "move" voxels in y direction from voxel with curridx

	if (move == 0) { return curridx; }
	int tempidx = curridx;
	if (move > 0) {
		for (int i = 0; i < move; i++) {
			tempidx = x[tempidx].neighbors[2];
			if (tempidx == -1) {
				break;
			}
		}
	}
	else {
		for (int i = 0; i < -1 * move; i++) {
			tempidx = x[tempidx].neighbors[3];
			if (tempidx == -1) {
				break;
			}
		}
	}

	return tempidx;
}

int moveZ(std::vector<designVariable> const &x, int curridx, int move) {
	//return the voxel id of voxel that is "move" voxels in z direction from voxel with curridx

	if (move == 0) { return curridx; }
	int tempidx = curridx;
	if (move > 0) {
		for (int i = 0; i < move; i++) {
			tempidx = x[tempidx].neighbors[0];
			if (tempidx == -1) {
				break;
			}
		}
	}
	else {
		for (int i = 0; i < -1 * move; i++) {
			tempidx = x[tempidx].neighbors[1];
			if (tempidx == -1) {
				break;
			}
		}
	}

	return tempidx;
}

double length(std::vector<double> a) {
	return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

double error(int n, double* x_new, double* x_old) {
	double sum = 0;
	for (int i = 0; i < n; i++) {
		sum += (x_new[i] - x_old[i])*(x_new[i] - x_old[i]);
	}

	return std::sqrt(sum);
}

void makeSurfaceMesh(std::vector<designVariable> const &x, std::vector<std::vector<double>> const &nodes) {

	std::ofstream SurfaceObjFile;

	SurfaceObjFile.open("mesh.obj");
	int nodeCounter = 0;

	for (int i = 0; i < x.size(); i++) {

		if (x[i].neighbors[0] == -1) { //add z+ indices and z+ face
			SurfaceObjFile << "v " << nodes[x[i].nodeIdx[3]][0] << " " << nodes[x[i].nodeIdx[3]][1] << " " << nodes[x[i].nodeIdx[3]][2] << std::endl;
			SurfaceObjFile << "v " << nodes[x[i].nodeIdx[6]][0] << " " << nodes[x[i].nodeIdx[6]][1] << " " << nodes[x[i].nodeIdx[6]][2] << std::endl;
			SurfaceObjFile << "v " << nodes[x[i].nodeIdx[7]][0] << " " << nodes[x[i].nodeIdx[7]][1] << " " << nodes[x[i].nodeIdx[7]][2] << std::endl;
			SurfaceObjFile << "v " << nodes[x[i].nodeIdx[5]][0] << " " << nodes[x[i].nodeIdx[5]][1] << " " << nodes[x[i].nodeIdx[5]][2] << std::endl;
			SurfaceObjFile << "f " << nodeCounter + 1 << " " << nodeCounter + 2 << " " << nodeCounter + 3 << std::endl;
			SurfaceObjFile << "f " << nodeCounter + 3 << " " << nodeCounter + 4 << " " << nodeCounter + 1 << std::endl;
			nodeCounter += 4;
		}

		if (x[i].neighbors[1] == -1) { //add z- indices and z- face		}
			SurfaceObjFile << "v " << nodes[x[i].nodeIdx[0]][0] << " " << nodes[x[i].nodeIdx[0]][1] << " " << nodes[x[i].nodeIdx[0]][2] << std::endl;
			SurfaceObjFile << "v " << nodes[x[i].nodeIdx[1]][0] << " " << nodes[x[i].nodeIdx[1]][1] << " " << nodes[x[i].nodeIdx[1]][2] << std::endl;
			SurfaceObjFile << "v " << nodes[x[i].nodeIdx[4]][0] << " " << nodes[x[i].nodeIdx[4]][1] << " " << nodes[x[i].nodeIdx[4]][2] << std::endl;
			SurfaceObjFile << "v " << nodes[x[i].nodeIdx[2]][0] << " " << nodes[x[i].nodeIdx[2]][1] << " " << nodes[x[i].nodeIdx[2]][2] << std::endl;
			SurfaceObjFile << "f " << nodeCounter + 1 << " " << nodeCounter + 2 << " " << nodeCounter + 3 << std::endl;
			SurfaceObjFile << "f " << nodeCounter + 3 << " " << nodeCounter + 4 << " " << nodeCounter + 1 << std::endl;
			nodeCounter += 4;
		}
		if (x[i].neighbors[2] == -1) { //add y+ indices and y+ face
			SurfaceObjFile << "v " << nodes[x[i].nodeIdx[2]][0] << " " << nodes[x[i].nodeIdx[2]][1] << " " << nodes[x[i].nodeIdx[2]][2] << std::endl;
			SurfaceObjFile << "v " << nodes[x[i].nodeIdx[4]][0] << " " << nodes[x[i].nodeIdx[4]][1] << " " << nodes[x[i].nodeIdx[4]][2] << std::endl;
			SurfaceObjFile << "v " << nodes[x[i].nodeIdx[7]][0] << " " << nodes[x[i].nodeIdx[7]][1] << " " << nodes[x[i].nodeIdx[7]][2] << std::endl;
			SurfaceObjFile << "v " << nodes[x[i].nodeIdx[5]][0] << " " << nodes[x[i].nodeIdx[5]][1] << " " << nodes[x[i].nodeIdx[5]][2] << std::endl;
			SurfaceObjFile << "f " << nodeCounter + 1 << " " << nodeCounter + 2 << " " << nodeCounter + 3 << std::endl;
			SurfaceObjFile << "f " << nodeCounter + 3 << " " << nodeCounter + 4 << " " << nodeCounter + 1 << std::endl;
			nodeCounter += 4;
		}

		if (x[i].neighbors[3] == -1) { //add y- indices and y- face
			SurfaceObjFile << "v " << nodes[x[i].nodeIdx[0]][0] << " " << nodes[x[i].nodeIdx[0]][1] << " " << nodes[x[i].nodeIdx[0]][2] << std::endl;
			SurfaceObjFile << "v " << nodes[x[i].nodeIdx[1]][0] << " " << nodes[x[i].nodeIdx[1]][1] << " " << nodes[x[i].nodeIdx[1]][2] << std::endl;
			SurfaceObjFile << "v " << nodes[x[i].nodeIdx[6]][0] << " " << nodes[x[i].nodeIdx[6]][1] << " " << nodes[x[i].nodeIdx[6]][2] << std::endl;
			SurfaceObjFile << "v " << nodes[x[i].nodeIdx[3]][0] << " " << nodes[x[i].nodeIdx[3]][1] << " " << nodes[x[i].nodeIdx[3]][2] << std::endl;
			SurfaceObjFile << "f " << nodeCounter + 1 << " " << nodeCounter + 2 << " " << nodeCounter + 3 << std::endl;
			SurfaceObjFile << "f " << nodeCounter + 3 << " " << nodeCounter + 4 << " " << nodeCounter + 1 << std::endl;			nodeCounter += 4;
		}
		if (x[i].neighbors[4] == -1) { //add x+ indices and x+ face
			SurfaceObjFile << "v " << nodes[x[i].nodeIdx[1]][0] << " " << nodes[x[i].nodeIdx[1]][1] << " " << nodes[x[i].nodeIdx[1]][2] << std::endl;
			SurfaceObjFile << "v " << nodes[x[i].nodeIdx[6]][0] << " " << nodes[x[i].nodeIdx[6]][1] << " " << nodes[x[i].nodeIdx[6]][2] << std::endl;
			SurfaceObjFile << "v " << nodes[x[i].nodeIdx[7]][0] << " " << nodes[x[i].nodeIdx[7]][1] << " " << nodes[x[i].nodeIdx[7]][2] << std::endl;
			SurfaceObjFile << "v " << nodes[x[i].nodeIdx[4]][0] << " " << nodes[x[i].nodeIdx[4]][1] << " " << nodes[x[i].nodeIdx[4]][2] << std::endl;
			SurfaceObjFile << "f " << nodeCounter + 1 << " " << nodeCounter + 2 << " " << nodeCounter + 3 << std::endl;
			SurfaceObjFile << "f " << nodeCounter + 3 << " " << nodeCounter + 4 << " " << nodeCounter + 1 << std::endl;
			nodeCounter += 4;
		}
		if (x[i].neighbors[5] == -1) { //add x- indices and x- face
			SurfaceObjFile << "v " << nodes[x[i].nodeIdx[0]][0] << " " << nodes[x[i].nodeIdx[0]][1] << " " << nodes[x[i].nodeIdx[0]][2] << std::endl;
			SurfaceObjFile << "v " << nodes[x[i].nodeIdx[3]][0] << " " << nodes[x[i].nodeIdx[3]][1] << " " << nodes[x[i].nodeIdx[3]][2] << std::endl;
			SurfaceObjFile << "v " << nodes[x[i].nodeIdx[5]][0] << " " << nodes[x[i].nodeIdx[5]][1] << " " << nodes[x[i].nodeIdx[5]][2] << std::endl;
			SurfaceObjFile << "v " << nodes[x[i].nodeIdx[2]][0] << " " << nodes[x[i].nodeIdx[2]][1] << " " << nodes[x[i].nodeIdx[2]][2] << std::endl;
			SurfaceObjFile << "f " << nodeCounter + 1 << " " << nodeCounter + 2 << " " << nodeCounter + 3 << std::endl;
			SurfaceObjFile << "f " << nodeCounter + 3 << " " << nodeCounter + 4 << " " << nodeCounter + 1 << std::endl;
			nodeCounter += 4;
		}
	}
}

void makeMesh(std::vector<designVariable> const &x, std::vector<std::vector<double>> const &nodes) {

	std::ofstream SurfaceObjFile;

	SurfaceObjFile.open("finalObject.obj");
	int nodeCounter = 0;

	for (int i = 0; i < x.size(); i++) {

		if (x[i].value <= 0.001) {
			continue;
		}

		SurfaceObjFile << "v " << nodes[x[i].nodeIdx[3]][0] << " " << nodes[x[i].nodeIdx[3]][1] << " " << nodes[x[i].nodeIdx[3]][2] << std::endl;
		SurfaceObjFile << "v " << nodes[x[i].nodeIdx[6]][0] << " " << nodes[x[i].nodeIdx[6]][1] << " " << nodes[x[i].nodeIdx[6]][2] << std::endl;
		SurfaceObjFile << "v " << nodes[x[i].nodeIdx[7]][0] << " " << nodes[x[i].nodeIdx[7]][1] << " " << nodes[x[i].nodeIdx[7]][2] << std::endl;
		SurfaceObjFile << "v " << nodes[x[i].nodeIdx[5]][0] << " " << nodes[x[i].nodeIdx[5]][1] << " " << nodes[x[i].nodeIdx[5]][2] << std::endl;
		SurfaceObjFile << "f " << nodeCounter + 1 << " " << nodeCounter + 2 << " " << nodeCounter + 3 << std::endl;
		SurfaceObjFile << "f " << nodeCounter + 3 << " " << nodeCounter + 4 << " " << nodeCounter + 1 << std::endl;
		nodeCounter += 4;


		SurfaceObjFile << "v " << nodes[x[i].nodeIdx[0]][0] << " " << nodes[x[i].nodeIdx[0]][1] << " " << nodes[x[i].nodeIdx[0]][2] << std::endl;
		SurfaceObjFile << "v " << nodes[x[i].nodeIdx[1]][0] << " " << nodes[x[i].nodeIdx[1]][1] << " " << nodes[x[i].nodeIdx[1]][2] << std::endl;
		SurfaceObjFile << "v " << nodes[x[i].nodeIdx[4]][0] << " " << nodes[x[i].nodeIdx[4]][1] << " " << nodes[x[i].nodeIdx[4]][2] << std::endl;
		SurfaceObjFile << "v " << nodes[x[i].nodeIdx[2]][0] << " " << nodes[x[i].nodeIdx[2]][1] << " " << nodes[x[i].nodeIdx[2]][2] << std::endl;
		SurfaceObjFile << "f " << nodeCounter + 1 << " " << nodeCounter + 2 << " " << nodeCounter + 3 << std::endl;
		SurfaceObjFile << "f " << nodeCounter + 3 << " " << nodeCounter + 4 << " " << nodeCounter + 1 << std::endl;
		nodeCounter += 4;

		SurfaceObjFile << "v " << nodes[x[i].nodeIdx[2]][0] << " " << nodes[x[i].nodeIdx[2]][1] << " " << nodes[x[i].nodeIdx[2]][2] << std::endl;
		SurfaceObjFile << "v " << nodes[x[i].nodeIdx[4]][0] << " " << nodes[x[i].nodeIdx[4]][1] << " " << nodes[x[i].nodeIdx[4]][2] << std::endl;
		SurfaceObjFile << "v " << nodes[x[i].nodeIdx[7]][0] << " " << nodes[x[i].nodeIdx[7]][1] << " " << nodes[x[i].nodeIdx[7]][2] << std::endl;
		SurfaceObjFile << "v " << nodes[x[i].nodeIdx[5]][0] << " " << nodes[x[i].nodeIdx[5]][1] << " " << nodes[x[i].nodeIdx[5]][2] << std::endl;
		SurfaceObjFile << "f " << nodeCounter + 1 << " " << nodeCounter + 2 << " " << nodeCounter + 3 << std::endl;
		SurfaceObjFile << "f " << nodeCounter + 3 << " " << nodeCounter + 4 << " " << nodeCounter + 1 << std::endl;
		nodeCounter += 4;


		SurfaceObjFile << "v " << nodes[x[i].nodeIdx[0]][0] << " " << nodes[x[i].nodeIdx[0]][1] << " " << nodes[x[i].nodeIdx[0]][2] << std::endl;
		SurfaceObjFile << "v " << nodes[x[i].nodeIdx[1]][0] << " " << nodes[x[i].nodeIdx[1]][1] << " " << nodes[x[i].nodeIdx[1]][2] << std::endl;
		SurfaceObjFile << "v " << nodes[x[i].nodeIdx[6]][0] << " " << nodes[x[i].nodeIdx[6]][1] << " " << nodes[x[i].nodeIdx[6]][2] << std::endl;
		SurfaceObjFile << "v " << nodes[x[i].nodeIdx[3]][0] << " " << nodes[x[i].nodeIdx[3]][1] << " " << nodes[x[i].nodeIdx[3]][2] << std::endl;
		SurfaceObjFile << "f " << nodeCounter + 1 << " " << nodeCounter + 2 << " " << nodeCounter + 3 << std::endl;
		SurfaceObjFile << "f " << nodeCounter + 3 << " " << nodeCounter + 4 << " " << nodeCounter + 1 << std::endl;			nodeCounter += 4;

		SurfaceObjFile << "v " << nodes[x[i].nodeIdx[1]][0] << " " << nodes[x[i].nodeIdx[1]][1] << " " << nodes[x[i].nodeIdx[1]][2] << std::endl;
		SurfaceObjFile << "v " << nodes[x[i].nodeIdx[6]][0] << " " << nodes[x[i].nodeIdx[6]][1] << " " << nodes[x[i].nodeIdx[6]][2] << std::endl;
		SurfaceObjFile << "v " << nodes[x[i].nodeIdx[7]][0] << " " << nodes[x[i].nodeIdx[7]][1] << " " << nodes[x[i].nodeIdx[7]][2] << std::endl;
		SurfaceObjFile << "v " << nodes[x[i].nodeIdx[4]][0] << " " << nodes[x[i].nodeIdx[4]][1] << " " << nodes[x[i].nodeIdx[4]][2] << std::endl;
		SurfaceObjFile << "f " << nodeCounter + 1 << " " << nodeCounter + 2 << " " << nodeCounter + 3 << std::endl;
		SurfaceObjFile << "f " << nodeCounter + 3 << " " << nodeCounter + 4 << " " << nodeCounter + 1 << std::endl;
		nodeCounter += 4;

		SurfaceObjFile << "v " << nodes[x[i].nodeIdx[0]][0] << " " << nodes[x[i].nodeIdx[0]][1] << " " << nodes[x[i].nodeIdx[0]][2] << std::endl;
		SurfaceObjFile << "v " << nodes[x[i].nodeIdx[3]][0] << " " << nodes[x[i].nodeIdx[3]][1] << " " << nodes[x[i].nodeIdx[3]][2] << std::endl;
		SurfaceObjFile << "v " << nodes[x[i].nodeIdx[5]][0] << " " << nodes[x[i].nodeIdx[5]][1] << " " << nodes[x[i].nodeIdx[5]][2] << std::endl;
		SurfaceObjFile << "v " << nodes[x[i].nodeIdx[2]][0] << " " << nodes[x[i].nodeIdx[2]][1] << " " << nodes[x[i].nodeIdx[2]][2] << std::endl;
		SurfaceObjFile << "f " << nodeCounter + 1 << " " << nodeCounter + 2 << " " << nodeCounter + 3 << std::endl;
		SurfaceObjFile << "f " << nodeCounter + 3 << " " << nodeCounter + 4 << " " << nodeCounter + 1 << std::endl;
		nodeCounter += 4;

	}
}

void createVegFile(double E, double ni, double density, std::vector<std::vector<double>> const &nodes, std::vector<designVariable> const &x) {

	std::ofstream myfile;
	myfile.open("points.veg");
	myfile << "*VERTICES\n";
	myfile << nodes.size() << " 3 0 0\n";
	for (int i = 0; i < nodes.size(); i++) {
		//	std::cout << i << std::endl;
		myfile << i + 1 << " " << nodes[i][0] << " " << nodes[i][1] << " " << nodes[i][2] << std::endl;
	}
	myfile << "*ELEMENTS\nCUBIC\n";
	myfile << x.size() << " 8 0 0\n";
	for (int i = 0; i < x.size(); i++) {
		myfile << i + 1 << " " << x[i].nodeIdx[0] + 1 << " " << x[i].nodeIdx[1] + 1 << " " << x[i].nodeIdx[4] + 1 << " " << x[i].nodeIdx[2] + 1 << " " << x[i].nodeIdx[3] + 1 << " " << x[i].nodeIdx[6] + 1 << " " << x[i].nodeIdx[7] + 1 << " " << x[i].nodeIdx[5] + 1 << std::endl;
	}

	//get material for each element
	for (int i = 0; i < x.size(); i++) {
		myfile << "*MATERIAL mat" << i + 1 << std::endl;
		myfile << "ENU, " << density << ", " << std::pow(x[i].value, 3)*E << ", " << ni << std::endl;
	}

	for (int i = 0; i < x.size(); i++) {
		myfile << "*SET set" << i + 1 << std::endl;
		myfile << x[i].varIdx + 1 << std::endl;
	}

	myfile << "*REGION" << std::endl;
	for (int i = 0; i < x.size(); i++) {
		myfile << "set" << i + 1 << ", mat" << i + 1 << std::endl;
	}

}

void calc_rho(std::vector<std::vector<double>> const &nodes, std::vector<designVariable>  &x, int NeighbourLayers, double dx, int numOfVoxelsX, int numOfVoxelsY, int numOfVoxelsZ) {
	//calculate rho for all x;
	double r0 = dx * NeighbourLayers*1.5;
	double wjx = 0;
	double sum_wj = 0;
	double wj = 0;
	int nx = numOfVoxelsX + 1;
	int ny = numOfVoxelsY + 1;
	int nz = numOfVoxelsZ + 1;
	int startVoxelZ;
	int startVoxelTemp;
	std::set<int> idxSet;

	int tempidx;
	for (int l = 0; l < x.size(); l++) {

		if (x[l].tunable == false) {
			continue;
		}

		int curridx = x[l].varIdx;

		for (int z = -NeighbourLayers; z <= NeighbourLayers; z++) {
			for (int y = -NeighbourLayers; y <= NeighbourLayers; y++) {
				for (int xAxis = -NeighbourLayers; xAxis <= NeighbourLayers; xAxis++) {
					if (xAxis == 0 && y == 0 && z == 0) { continue; }
					tempidx = moveZ(x, curridx, z);
					if (tempidx == -1) {
						continue;
					}
					tempidx = moveY(x, tempidx, y);
					if (tempidx == -1) {
						continue;
					}
					tempidx = moveX(x, tempidx, xAxis);
					if (tempidx == -1) {
						continue;
					}

					idxSet.insert(tempidx);
				}
			}
		}

		//calc  rho 
		for (auto i = idxSet.begin(); i != idxSet.end(); i++) {
			designVariable temp = x[*i];
			wj = (length(temp.position) - r0) / r0;
			if (wj > 0) {
				wjx += wj * temp.value;
				sum_wj += wj;
			}
		}

		x[l].rho = wjx / sum_wj;

		idxSet.clear();
	}
}

void calc_sensitivty_filter(double* dc_filtered, double* dc, std::vector<std::vector<double>> const &nodes, std::vector<designVariable>  &x, int NeighbourLayers, double dx, int numOfVoxelsX, int numOfVoxelsY, int numOfVoxelsZ) {
	//calculate rho for all x;
	double r0 = dx * NeighbourLayers*1.5;
	double wjx = 0;
	double sum_wj = 0;
	double wj = 0;
	int nx = numOfVoxelsX + 1;
	int ny = numOfVoxelsY + 1;
	int nz = numOfVoxelsZ + 1;
	int startVoxelZ;
	int startVoxelTemp;
	std::set<int> idxSet;

	int tempidx;
	for (int l = 0; l < x.size(); l++) {

		if (x[l].tunable == false) {
			continue;
		}

		int curridx = x[l].varIdx;

		for (int z = -NeighbourLayers; z <= NeighbourLayers; z++) {
			for (int y = -NeighbourLayers; y <= NeighbourLayers; y++) {
				for (int xAxis = -NeighbourLayers; xAxis <= NeighbourLayers; xAxis++) {
					if (xAxis == 0 && y == 0 && z == 0) { continue; }
					tempidx = moveZ(x, curridx, z);
					if (tempidx == -1) {
						continue;
					}
					tempidx = moveY(x, tempidx, y);
					if (tempidx == -1) {
						continue;
					}
					tempidx = moveX(x, tempidx, xAxis);
					if (tempidx == -1) {
						continue;
					}

					idxSet.insert(tempidx);
				}
			}
		}

		//calc  rho 
		for (auto i = idxSet.begin(); i != idxSet.end(); i++) {
			designVariable temp = x[*i];
			wj = (length(temp.position) - r0) / r0;
			if (wj > 0) {
				wjx += wj * temp.value*dc[temp.varIdx];
				sum_wj += wj;
			}
		}

		dc_filtered[l] = (1 / (x[l].value*sum_wj))*wjx;

		idxSet.clear();
		sum_wj = 0;
		wjx = 0;
	}
}


int main() {

	//voxel parameters
	int numOfVoxelsX = 30;
	int numOfVoxelsY = 60;
	int numOfVoxelsZ = 1;
	float dx = 1;//size - irrelevant?
	std::vector<designVariable> x;

	designVariable temp_x; //helper variable for voxel creation
	int idx = 0; //helper for valriable voxel creation
	//number of nodes
	int nx = numOfVoxelsX + 1;
	int ny = numOfVoxelsY + 1;
	int nz = numOfVoxelsZ + 1;
	std::vector<std::vector<double>> nodes; //node locations
	//create array that connects node indices to coordinates.
	nodes.resize(nx*ny*nz);
	//tested and seems to work
	std::cout << "making voxels" << std::endl;
	for (int kk = 0; kk < numOfVoxelsZ; kk++) { //z axis
		for (int jj = 0; jj < numOfVoxelsY; jj++) { // y axis
			for (int ii = 0; ii < numOfVoxelsX; ii++) { //x axis
				nodes[kk * nx*ny + jj * nx + ii] = std::vector<double>{ ii*dx, jj*dx, kk*dx };
				nodes[kk * nx*ny + jj * nx + ii + 1] = std::vector<double>{ (ii + 1)*dx, jj*dx, kk*dx };
				nodes[kk * nx*ny + jj * nx + ii] = std::vector<double>{ ii*dx, jj*dx, kk*dx };
				nodes[kk * nx*ny + (jj + 1) * nx + ii] = std::vector<double>{ ii*dx, (jj + 1)*dx, kk*dx };
				nodes[(kk + 1) * nx*ny + jj * nx + ii] = std::vector<double>{ ii*dx, jj*dx, (kk + 1)*dx };
				nodes[kk * nx*ny + (jj + 1) * nx + ii + 1] = std::vector<double>{ (ii + 1)*dx, (jj + 1)*dx, kk*dx };
				nodes[(kk + 1) * nx*ny + (jj + 1) * nx + ii] = std::vector<double>{ ii*dx, (jj + 1)*dx, (kk + 1)*dx };
				nodes[(kk + 1) * nx*ny + jj * nx + ii + 1] = std::vector<double>{ (ii + 1)*dx, jj*dx, (kk + 1)*dx };
				nodes[(kk + 1) * nx*ny + (jj + 1) * nx + ii + 1] = std::vector<double>{ (ii + 1)*dx, (jj + 1)*dx, (kk + 1)*dx };
				//idx of neighbour voxels
				temp_x.neighbors[0] = kk == numOfVoxelsZ - 1 ? -1 : idx + (ny - 1)*(nx - 1); //z+
				temp_x.neighbors[1] = kk == 0 ? -1 : idx - (ny - 1) * (nx - 1); //z-
				temp_x.neighbors[2] = jj == numOfVoxelsY - 1 ? -1 : idx + nx - 1;
				temp_x.neighbors[3] = jj == 0 ? -1 : idx - (nx - 1);
				temp_x.neighbors[4] = ii == numOfVoxelsX - 1 ? -1 : idx + 1;//x+
				temp_x.neighbors[5] = ii == 0 ? -1 : idx - 1;//x-
				//position of voxel
				temp_x.position = std::vector<double>{ (ii + 0.5)*dx, (jj + 0.5)*dx, (kk + 0.5)*dx };
				//	std::cout << nodes[kk * nx*ny + jj * nx + ii][0]<<" "<< nodes[kk * nx*ny + jj * nx + ii][1] <<" "<< nodes[kk * nx*ny + jj * nx + ii][2] << std::endl;
					//std::cout << length(temp_x.position) << std::endl;
				temp_x.varIdx = idx;
				idx++;
				//voxel's node's indices
				temp_x.nodeIdx[0] = kk * nx*ny + jj * nx + ii;
				temp_x.nodeIdx[1] = kk * nx*ny + jj * nx + ii + 1;
				temp_x.nodeIdx[2] = kk * nx*ny + (jj + 1) * nx + ii;
				temp_x.nodeIdx[3] = (kk + 1) * nx*ny + jj * nx + ii;
				temp_x.nodeIdx[4] = kk * nx*ny + (jj + 1) * nx + ii + 1;
				temp_x.nodeIdx[5] = (kk + 1) * nx*ny + (jj + 1) * nx + ii;
				temp_x.nodeIdx[6] = (kk + 1) * nx*ny + jj * nx + ii + 1;
				temp_x.nodeIdx[7] = (kk + 1) * nx*ny + (jj + 1) * nx + ii + 1;

				x.push_back(temp_x);
			}
		}
	}

	std::cout << "done making voxels" << std::endl;

	int NeighbourLayers = 1; //amount of layers for smoothing

	//forces - get nodeidx and force vector
	int numOfForces = 1;//*ny;
	int numOfUntunableParams = 0;
	std::vector<std::pair<int, std::vector<double>>> force;
	force.resize(numOfForces * 2);
	std::cout << "num of forces: " << force.size() << std::endl;
	int curridx = numOfVoxelsX * numOfVoxelsY*numOfVoxelsZ - 1;
	force[0].first = x[curridx].nodeIdx[4]; force[0].second = { 1,0,0 };
	force[1].first = x[curridx].nodeIdx[7]; force[1].second = { 1,0,0 };
	for (int i = 1; i < numOfForces; i++) {
		numOfUntunableParams++;
		force[2 * i].first = x[curridx].nodeIdx[2]; //last voxel, corner node
		force[2 * i].second.resize(3); force[2 * i].second = { 1,0,0 };
		force[2 * i + 1].first = x[curridx].nodeIdx[5];
		force[2 * i + 1].second.resize(3); force[2 * i + 1].second = { 1,0,0 };
		//x[curridx].tunable = false;
		x[curridx].value = 0.5;
		curridx = moveX(x, curridx, -1);
	}

	//constarints - get nodexidx and add to set
	int numofConstraints = nx * ny;
	std::set<int> constraints;
	curridx = 0;
	for (int ix = 0; ix < numOfVoxelsX; ix++) {
		for (int iz = 0; iz < numOfVoxelsZ; iz++) {
			int temp = moveX(x, curridx, ix);
			temp = moveZ(x, temp, iz);
			//	x[temp].tunable = false;
			x[temp].value = 0.5;
			constraints.insert(x[temp].nodeIdx[0]);
			constraints.insert(x[temp].nodeIdx[1]);
			constraints.insert(x[temp].nodeIdx[6]);
			constraints.insert(x[temp].nodeIdx[3]);
		}
	}

	//map to free DOF for FEM

	std::vector<int> freeDOF;
	freeDOF.resize(nodes.size());
	int ind = 0; //holds number of free nodes
	for (int i = 0; i < freeDOF.size(); i++) {

		if (constraints.find(i) == constraints.end()) {

			freeDOF[i] = ind;
			ind++;
		}
		else {

			freeDOF[i] = -1;
		}

	}

	int numOfFixed = nodes.size() - ind; //number of fixed nodes

	//initialize convergence criterion variable with large value (force at least one iteration)
	double change = 1000;

	//constraints

	//base material properties
	double ni = 0.3;
	double E = 1;
	double density = 2800;

	//algorithm constants
	double epsilon = 0.0001;
	double V_allowed = 0.3*numOfVoxelsX*numOfVoxelsY*numOfVoxelsZ;
	std::cout << V_allowed << std::endl;
	int power = 3; //penalty 


	Eigen::MatrixXd KE(24, 24); //assumes x[0] value of node is 1 at the start	

	KE = K_mat();


	//visualization using libigl

	Eigen::MatrixXd V;
	Eigen::MatrixXi F;

	Eigen::MatrixXd p(force.size(), 3), p2(constraints.size(), 3);
	Eigen::RowVector3d c(1, 0, 0), c2(0, 1, 0);

	for (int jj = 0; jj < force.size(); jj++) {
		p(jj, 0) = nodes[force[jj].first][0]; p(jj, 1) = nodes[force[jj].first][1]; p(jj, 2) = nodes[force[jj].first][2];
	}

	int currIdx = 0;
	for (auto jj = constraints.begin(); jj != constraints.end(); jj++) {
		p2(currIdx, 0) = nodes[*jj][0]; p2(currIdx, 1) = nodes[*jj][1]; p2(currIdx, 2) = nodes[*jj][2];
		currIdx++;
	}


	//show constraints and forces
	// Plot the mesh
	makeSurfaceMesh(x, nodes);
	igl::readOBJ("mesh.obj", V, F);
	igl::opengl::glfw::Viewer viewer;
	viewer.data().add_points(p, c);
	viewer.data().add_points(p2, c2);
	viewer.data().set_mesh(V, F);
	viewer.launch();

	//initialize mma parameters
	double* dc = new double[x.size()];
	double* xnew = new double[x.size()];
	double* dg = new double[x.size()];
	double* g = new double[1];
	double *xmin = new double[x.size()];
	double *xmax = new double[x.size()];
	double* dc_filtered = new double[x.size()];

	MMASolver* mma = new MMASolver(x.size(), 1);

	g[0] = -V_allowed;
	for (int i = 0; i < x.size(); i++) {
		dg[i] = 1;
		xmin[i] = 0.001;
		xmax[i] = 1;
		g[0] += x[i].value;
		xnew[i] = x[i].value;
	}

	calc_rho(nodes, x, 1, dx, numOfVoxelsX, numOfVoxelsY, numOfVoxelsZ);

	while (change > 0.01) {

		Eigen::VectorXd f(3 * ind), u(nodes.size() * 3);

		for (int j = 0; j < 3 * ind; j++) {
			f(j) = 0;
		}
		for (int j = 0; j < numOfForces; j++) {
			int tempForceid = freeDOF[force[j].first] * 3;
			f(tempForceid + 0) = force[j].second[0];
			f(tempForceid + 1) = force[j].second[1];
			f(tempForceid + 2) = force[j].second[2];
		}

		//	std::cout << "done with forces! " << std::endl;

		auto K_global = build_global_stiffness_matrix(x, nodes.size(), KE, power, freeDOF, ind);

		auto u_free = do_FEM(K_global, f);

		//	std::cout << "calculate free u" << std::endl;

		//auto DenseK = K_global.toDense();

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

		std::cout << "success in simulation" << std::endl;

		std::cout << "calculating derivatives" << std::endl;
		for (int i = 0; i < x.size(); i++) {


			Eigen::VectorXd displacement(24), f_el(24);
			displacement(0) = u[x[i].nodeIdx[0] * 3]; displacement(1) = u[x[i].nodeIdx[0] * 3 + 1]; displacement(2) = u[x[i].nodeIdx[0] * 3 + 2];
			displacement(3) = u[x[i].nodeIdx[1] * 3]; displacement(4) = u[x[i].nodeIdx[1] * 3 + 1]; displacement(5) = u[x[i].nodeIdx[1] * 3 + 2];
			displacement(6) = u[x[i].nodeIdx[2] * 3]; displacement(7) = u[x[i].nodeIdx[2] * 3 + 1]; displacement(8) = u[x[i].nodeIdx[2] * 3 + 2];
			displacement(9) = u[x[i].nodeIdx[4] * 3]; displacement(10) = u[x[i].nodeIdx[4] * 3 + 1]; displacement(11) = u[x[i].nodeIdx[4] * 3 + 2];
			displacement(12) = u[x[i].nodeIdx[3] * 3]; displacement(13) = u[x[i].nodeIdx[3] * 3 + 1]; displacement(14) = u[x[i].nodeIdx[3] * 3 + 2];
			displacement(15) = u[x[i].nodeIdx[6] * 3]; displacement(16) = u[x[i].nodeIdx[6] * 3 + 1]; displacement(17) = u[x[i].nodeIdx[6] * 3 + 2];
			displacement(18) = u[x[i].nodeIdx[5] * 3]; displacement(19) = u[x[i].nodeIdx[5] * 3 + 1]; displacement(20) = u[x[i].nodeIdx[5] * 3 + 2];
			displacement(21) = u[x[i].nodeIdx[7] * 3]; displacement(22) = u[x[i].nodeIdx[7] * 3 + 1]; displacement(23) = u[x[i].nodeIdx[7] * 3 + 2];

			dc[i] = -power * std::pow(x[i].value, power - 1)* displacement.transpose()*KE*displacement;

			//	std::cout << "element " << i << " derivative: " << dc[i] << std::endl;
		}

		//calculate sensitivty filter
	//	calc_sensitivty_filter(dc_filtered, dc, nodes, x, 2, dx, numOfVoxelsX, numOfVoxelsY, numOfVoxelsZ);


		mma->Update(xnew, dc, g, dg, xmin, xmax);

		change = -1;
		g[0] = -V_allowed;
		for (int i = 0; i < x.size(); i++) {
			change = Max(change, std::abs(x[i].value - xnew[i]));
			x[i].value = xnew[i];
		}
		//	std::cout <<"element "<<i<<": "<< x[i].value << std::endl;

		calc_rho(nodes, x, 1, dx, numOfVoxelsX, numOfVoxelsY, numOfVoxelsZ);

		for (int i = 0; i < x.size(); i++) {

			g[0] += x[i].value;
		}

		std::cout << g[0] << std::endl;

	}

	delete[] xnew;
	delete[] dc;
	delete[] xmin;
	delete[] xmax;
	delete[] g;
	delete[] dg;

	makeMesh(x, nodes);

}

