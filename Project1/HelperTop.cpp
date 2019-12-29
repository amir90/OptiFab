#pragma once
#include "HelperTop.h"
#include<iostream>
#include<fstream>
#include<set>


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

double length(std::vector<double> a) {
	return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

int moveX(std::vector<designVariable> const &x, int curridx, int move) {

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

int turn_untunable(std::vector<std::vector<double>> const &nodes, std::vector<designVariable>  &x, int NeighbourLayers, double dx, std::vector<std::pair<int,std::vector<double>>>& force) {

	int numOfUntunable = 0;
	int startVoxelZ;
	int startVoxelTemp;
	std::set<int> forceIdxSet;
	std::vector<int> voxelIdx;

	int tempidx;
	for (auto l = force.begin(); l != force.end(); l++) { //for all forces
		forceIdxSet.insert(l->first);
	}

	std::cout << "inside untunable vars" << std::endl;


	for (int l = 0; l < x.size(); l++) { //for all voxels
		for (int j = 0; j < 8; j++) {
			if (forceIdxSet.find(x[l].nodeIdx[j]) != forceIdxSet.end()) { 
				voxelIdx.push_back(l);
				continue; //voxel known to be in vector - continue to next voxel
			}
		}
	}

	for (int l = 0; l < voxelIdx.size(); l++) {
		int curridx = voxelIdx[l];

		for (int z = -NeighbourLayers; z <= NeighbourLayers; z++) {
			for (int y = -NeighbourLayers; y <= NeighbourLayers; y++) {
				for (int xAxis = -NeighbourLayers; xAxis <= NeighbourLayers; xAxis++) {
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
					if (x[tempidx].tunable == true) {
						x[tempidx].tunable = false;
						x[tempidx].value = 1;
						x[tempidx].rho = 1;
						numOfUntunable++;
					}
				}
			}
		}
	}

	return numOfUntunable;
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
			SurfaceObjFile << "f " << nodeCounter + 3 << " " << nodeCounter + 2 << " " << nodeCounter + 1 << std::endl;
			SurfaceObjFile << "f " << nodeCounter + 1 << " " << nodeCounter + 4 << " " << nodeCounter + 3 << std::endl;
			nodeCounter += 4;

		}
		if (x[i].neighbors[2] == -1) { //add y+ indices and y+ face
			SurfaceObjFile << "v " << nodes[x[i].nodeIdx[2]][0] << " " << nodes[x[i].nodeIdx[2]][1] << " " << nodes[x[i].nodeIdx[2]][2] << std::endl;
			SurfaceObjFile << "v " << nodes[x[i].nodeIdx[4]][0] << " " << nodes[x[i].nodeIdx[4]][1] << " " << nodes[x[i].nodeIdx[4]][2] << std::endl;
			SurfaceObjFile << "v " << nodes[x[i].nodeIdx[7]][0] << " " << nodes[x[i].nodeIdx[7]][1] << " " << nodes[x[i].nodeIdx[7]][2] << std::endl;
			SurfaceObjFile << "v " << nodes[x[i].nodeIdx[5]][0] << " " << nodes[x[i].nodeIdx[5]][1] << " " << nodes[x[i].nodeIdx[5]][2] << std::endl;
			SurfaceObjFile << "f " << nodeCounter + 3 << " " << nodeCounter + 2 << " " << nodeCounter + 1 << std::endl;
			SurfaceObjFile << "f " << nodeCounter + 1 << " " << nodeCounter + 4 << " " << nodeCounter + 3 << std::endl;
			nodeCounter += 4;
		}

		if (x[i].neighbors[3] == -1) { //add y- indices and y- face
			SurfaceObjFile << "v " << nodes[x[i].nodeIdx[0]][0] << " " << nodes[x[i].nodeIdx[0]][1] << " " << nodes[x[i].nodeIdx[0]][2] << std::endl;
			SurfaceObjFile << "v " << nodes[x[i].nodeIdx[1]][0] << " " << nodes[x[i].nodeIdx[1]][1] << " " << nodes[x[i].nodeIdx[1]][2] << std::endl;
			SurfaceObjFile << "v " << nodes[x[i].nodeIdx[6]][0] << " " << nodes[x[i].nodeIdx[6]][1] << " " << nodes[x[i].nodeIdx[6]][2] << std::endl;
			SurfaceObjFile << "v " << nodes[x[i].nodeIdx[3]][0] << " " << nodes[x[i].nodeIdx[3]][1] << " " << nodes[x[i].nodeIdx[3]][2] << std::endl;
			SurfaceObjFile << "f " << nodeCounter + 1 << " " << nodeCounter + 2 << " " << nodeCounter + 3 << std::endl;
			SurfaceObjFile << "f " << nodeCounter + 3 << " " << nodeCounter + 4 << " " << nodeCounter + 1 << std::endl;
			nodeCounter += 4;
		}

		if (x[i].neighbors[4] == -1) { //add x+ indices and x+ face
			SurfaceObjFile << "v " << nodes[x[i].nodeIdx[1]][0] << " " << nodes[x[i].nodeIdx[1]][1] << " " << nodes[x[i].nodeIdx[1]][2] << std::endl;
			SurfaceObjFile << "v " << nodes[x[i].nodeIdx[6]][0] << " " << nodes[x[i].nodeIdx[6]][1] << " " << nodes[x[i].nodeIdx[6]][2] << std::endl;
			SurfaceObjFile << "v " << nodes[x[i].nodeIdx[7]][0] << " " << nodes[x[i].nodeIdx[7]][1] << " " << nodes[x[i].nodeIdx[7]][2] << std::endl;
			SurfaceObjFile << "v " << nodes[x[i].nodeIdx[4]][0] << " " << nodes[x[i].nodeIdx[4]][1] << " " << nodes[x[i].nodeIdx[4]][2] << std::endl;
			SurfaceObjFile << "f " << nodeCounter + 3 << " " << nodeCounter + 2 << " " << nodeCounter + 1 << std::endl;
			SurfaceObjFile << "f " << nodeCounter + 1 << " " << nodeCounter + 4 << " " << nodeCounter + 3 << std::endl;
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

	SurfaceObjFile.close();
}

void makeSurfaceMesh2(std::vector<designVariable> const &x, std::vector<std::vector<double>> const &nodes) {

	std::ofstream SurfaceObjFile2;

	SurfaceObjFile2.open("mesh2.obj");
	int nodeCounter = 0;

	//need to set up a node to surface node adapter array

	std::vector<int> surfaceNodeAdapter(nodes.size(), -1);
	int counter = 0;

	for (int i = 0; i < x.size(); i++) {

		std::vector<int> indices = { 3,6,7,5 };

		if (x[i].neighbors[0] == -1) {
			for (int j = 0; j < indices.size(); j++) {
				int k = indices[j];
				if (surfaceNodeAdapter[x[i].nodeIdx[k]] == -1) {
					SurfaceObjFile2 << "v " << nodes[x[i].nodeIdx[k]][0] << " " << nodes[x[i].nodeIdx[k]][1] << " " << nodes[x[i].nodeIdx[k]][2] << std::endl;
					surfaceNodeAdapter[x[i].nodeIdx[k]] = counter++;
				}
			}
		}

		indices = { 0,1,4,2 };

		if (x[i].neighbors[1] == -1) {
			for (int j = 0; j < indices.size(); j++) {
				int k = indices[j];
				if (surfaceNodeAdapter[x[i].nodeIdx[k]] == -1) {
					SurfaceObjFile2 << "v " << nodes[x[i].nodeIdx[k]][0] << " " << nodes[x[i].nodeIdx[k]][1] << " " << nodes[x[i].nodeIdx[k]][2] << std::endl;
					surfaceNodeAdapter[x[i].nodeIdx[k]] = counter++;
				}
			}
		}

		indices = { 2,4,7,5 };

		if (x[i].neighbors[2] == -1) {
			for (int j = 0; j < indices.size(); j++) {
				int k = indices[j];
				if (surfaceNodeAdapter[x[i].nodeIdx[k]] == -1) {
					SurfaceObjFile2 << "v " << nodes[x[i].nodeIdx[k]][0] << " " << nodes[x[i].nodeIdx[k]][1] << " " << nodes[x[i].nodeIdx[k]][2] << std::endl;
					surfaceNodeAdapter[x[i].nodeIdx[k]] = counter++;
				}
			}
		}

		indices = { 0,1,6,3 };

		if (x[i].neighbors[3] == -1) {
			for (int j = 0; j < indices.size(); j++) {
				int k = indices[j];
				if (surfaceNodeAdapter[x[i].nodeIdx[k]] == -1) {
					SurfaceObjFile2 << "v " << nodes[x[i].nodeIdx[k]][0] << " " << nodes[x[i].nodeIdx[k]][1] << " " << nodes[x[i].nodeIdx[k]][2] << std::endl;
					surfaceNodeAdapter[x[i].nodeIdx[k]] = counter++;
				}
			}
		}

		indices = { 1,6,7,4 };

		if (x[i].neighbors[4] == -1) {
			for (int j = 0; j < indices.size(); j++) {
				int k = indices[j];
				if (surfaceNodeAdapter[x[i].nodeIdx[k]] == -1) {
					SurfaceObjFile2 << "v " << nodes[x[i].nodeIdx[k]][0] << " " << nodes[x[i].nodeIdx[k]][1] << " " << nodes[x[i].nodeIdx[k]][2] << std::endl;
					surfaceNodeAdapter[x[i].nodeIdx[k]] = counter++;
				}
			}
		}

		indices = { 0,3,5,2 };


		if (x[i].neighbors[5] == -1) {
			for (int j = 0; j < indices.size(); j++) {
				int k = indices[j];
				if (surfaceNodeAdapter[x[i].nodeIdx[k]] == -1) {
					SurfaceObjFile2 << "v " << nodes[x[i].nodeIdx[k]][0] << " " << nodes[x[i].nodeIdx[k]][1] << " " << nodes[x[i].nodeIdx[k]][2] << std::endl;
					surfaceNodeAdapter[x[i].nodeIdx[k]] = counter++;
				}
			}
		}
	}

	for (int i = 0; i < x.size(); i++) {
		//{ 3,6,7,5 }
		//SurfaceObjFile << "f " << nodeCounter + 1 << " " << nodeCounter + 2 << " " << nodeCounter + 3 << std::endl;
		//SurfaceObjFile << "f " << nodeCounter + 3 << " " << nodeCounter + 4 << " " << nodeCounter + 1 << std::endl;
		if (x[i].neighbors[0] == -1) { //add z+ indices and z+ face
			SurfaceObjFile2 << "f " << surfaceNodeAdapter[x[i].nodeIdx[3]] + 1 << " " << surfaceNodeAdapter[x[i].nodeIdx[6]] + 1 << " " << surfaceNodeAdapter[x[i].nodeIdx[7]] + 1 << std::endl;
			SurfaceObjFile2 << "f " << surfaceNodeAdapter[x[i].nodeIdx[7]] + 1 << " " << surfaceNodeAdapter[x[i].nodeIdx[5]] + 1 << " " << surfaceNodeAdapter[x[i].nodeIdx[3]] + 1 << std::endl;
		}
		//{0,1,4,2}
		//SurfaceObjFile << "f " << nodeCounter + 3 << " " << nodeCounter + 2 << " " << nodeCounter + 1 << std::endl;
		//SurfaceObjFile << "f " << nodeCounter + 1 << " " << nodeCounter + 4 << " " << nodeCounter + 3 << std::endl;
		if (x[i].neighbors[1] == -1) { //add z- indices and z- face		}
			SurfaceObjFile2 << "f " << surfaceNodeAdapter[x[i].nodeIdx[4]] + 1 << " " << surfaceNodeAdapter[x[i].nodeIdx[1]] + 1 << " " << surfaceNodeAdapter[x[i].nodeIdx[0]] + 1 << std::endl;
			SurfaceObjFile2 << "f " << surfaceNodeAdapter[x[i].nodeIdx[0]] + 1 << " " << surfaceNodeAdapter[x[i].nodeIdx[2]] + 1 << " " << surfaceNodeAdapter[x[i].nodeIdx[4]] + 1 << std::endl;

		}
		//{2,4,7,5}
		//SurfaceObjFile << "f " << nodeCounter + 1 << " " << nodeCounter + 2 << " " << nodeCounter + 3 << std::endl;
		//SurfaceObjFile << "f " << nodeCounter + 3 << " " << nodeCounter + 4 << " " << nodeCounter + 1 << std::endl;
		if (x[i].neighbors[2] == -1) { //add y+ indices and y+ face
			SurfaceObjFile2 << "f " << surfaceNodeAdapter[x[i].nodeIdx[7]] + 1 << " " << surfaceNodeAdapter[x[i].nodeIdx[4]] + 1 << " " << surfaceNodeAdapter[x[i].nodeIdx[2]] + 1 << std::endl;
			SurfaceObjFile2 << "f " << surfaceNodeAdapter[x[i].nodeIdx[2]] + 1 << " " << surfaceNodeAdapter[x[i].nodeIdx[5]] + 1 << " " << surfaceNodeAdapter[x[i].nodeIdx[7]] + 1 << std::endl;
		}
		//{0,1,6,3}
		if (x[i].neighbors[3] == -1) { //add y- indices and y- face
			SurfaceObjFile2 << "f " << surfaceNodeAdapter[x[i].nodeIdx[0]] + 1 << " " << surfaceNodeAdapter[x[i].nodeIdx[1]] + 1 << " " << surfaceNodeAdapter[x[i].nodeIdx[6]] + 1 << std::endl;
			SurfaceObjFile2 << "f " << surfaceNodeAdapter[x[i].nodeIdx[6]] + 1 << " " << surfaceNodeAdapter[x[i].nodeIdx[3]] + 1 << " " << surfaceNodeAdapter[x[i].nodeIdx[0]] + 1 << std::endl;
		}
		//{1,6,7,4}
		if (x[i].neighbors[4] == -1) { //add x+ indices and x+ face
			SurfaceObjFile2 << "f " << surfaceNodeAdapter[x[i].nodeIdx[7]] + 1 << " " << surfaceNodeAdapter[x[i].nodeIdx[6]] + 1 << " " << surfaceNodeAdapter[x[i].nodeIdx[1]] + 1 << std::endl;
			SurfaceObjFile2 << "f " << surfaceNodeAdapter[x[i].nodeIdx[1]] + 1 << " " << surfaceNodeAdapter[x[i].nodeIdx[4]] + 1 << " " << surfaceNodeAdapter[x[i].nodeIdx[7]] + 1 << std::endl;
		}
		//{0,3,5,2}
		if (x[i].neighbors[5] == -1) { //add x- indices and x- face
			SurfaceObjFile2 << "f " << surfaceNodeAdapter[x[i].nodeIdx[0]] + 1 << " " << surfaceNodeAdapter[x[i].nodeIdx[3]] + 1 << " " << surfaceNodeAdapter[x[i].nodeIdx[5]] + 1 << std::endl;
			SurfaceObjFile2 << "f " << surfaceNodeAdapter[x[i].nodeIdx[5]] + 1 << " " << surfaceNodeAdapter[x[i].nodeIdx[2]] + 1 << " " << surfaceNodeAdapter[x[i].nodeIdx[0]] + 1 << std::endl;
		}
	}
	SurfaceObjFile2.close();
}

void makeSurfaceMesh3(std::vector<designVariable> const &x, std::vector<std::vector<double>> const &nodes) {
	//writes mesh to obj file
	//inputs:
	//	x - vector of mesh cubes
	// nodes - vector of position of nodes

	std::ofstream SurfaceObjFile3;
	SurfaceObjFile3.open("mesh.obj");
	int nodeCounter = 0;

	//need to set up a node to surface node adapter array

	std::vector<int> surfaceNodeAdapter(nodes.size(), -1);
	int counter = 0;

	for (int i = 0; i < nodes.size(); i++) {
		SurfaceObjFile3 << "v " << nodes[i][0] << " " << nodes[i][1] << " " << nodes[i][2] << std::endl;

	}


	for (int i = 0; i < x.size(); i++) {
		//{ 3,6,7,5 }
		//SurfaceObjFile << "f " << nodeCounter + 1 << " " << nodeCounter + 2 << " " << nodeCounter + 3 << std::endl;
		//SurfaceObjFile << "f " << nodeCounter + 3 << " " << nodeCounter + 4 << " " << nodeCounter + 1 << std::endl;
		SurfaceObjFile3 << "f " << x[i].nodeIdx[3] + 1 << " " << x[i].nodeIdx[6] + 1 << " " << x[i].nodeIdx[7] + 1 << std::endl;
		SurfaceObjFile3 << "f " << x[i].nodeIdx[7] + 1 << " " << x[i].nodeIdx[5] + 1 << " " << x[i].nodeIdx[3] + 1 << std::endl;

		//{0,1,4,2}
		//SurfaceObjFile << "f " << nodeCounter + 3 << " " << nodeCounter + 2 << " " << nodeCounter + 1 << std::endl;
		//SurfaceObjFile << "f " << nodeCounter + 1 << " " << nodeCounter + 4 << " " << nodeCounter + 3 << std::endl;
		SurfaceObjFile3 << "f " << x[i].nodeIdx[4] + 1 << " " << x[i].nodeIdx[1] + 1 << " " << x[i].nodeIdx[0] + 1 << std::endl;
		SurfaceObjFile3 << "f " << x[i].nodeIdx[0] + 1 << " " << x[i].nodeIdx[2] + 1 << " " << x[i].nodeIdx[4] + 1 << std::endl;


		//{2,4,7,5}
		//SurfaceObjFile << "f " << nodeCounter + 1 << " " << nodeCounter + 2 << " " << nodeCounter + 3 << std::endl;
		//SurfaceObjFile << "f " << nodeCounter + 3 << " " << nodeCounter + 4 << " " << nodeCounter + 1 << std::endl;
		SurfaceObjFile3 << "f " << x[i].nodeIdx[7] + 1 << " " << x[i].nodeIdx[4] + 1 << " " << x[i].nodeIdx[2] + 1 << std::endl;
		SurfaceObjFile3 << "f " << x[i].nodeIdx[2] + 1 << " " << x[i].nodeIdx[5] + 1 << " " << x[i].nodeIdx[7] + 1 << std::endl;

		//{0,1,6,3}
		SurfaceObjFile3 << "f " << x[i].nodeIdx[0] + 1 << " " << x[i].nodeIdx[1] + 1 << " " << x[i].nodeIdx[6] + 1 << std::endl;
		SurfaceObjFile3 << "f " << x[i].nodeIdx[6] + 1 << " " << x[i].nodeIdx[3] + 1 << " " << x[i].nodeIdx[0] + 1 << std::endl;

		//{1,6,7,4}
		SurfaceObjFile3 << "f " << x[i].nodeIdx[7] + 1 << " " << x[i].nodeIdx[6] + 1 << " " << x[i].nodeIdx[1] + 1 << std::endl;
		SurfaceObjFile3 << "f " << x[i].nodeIdx[1] + 1 << " " << x[i].nodeIdx[4] + 1 << " " << x[i].nodeIdx[7] + 1 << std::endl;

		//{0,3,5,2}
		SurfaceObjFile3 << "f " << x[i].nodeIdx[0] + 1 << " " << x[i].nodeIdx[3] + 1 << " " << x[i].nodeIdx[5] + 1 << std::endl;
		SurfaceObjFile3 << "f " << x[i].nodeIdx[5] + 1 << " " << x[i].nodeIdx[2] + 1 << " " << x[i].nodeIdx[0] + 1 << std::endl;

	}
	SurfaceObjFile3.close();
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

void calc_rho(std::vector<std::vector<double>> const &nodes, std::vector<designVariable>  &x, int NeighbourLayers, double dx, int nVx, int nVy, int nVz, float r0) {
	//calculate rho for all x;

	for (int l = 0; l < x.size(); l++) { //for all voxels
		x[l].rho = x[l].value;
	}
	return;

		std::cout << "begin calculating rho" << std::endl;
		double wjx = 0;
		double sum_wj = 0;
		double wj = 0;
		int numOfVoxels = nVx * nVy*nVz;
		int nx = nVx + 1;
		int ny = nVy + 1;
		int nz = nVz + 1;
		int startVoxelZ;
		int startVoxelTemp;
		std::set<int> idxSet;

		int tempidx;
		for (int l = 0; l < x.size(); l++) { //for all voxels

			if (x[l].tunable==false) { continue; }

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

			if (l == 0) {

				std::cout <<" idxSet size for first voxel: "<< idxSet.size() << std::endl;
				for (auto k = idxSet.begin(); k != idxSet.end(); k++) {
					std::cout << " idxSet for first voxel: " << *k << std::endl;
				}
			}

			//calc  rho
			for (auto i = idxSet.begin(); i != idxSet.end(); i++) {
	//			std::cout << "element " << l << " node: " << *i << std::endl;
				designVariable temp = x[*i];
				wj = (length({ temp.position[0] - x[l].position[0],temp.position[1] - x[l].position[1],temp.position[2] - x[l].position[2] }) - r0) / r0;
				if (l == 0) {
					std::cout << "r0: " << r0 << std::endl;
					std::cout << "wj: " << wj<< std::endl;
				}
				if (wj < 0) {
					x[*i].influencesVoxels.insert(l);
					wjx += wj * temp.value;
				sum_wj += wj;
				}
			}
			x[l].sum_wj = std::abs(sum_wj);
		//std::cout << wjx/sum_wj << std::endl;
			x[l].rho = std::abs(wjx / sum_wj);

			idxSet.clear();
			sum_wj = 0;
			wjx = 0;
		}
		std::cout << "finished calculating rho" << std::endl;

}

void create_mesh(std::vector<std::vector<double>>& nodes, std::vector<designVariable> &x, int numOfVoxelsX, int numOfVoxelsY, int numOfVoxelsZ, float dx) {

	nodes.clear();
	x.clear();

	int idx = 0; //helper for variable voxel creation
	//number of nodes
	int nx = numOfVoxelsX + 1;
	int ny = numOfVoxelsY + 1;
	int nz = numOfVoxelsZ + 1;

	nodes.resize(nx*ny*nz);

	designVariable temp_x; //helper variable for voxel creation

	//voxel initialization - TODO: test on simple configuration
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
				temp_x.neighbors[2] = jj == numOfVoxelsY - 1 ? -1 : idx + nx - 1; //y+
				temp_x.neighbors[3] = jj == 0 ? -1 : idx - (nx - 1); //y-
				temp_x.neighbors[4] = ii == numOfVoxelsX - 1 ? -1 : idx + 1;//x+
				temp_x.neighbors[5] = ii == 0 ? -1 : idx - 1;//x-
				//position of voxel
				temp_x.position = std::vector<double>{ (ii + 0.5)*dx, (jj + 0.5)*dx, (kk + 0.5)*dx };
				//	std::cout << nodes[kk * nx*ny + jj * nx + ii][0]<<" "<< nodes[kk * nx*ny + jj * nx + ii][1] <<" "<< nodes[kk * nx*ny + jj * nx + ii][2] << std::endl;
					//std::cout << length(temp_x.position) << std::endl;
				temp_x.varIdx = idx;
				idx++;
				//voxel's node's indices
				temp_x.nodeIdx[0] = kk * nx*ny + jj * nx + ii; //local node 0
				temp_x.nodeIdx[1] = kk * nx*ny + jj * nx + ii + 1; //local node 1
				temp_x.nodeIdx[2] = kk * nx*ny + (jj + 1) * nx + ii;// local node 2
				temp_x.nodeIdx[3] = (kk + 1) * nx*ny + jj * nx + ii; //local node 4
				temp_x.nodeIdx[4] = kk * nx*ny + (jj + 1) * nx + ii + 1; //local node 3
				temp_x.nodeIdx[5] = (kk + 1) * nx*ny + (jj + 1) * nx + ii; //local node 6
				temp_x.nodeIdx[6] = (kk + 1) * nx*ny + jj * nx + ii + 1; //local node 5
				temp_x.nodeIdx[7] = (kk + 1) * nx*ny + (jj + 1) * nx + ii + 1; //local node 

				x.push_back(temp_x);
			}
		}
	}

}