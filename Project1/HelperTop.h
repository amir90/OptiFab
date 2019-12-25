#pragma once
#include <Eigen\Core>
#include <Eigen\Dense>
#include <Eigen/Sparse>

struct designVariable {
	int nodeIdx[8]; //each cell has eight nodes
	int varIdx; //index of design variable
	double value = 0.5; // x value
	int neighbors[6]; //idx of neighbors
	std::vector<double> position;
	std::vector<int> rho_voxInd; //voxels which influence this voxel's rho
	bool tunable = true; //integrate with voxel tunable parameter?
	double rho = 0.5;//weighted value of element
	Eigen::VectorXd displacements = Eigen::VectorXd::Zero(6);
	//double displacements[24];
	std::vector<int> influencesVoxels; //voxels which this voxel influences
	double sum_wj = 0;
};

double Min(double d1, double d2);
double Max(double d1, double d2);
int Min(int d1, int d2);
int Max(int d1, int d2);
double length(std::vector<double> a);
int moveX(std::vector<designVariable> const &x, int curridx, int move);
int moveY(std::vector<designVariable> const &x, int curridx, int move);
int moveZ(std::vector<designVariable> const &x, int curridx, int move);
int turn_untunable(std::vector<std::vector<double>> const &nodes, std::vector<designVariable>  &x, int NeighbourLayers, double dx, std::vector<std::pair<int, std::vector<double>>>& force);
void makeSurfaceMesh(std::vector<designVariable> const &x, std::vector<std::vector<double>> const &nodes);
void makeSurfaceMesh2(std::vector<designVariable> const &x, std::vector<std::vector<double>> const &nodes);
void makeSurfaceMesh3(std::vector<designVariable> const &x, std::vector<std::vector<double>> const &nodes);
void makeMesh(std::vector<designVariable> const &x, std::vector<std::vector<double>> const &nodes);
void calc_rho(std::vector<std::vector<double>> const &nodes, std::vector<designVariable>  &x, int NeighbourLayers, double dx, int nVx, int nVy, int nVz, float r0);
void create_mesh(std::vector<std::vector<double>>& nodes, std::vector<designVariable> &x, int numOfVoxelsX, int numOfVoxelsY, int numOfVoxelsZ, float dx);
